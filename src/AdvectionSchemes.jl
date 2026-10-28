module AdvectionSchemes

export RiverMouthAdvection

using Adapt
using Oceananigans: Oceananigans, Center
using Oceananigans.Architectures: Architectures, on_architecture
using Oceananigans.Advection: Advection, AbstractAdvectionScheme, BoundsPreservingWENO, UpwindBiased,
    LeftBias, RightBias, _biased_interpolate_xᶠᵃᵃ, _biased_interpolate_yᵃᶠᵃ, _biased_interpolate_zᵃᵃᶠ,
    rescaled_reconstruction, upwind_biased_product
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: Field, set!
using Oceananigans.Grids: AbstractGrid
using Oceananigans.Operators: Axᶠᶜᶜ, Ayᶜᶠᶜ, Azᶜᶜᶠ
using Oceananigans.TimeSteppers: TimeSteppers

"""
    RiverMouthAdvection(; interior, width_cells)

A tracer advection config: first-order upwind for the horizontal fluxes near the river mouths,
`interior` everywhere else.

A river holds its mouth columns near fresh against salt water a cell away, and a high-order
reconstruction across that front overshoots both ways, creating water saltier or colder than any
source that the mouth keeps feeding. First-order upwind cannot overshoot, at the price of numerical
diffusion of about |u|Δx/2, which this keeps to the mouths. Every face carries one flux shared by
both cells, so the tracer stays conserved; vertical fluxes are `interior`'s everywhere.

Opt-in, and only for a setup with rivers. Applies to the tracers it is given for, so
`tracer_advection = (; T = RiverMouthAdvection(...), S = RiverMouthAdvection(...), ...)` limits it
to temperature and salinity.

# Fields
- `interior`: the scheme everywhere else, any Oceananigans tracer advection scheme, a bounded
  `WENO` included.
- `width_cells`: the zone around every river column, in cells; a horizontal face is upwinded when
  either cell beside it lies within `width_cells` of a river column in both `i` and `j`.
"""
struct RiverMouthAdvection{I}
    interior::I
    width_cells::Int
end

RiverMouthAdvection(; interior, width_cells) = RiverMouthAdvection(interior, Int(width_cells))

"""
    RiverMouthScheme(interior, upwind, zone)

The scheme a `RiverMouthAdvection` builds once the river columns are known: `upwind` on the
horizontal faces beside a cell where `zone` is one, `interior` on every other face. `zone` is a 2D
field of zeros and ones, its halos filled so a face on the domain edge reads it too.
"""
struct RiverMouthScheme{N, FT, TD, I, U, Z} <: AbstractAdvectionScheme{N, FT, TD}
    interior::I
    upwind::U
    zone::Z
    RiverMouthScheme{N, FT, TD}(interior::I, upwind::U, zone::Z) where {N, FT, TD, I, U, Z} =
        new{N, FT, TD, I, U, Z}(interior, upwind, zone)
end

function RiverMouthScheme(interior, upwind, zone)
    N = max(Oceananigans.Grids.required_halo_size_x(interior), Oceananigans.Grids.required_halo_size_y(interior),
            Oceananigans.Grids.required_halo_size_z(interior))
    TD = typeof(TimeSteppers.time_discretization(interior))
    return RiverMouthScheme{N, eltype(interior), TD}(interior, upwind, zone)
end

"""
    model_tracer_advection(advection, grid, river_columns)

The tracer advection `coupled_model` gives `HydrostaticFreeSurfaceModel`.

The fallback is the identity, so any Oceananigans scheme a setup writes passes through unchanged, and
a `NamedTuple` of per-tracer schemes is resolved entry by entry. A `RiverMouthAdvection` builds its
zone from `river_columns`, every column a river enters through as `(; target_i, target_j)`.
"""
model_tracer_advection(advection, grid, river_columns) = advection

model_tracer_advection(advection::NamedTuple, grid, river_columns) =
    map(scheme -> model_tracer_advection(scheme, grid, river_columns), advection)

function model_tracer_advection(config::RiverMouthAdvection, grid, river_columns)
    isnothing(river_columns) && error(
        "RiverMouthAdvection upwinds around the columns the rivers enter through, but this setup " *
        "names no rivers. Name a `river_config` on the `FjordConfig`, or drop RiverMouthAdvection.",
    )
    upwind = UpwindBiased(eltype(grid); order = 1)
    return RiverMouthScheme(config.interior, upwind, river_mouth_zone(grid, river_columns, config.width_cells))
end

"""
    river_mouth_zone(grid, river_columns, width_cells)

A 2D field on `grid`, one within `width_cells` of a river column in both `i` and `j`, zero elsewhere.
"""
function river_mouth_zone(grid, river_columns, width_cells)
    Nx, Ny = size(grid, 1), size(grid, 2)
    mask = zeros(eltype(grid), Nx, Ny, 1)
    for (i, j) in zip(river_columns.target_i, river_columns.target_j)
        mask[max(1, i - width_cells):min(Nx, i + width_cells), max(1, j - width_cells):min(Ny, j + width_cells), 1] .= 1
    end
    zone = Field{Center, Center, Nothing}(grid)
    set!(zone, mask)
    fill_halo_regions!(zone)
    return zone
end

#####
##### The Oceananigans interface
#####

Base.summary(scheme::RiverMouthScheme) =
    string("RiverMouthScheme(interior=", summary(scheme.interior), ", near the mouths=", summary(scheme.upwind), ")")

Base.show(io::IO, scheme::RiverMouthScheme) = print(io, summary(scheme))

TimeSteppers.time_discretization(scheme::RiverMouthScheme) = TimeSteppers.time_discretization(scheme.interior)

Adapt.adapt_structure(to, scheme::RiverMouthScheme{N, FT, TD}) where {N, FT, TD} =
    RiverMouthScheme{N, FT, TD}(Adapt.adapt(to, scheme.interior), Adapt.adapt(to, scheme.upwind),
                                Adapt.adapt(to, scheme.zone))

Architectures.on_architecture(arch, scheme::RiverMouthScheme{N, FT, TD}) where {N, FT, TD} =
    RiverMouthScheme{N, FT, TD}(on_architecture(arch, scheme.interior), on_architecture(arch, scheme.upwind),
                                on_architecture(arch, scheme.zone))

# The order is adapted and the limiter field allocated on `interior`, which is what reads them.
Advection.adapt_advection_order(scheme::RiverMouthScheme, grid::AbstractGrid) =
    RiverMouthScheme(Advection.adapt_advection_order(scheme.interior, grid), scheme.upwind, scheme.zone)

Advection.materialize_advection(scheme::RiverMouthScheme, grid) =
    RiverMouthScheme(Advection.materialize_advection(scheme.interior, grid),
                     Advection.materialize_advection(scheme.upwind, grid), scheme.zone)

Advection.update_advection!(scheme::RiverMouthScheme, model, tracer) =
    Advection.update_advection!(scheme.interior, model, tracer)

#####
##### Fluxes
#####

@inline in_zone_x(i, j, zone) = @inbounds (zone[i - 1, j, 1] + zone[i, j, 1]) > 0
@inline in_zone_y(i, j, zone) = @inbounds (zone[i, j - 1, 1] + zone[i, j, 1]) > 0

@inline Advection._advective_tracer_flux_x(i, j, k, grid, scheme::RiverMouthScheme, U, c) =
    ifelse(in_zone_x(i, j, scheme.zone),
           Advection._advective_tracer_flux_x(i, j, k, grid, scheme.upwind, U, c),
           interior_flux_x(i, j, k, grid, scheme.interior, U, c))

@inline Advection._advective_tracer_flux_y(i, j, k, grid, scheme::RiverMouthScheme, V, c) =
    ifelse(in_zone_y(i, j, scheme.zone),
           Advection._advective_tracer_flux_y(i, j, k, grid, scheme.upwind, V, c),
           interior_flux_y(i, j, k, grid, scheme.interior, V, c))

@inline Advection._advective_tracer_flux_z(i, j, k, grid, scheme::RiverMouthScheme, W, c) =
    interior_flux_z(i, j, k, grid, scheme.interior, W, c)

@inline interior_flux_x(i, j, k, grid, scheme, U, c) =
    Advection._advective_tracer_flux_x(i, j, k, grid, scheme, U, c)
@inline interior_flux_y(i, j, k, grid, scheme, V, c) =
    Advection._advective_tracer_flux_y(i, j, k, grid, scheme, V, c)
@inline interior_flux_z(i, j, k, grid, scheme, W, c) =
    Advection._advective_tracer_flux_z(i, j, k, grid, scheme, W, c)

# A bounded WENO limits its reconstructions only inside its own `div_Uc`, so wrapping it would drop
# the bounds. These are the same face fluxes, one face at a time: each side's reconstruction rescaled
# by the limiter of the cell it was made in.
@inline function interior_flux_x(i, j, k, grid, scheme::BoundsPreservingWENO, U, c)
    θ = scheme.bounds.limiter
    ĉᴸ = _biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, LeftBias, c)
    ĉᴿ = _biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, RightBias, c)
    cᴸ = rescaled_reconstruction(ĉᴸ, i - 1, j, k, grid, θ, c)
    cᴿ = rescaled_reconstruction(ĉᴿ, i, j, k, grid, θ, c)
    return Axᶠᶜᶜ(i, j, k, grid) * upwind_biased_product(@inbounds(U[i, j, k]), cᴸ, cᴿ)
end

@inline function interior_flux_y(i, j, k, grid, scheme::BoundsPreservingWENO, V, c)
    θ = scheme.bounds.limiter
    ĉᴸ = _biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, LeftBias, c)
    ĉᴿ = _biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, RightBias, c)
    cᴸ = rescaled_reconstruction(ĉᴸ, i, j - 1, k, grid, θ, c)
    cᴿ = rescaled_reconstruction(ĉᴿ, i, j, k, grid, θ, c)
    return Ayᶜᶠᶜ(i, j, k, grid) * upwind_biased_product(@inbounds(V[i, j, k]), cᴸ, cᴿ)
end

@inline function interior_flux_z(i, j, k, grid, scheme::BoundsPreservingWENO, W, c)
    θ = scheme.bounds.limiter
    ĉᴸ = _biased_interpolate_zᵃᵃᶠ(i, j, k, grid, scheme, LeftBias, c)
    ĉᴿ = _biased_interpolate_zᵃᵃᶠ(i, j, k, grid, scheme, RightBias, c)
    cᴸ = rescaled_reconstruction(ĉᴸ, i, j, k - 1, grid, θ, c)
    cᴿ = rescaled_reconstruction(ĉᴿ, i, j, k, grid, θ, c)
    return Azᶜᶜᶠ(i, j, k, grid) * upwind_biased_product(@inbounds(W[i, j, k]), cᴸ, cᴿ)
end

end # module
