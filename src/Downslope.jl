module Downslope

export DownslopeFlow

using Oceananigans: Oceananigans, Callback, Center, IterationInterval, fields
using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: interior
using Oceananigans.Grids: inactive_node, znode
using Oceananigans.Operators: Δxᶜᶠᶜ, Δxᶠᶜᶜ, Δyᶜᶠᶜ, Δyᶠᶜᶜ, Δzᶜᶜᶜ, Vᶜᶜᶜ
using Oceananigans.Utils: launch!
using KernelAbstractions: @kernel, @index, @Const
using SeawaterPolynomials: ρ′

using ..Configs: AbstractCallbackConfig, AbstractRunConfig
using ..Simulations: Simulations

"""
    DownslopeFlow(; name, tracers, friction_timescale)

A downslope-flow parameterization for z-level grids, after Campin & Goosse (1999) in the form of
MITgcm's `down_slope` package: dense water on a step goes down the slope to its own level instead
of overturning the deeper column beside it.

A z-level model carries dense water off a step sideways, onto lighter water in the next column, and
the vertical closure then mixes that column to the bottom. Wherever a column's bottom cell is denser
than the water at the same level in a deeper neighbouring column, this moves it, every time step,
along a closed loop: from that bottom cell to the deepest level of the neighbour at which it is still
denser (its neutral level), up the neighbour from there to the bottom cell's level, and back across.
The transport is MITgcm's

    U = Δy Δz s (g Δρ / ρ₀) friction_timescale

with `Δy` the face's width, `Δz` the bottom cell's thickness, `s` the bottom slope between the two
columns and `Δρ` the density excess at the bottom cell's level, both waters evaluated there with the
model's equation of state. Each step is capped so no cell of a loop exchanges more than a quarter of
its volume, since a cell can sit on up to four loops. The loop is closed and upstream, so every
tracer is conserved and stays within the values it started the step with.

A callback rather than a tracer forcing because it needs the density of one column's water at
another column's depth: it runs after every time step on the simulation's `Δt`.

# Fields
- `name`: the key it occupies in `simulation.callbacks`.
- `tracers`: the tracers moved, e.g. `(:T, :S)`. A closure's own tracer, such as CATKE's `e`, is
  left out: it is turbulence, not water.
- `friction_timescale`: `μ⁻¹` in MITgcm's notation (`DWNSLP_rec_mu`), in seconds. It sets the
  speed of the downslope flow, `s g Δρ / ρ₀ × friction_timescale`.
"""
struct DownslopeFlow{N} <: AbstractCallbackConfig
    name::Symbol
    tracers::NTuple{N,Symbol}
    friction_timescale::Float64
end

function DownslopeFlow(; name, tracers, friction_timescale)
    friction_timescale > 0 || throw(
        ArgumentError("A DownslopeFlow's `friction_timescale` must be positive, got $friction_timescale."),
    )
    isempty(tracers) && throw(ArgumentError("A DownslopeFlow moves no tracers; name at least one."))

    return DownslopeFlow(Symbol(name), Tuple(Symbol.(tracers)), Float64(friction_timescale))
end

"""
    DownslopeStep

The callable `attach_callback!` puts on the simulation: the parameters, the buffers the two passes
of `downslope_step!` share, and the model's equation of state.
"""
struct DownslopeStep{N,FT,E,U,K,I}
    tracers::NTuple{N,Symbol}
    friction_timescale::FT
    gravitational_acceleration::FT
    reference_density::FT
    equation_of_state::E
    transport_x::U
    source_level_x::K
    neutral_level_x::K
    transport_y::U
    source_level_y::K
    neutral_level_y::K
    increment::I
end

"""
    DownslopeStep(ocean_model, tracers, friction_timescale)

Allocate the face buffers on the model's device and read gravity, the reference density and the
equation of state from its `SeawaterBuoyancy`.
"""
function DownslopeStep(ocean_model, tracers, friction_timescale)
    grid = ocean_model.grid
    FT = eltype(grid)
    Nx, Ny, Nz = size(grid)
    arch = architecture(grid)

    for name in tracers
        name in keys(ocean_model.tracers) || throw(
            ArgumentError("DownslopeFlow names tracer :$name, which the model does not carry."),
        )
    end

    buoyancy = ocean_model.buoyancy.formulation
    equation_of_state = buoyancy.equation_of_state

    return DownslopeStep(
        tracers,
        convert(FT, friction_timescale),
        convert(FT, buoyancy.gravitational_acceleration),
        convert(FT, equation_of_state.reference_density),
        equation_of_state,
        on_architecture(arch, zeros(FT, Nx + 1, Ny)),
        on_architecture(arch, ones(Int32, Nx + 1, Ny)),
        on_architecture(arch, ones(Int32, Nx + 1, Ny)),
        on_architecture(arch, zeros(FT, Nx, Ny + 1)),
        on_architecture(arch, ones(Int32, Nx, Ny + 1)),
        on_architecture(arch, ones(Int32, Nx, Ny + 1)),
        on_architecture(arch, zeros(FT, Nx, Ny, Nz)),
    )
end

(step::DownslopeStep)(simulation) = downslope_step!(simulation.model.ocean.model, step, simulation.Δt)

"""
    downslope_step!(ocean_model, step, Δt)

One step of the downslope flow on `ocean_model`'s tracers: the transport across every face, then
each tracer's increment from those transports, all read from the state before the step.
"""
function downslope_step!(ocean_model, step::DownslopeStep, Δt)
    grid = ocean_model.grid
    Nx, Ny, Nz = size(grid)
    arch = architecture(grid)
    FT = eltype(grid)
    T, S = ocean_model.tracers.T, ocean_model.tracers.S
    bottom = grid.immersed_boundary.bottom_height
    parameters = (step.friction_timescale, step.gravitational_acceleration, step.reference_density, convert(FT, Δt))

    launch!(arch, grid, (Nx + 1, Ny), _downslope_faces_x!, step.transport_x, step.source_level_x,
            step.neutral_level_x, grid, T, S, bottom, step.equation_of_state, parameters)
    launch!(arch, grid, (Nx, Ny + 1), _downslope_faces_y!, step.transport_y, step.source_level_y,
            step.neutral_level_y, grid, T, S, bottom, step.equation_of_state, parameters)

    for name in step.tracers
        tracer = getproperty(ocean_model.tracers, name)
        launch!(arch, grid, :xyz, _downslope_increment!, step.increment, tracer, grid,
                step.transport_x, step.source_level_x, step.neutral_level_x,
                step.transport_y, step.source_level_y, step.neutral_level_y, convert(FT, Δt))
        interior(tracer) .+= step.increment
        # With the clock and fields, as the model fills them: an open boundary's value is a
        # time series it reads at the model's time.
        fill_halo_regions!(tracer, ocean_model.clock, fields(ocean_model))
    end

    return nothing
end

"""
    bottom_level(i, j, grid)

The lowest active level of column `(i, j)`, or `Nz + 1` for a land column.
"""
@inline function bottom_level(i, j, grid)
    Nz = size(grid, 3)
    level = Nz + 1
    for k in Nz:-1:1
        level = ifelse(inactive_node(i, j, k, grid, Center(), Center(), Center()), level, k)
    end
    return level
end

"""
    face_transport(i₁, j₁, i₂, j₂, width, distance, grid, T, S, bottom, eos, parameters)

The downslope transport across the face between columns 1 and 2, in m³/s, positive when column 1 is
the shallow one; the shallow column's bottom level; and the deep column's neutral level for that
water. Zero where the columns share a bottom level, the shallow water is not denser, or it is not
denser one level down.
"""
@inline function face_transport(i₁, j₁, i₂, j₂, width, distance, grid, T, S, bottom, eos, parameters)
    friction_timescale, g, ρ₀, Δt = parameters
    Nz = size(grid, 3)
    c = Center()

    bottom₁, bottom₂ = bottom_level(i₁, j₁, grid), bottom_level(i₂, j₂, grid)
    first_shallow = bottom₁ > bottom₂
    iₛ, jₛ = ifelse(first_shallow, i₁, i₂), ifelse(first_shallow, j₁, j₂)
    iₙ, jₙ = ifelse(first_shallow, i₂, i₁), ifelse(first_shallow, j₂, j₁)
    deep_bottom = min(bottom₁, bottom₂)
    kₛ = min(max(bottom₁, bottom₂), Nz)

    Tₛ, Sₛ = T[iₛ, jₛ, kₛ], S[iₛ, jₛ, kₛ]
    zₛ = znode(iₛ, jₛ, kₛ, grid, c, c, c)
    Δρ = ρ′(Tₛ, Sₛ, zₛ, eos) - ρ′(T[iₙ, jₙ, kₛ], S[iₙ, jₙ, kₛ], zₛ, eos)

    neutral = kₛ
    sinking = true
    for k in (kₛ - 1):-1:1
        z = znode(iₙ, jₙ, k, grid, c, c, c)
        sinking = sinking & (k >= deep_bottom) & (ρ′(Tₛ, Sₛ, z, eos) > ρ′(T[iₙ, jₙ, k], S[iₙ, jₙ, k], z, eos))
        neutral = ifelse(sinking, k, neutral)
    end

    smallest_volume = Vᶜᶜᶜ(iₛ, jₛ, kₛ, grid)
    for k in 1:Nz
        on_loop = (k >= neutral) & (k <= kₛ)
        smallest_volume = ifelse(on_loop, min(smallest_volume, Vᶜᶜᶜ(iₙ, jₙ, k, grid)), smallest_volume)
    end

    slope = (bottom[iₛ, jₛ, 1] - bottom[iₙ, jₙ, 1]) / distance
    transport = width * Δzᶜᶜᶜ(iₛ, jₛ, kₛ, grid) * slope * g * Δρ / ρ₀ * friction_timescale
    transport = min(transport, smallest_volume / (4 * Δt))

    active = (bottom₁ != bottom₂) & (max(bottom₁, bottom₂) <= Nz) & (Δρ > 0) & (neutral < kₛ)
    transport = ifelse(active, transport, zero(transport))

    return ifelse(first_shallow, transport, -transport), kₛ, neutral
end

@kernel function _downslope_faces_x!(transport, source_level, neutral_level, grid, @Const(T), @Const(S),
                                     @Const(bottom), eos, parameters)
    i, j = @index(Global, NTuple)
    Nx = size(grid, 1)
    i₁ = max(i - 1, 1)
    i₂ = min(i, Nx)
    U, kₛ, kₙ = face_transport(i₁, j, i₂, j, Δyᶠᶜᶜ(i₂, j, 1, grid), Δxᶠᶜᶜ(i₂, j, 1, grid),
                               grid, T, S, bottom, eos, parameters)
    on_boundary = (i == 1) | (i == Nx + 1)
    @inbounds transport[i, j] = ifelse(on_boundary, zero(U), U)
    @inbounds source_level[i, j] = kₛ
    @inbounds neutral_level[i, j] = kₙ
end

@kernel function _downslope_faces_y!(transport, source_level, neutral_level, grid, @Const(T), @Const(S),
                                     @Const(bottom), eos, parameters)
    i, j = @index(Global, NTuple)
    Ny = size(grid, 2)
    j₁ = max(j - 1, 1)
    j₂ = min(j, Ny)
    U, kₛ, kₙ = face_transport(i, j₁, i, j₂, Δxᶜᶠᶜ(i, j₂, 1, grid), Δyᶜᶠᶜ(i, j₂, 1, grid),
                               grid, T, S, bottom, eos, parameters)
    on_boundary = (j == 1) | (j == Ny + 1)
    @inbounds transport[i, j] = ifelse(on_boundary, zero(U), U)
    @inbounds source_level[i, j] = kₛ
    @inbounds neutral_level[i, j] = kₙ
end

"""
    loop_flux(U, kₛ, kₙ, this_first, C, i, j, k, iₒ, jₒ)

What one face's loop brings into cell `(i, j, k)`, in tracer units × m³/s. `this_first` says
whether this cell's column is the face's column 1, `(iₒ, jₒ)` is the other column. As the shallow
column's bottom cell, the cell trades its water for the deep column's at the same level; as one of
the deep column's cells from the neutral level up, it takes the water from below it, or at the
neutral level the shallow column's.
"""
@inline function loop_flux(U, kₛ, kₙ, this_first, C, i, j, k, iₒ, jₒ)
    magnitude = abs(U)
    this_shallow = this_first == (U > 0)
    from_shallow = C[iₒ, jₒ, kₛ]
    from_below = C[i, j, max(k - 1, 1)]
    as_shallow = ifelse(k == kₛ, magnitude * (from_shallow - C[i, j, k]), zero(magnitude))
    on_loop = (k >= kₙ) & (k <= kₛ)
    upstream = ifelse(k == kₙ, from_shallow, from_below)
    as_deep = ifelse(on_loop, magnitude * (upstream - C[i, j, k]), zero(magnitude))
    return ifelse(this_shallow, as_shallow, as_deep)
end

@kernel function _downslope_increment!(increment, @Const(C), grid, @Const(transport_x), @Const(source_x),
                                       @Const(neutral_x), @Const(transport_y), @Const(source_y),
                                       @Const(neutral_y), Δt)
    i, j, k = @index(Global, NTuple)
    Nx, Ny, _ = size(grid)
    west = loop_flux(transport_x[i, j], source_x[i, j], neutral_x[i, j], false, C, i, j, k, max(i - 1, 1), j)
    east = loop_flux(transport_x[i + 1, j], source_x[i + 1, j], neutral_x[i + 1, j], true, C, i, j, k, min(i + 1, Nx), j)
    south = loop_flux(transport_y[i, j], source_y[i, j], neutral_y[i, j], false, C, i, j, k, i, max(j - 1, 1))
    north = loop_flux(transport_y[i, j + 1], source_y[i, j + 1], neutral_y[i, j + 1], true, C, i, j, k, i, min(j + 1, Ny))
    flux = west + east + south + north
    dry = inactive_node(i, j, k, grid, Center(), Center(), Center())
    @inbounds increment[i, j, k] = ifelse(dry, zero(flux), Δt * flux / Vᶜᶜᶜ(i, j, k, grid))
end

"""
    attach_callback!(simulation, callback::DownslopeFlow, config)

Attach the downslope flow to the coupled simulation under `callback.name`, every iteration.
"""
function Simulations.attach_callback!(simulation, callback::DownslopeFlow, config::AbstractRunConfig)
    step = DownslopeStep(simulation.model.ocean.model, callback.tracers, callback.friction_timescale)
    simulation.callbacks[callback.name] = Callback(step, IterationInterval(1))
    return simulation
end

end
