module BoundaryConditions

export air_sea_flux_boundary_conditions,
    quadratic_bottom_drag_boundary_conditions,
    top_bottom_boundary_conditions,
    boundary_condition_sides,
    field_boundary_conditions,
    AirSeaFluxes,
    QuadraticBottomDrag,
    OpenLateralBoundaryFromData,
    MergedBoundaryConditions

using Oceananigans
using Oceananigans.BoundaryConditions: FluxBoundaryCondition
using Oceananigans.Fields: ZeroField
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryCondition
using NumericalEarth.Oceans:
    u_quadratic_bottom_drag,
    v_quadratic_bottom_drag,
    u_immersed_bottom_drag,
    v_immersed_bottom_drag,
    build_tracer_top_bc
using ..Configs: AbstractBoundaryConditionConfig, AbstractBoundaryConditionSetConfig
using ..Utils: recursive_merge

"""
    air_sea_flux_boundary_conditions(grid)

The air-sea exchange surface: wind stress on `u` and `v`, and heat and salt flux on `T` and `S`,
each as a fresh flux field the coupled model writes into.

Only `T` and `S` get a top condition, deliberately: NumericalEarth assembles air-sea fluxes for heat
and salt specifically, so a third tracer has no exchange to write into and no `build_tracer_top_bc`
method to build one with.
"""
function air_sea_flux_boundary_conditions(grid)
    top_zonal_momentum_flux = τx = Field{Face,Center,Nothing}(grid)
    top_meridional_momentum_flux = τy = Field{Center,Face,Nothing}(grid)
    top_ocean_heat_flux = Jᵀ = Field{Center,Center,Nothing}(grid)
    top_salt_flux = Jˢ = Field{Center,Center,Nothing}(grid)

    # The tracer top conditions wrap their flux field in a `FreshwaterExchange` because
    # NumericalEarth's `net_fluxes` reads the exchange back out of the boundary condition to get the
    # freshwater volume flux and its heat content. A bare `FluxBoundaryCondition(Jᵀ)` has nothing to
    # read, so the coupled model raises a `MethodError` the first time it assembles fluxes.
    #
    # `Jʷ` stays zero here: NumericalEarth turns freshwater into a volume change through a
    # free-surface forcing it adds only on a mutable-z grid, and this grid's z faces are static.
    # So the surface T and S fluxes are still exactly `Jᵀ` and `Jˢ`.
    top_freshwater_volume_flux = Jʷ = Field{Center,Center,Nothing}(grid)
    freshwater_heat_content = Field{Center,Center,Nothing}(grid)
    freshwater_salt_content = ZeroField()

    return (
        u = (top = FluxBoundaryCondition(τx),),
        v = (top = FluxBoundaryCondition(τy),),
        T = (top = build_tracer_top_bc(Jᵀ, Jʷ, freshwater_heat_content, nothing, :T),),
        S = (top = build_tracer_top_bc(Jˢ, Jʷ, freshwater_salt_content, nothing, :S),),
    )
end

"""
    quadratic_bottom_drag_boundary_conditions(coefficient)

Quadratic drag on `u` and `v` at the bottom, at the given dimensionless drag `coefficient`.

Needs neither the grid nor the forcing — all four of `u_quadratic_bottom_drag`,
`v_quadratic_bottom_drag`, `u_immersed_bottom_drag` and `v_immersed_bottom_drag` are discrete-form
functions taking the coefficient as their parameter — which is exactly why it is separable from the
surface half. NumericalEarth's drag also adds an unresolved background speed `ub` in quadrature to the
resolved one; it is zero here, which is the plain `-μ |u| u` drag.

**Both halves are required, and the `bottom` one alone does nothing on a fjord.** A `bottom` condition
acts at the *underlying* grid's floor, and `u_quadratic_bottom_drag` hardcodes that level's index,
reading `Φ.u[i, j, 1]`. On an `ImmersedBoundaryGrid` whose deepest column stops short of `k = 1` —
which every one of these setups is, `z_faces` reaching well below the deepest sounding so no column is
clipped — that cell is immersed in the entire domain and the flux is discarded. The seabed the fjord
actually has is the *immersed* boundary, and it takes an `ImmersedBoundaryCondition` to put a stress
on it; without one, Oceananigans' default there is free slip and the water column slides over the bed.

Only `bottom` inside the `ImmersedBoundaryCondition`, matching what NumericalEarth's own
`ocean_simulation` builds: immersed side walls stay free-slip, which is a separate choice.
"""
function quadratic_bottom_drag_boundary_conditions(coefficient)
    parameters = (μ = coefficient, ub = zero(coefficient))
    return (
        u = (
            bottom = FluxBoundaryCondition(u_quadratic_bottom_drag; discrete_form=true, parameters),
            immersed = ImmersedBoundaryCondition(
                bottom = FluxBoundaryCondition(u_immersed_bottom_drag; discrete_form=true, parameters),
            ),
        ),
        v = (
            bottom = FluxBoundaryCondition(v_quadratic_bottom_drag; discrete_form=true, parameters),
            immersed = ImmersedBoundaryCondition(
                bottom = FluxBoundaryCondition(v_immersed_bottom_drag; discrete_form=true, parameters),
            ),
        ),
    )
end

"""
    top_bottom_boundary_conditions(; grid, bottom_drag_coefficient)

The air-sea surface and the drag bottom together, as one nested named tuple.

Kept as the merge of `air_sea_flux_boundary_conditions` and
`quadratic_bottom_drag_boundary_conditions` rather than as their common ancestor: the two halves
share no computation, and a setup usually wants them as two independent configs it can swap one at a
time. This remains the one place the whole `(u, v, T, S)` top-and-bottom shape is stated in a single
expression, which is what the tests assert against.
"""
top_bottom_boundary_conditions(; grid, bottom_drag_coefficient) = recursive_merge(
    air_sea_flux_boundary_conditions(grid),
    quadratic_bottom_drag_boundary_conditions(bottom_drag_coefficient),
)

"""
    boundary_condition_sides(config, grid, forcing, boundary_config, tracers, boundaries)
    boundary_condition_sides(config, grid, forcing, boundary_config, tracers)

One `AbstractBoundaryConditionConfig`'s contribution, as a nested named tuple keyed by field and
then by side: `(u = (top = …, bottom = …), T = (top = …,))`.

The per-piece hook. Named for what it returns — the per-field *sides* that
`field_boundary_conditions` materializes — rather than `boundary_conditions`, which is what it used
to be called. That name shadowed `Oceananigans.Fields.boundary_conditions`, so it was the one
FjordSim hook that could not be re-exported and had to be reached as
`FjordSim.BoundaryConditions.boundary_conditions`; and any caller binding a local of that name
shadowed the hook it was trying to call.

`boundary_config` is the setup's `AbstractBoundaryDataConfig`, or `nothing` — what a piece reads
`open_edges` from. It was the *forcing* config while the boundary dataset and the open edge both hung
off that one; the arity is unchanged, so an out-of-tree piece that ignores this argument, as both
built-in pieces below do, is unaffected. One that reads it now receives the boundary config.

`boundaries` is the prepared open-boundary state — the `NamedTuple` of reduced `FieldTimeSeries`
`FjordSim.BoundaryData.boundary_series` returns, or `nothing` for a setup with no boundary data. Passed
in rather than read here for the same reason `forcing` is: `build_simulation` is where the
`start_date` every time axis is zeroed at is known, and a config that read the file itself would
have to be told that instant a second time.

The five-argument form is the fallback, and forwards to the six-argument one *without* `boundaries`,
so a boundary config that ignores the prepared state — the built-in `AirSeaFluxes` and
`QuadraticBottomDrag`, or an out-of-tree config with a constant exterior — implements the shorter
signature and never sees it.

Declared with no six-argument method beyond that forwarding, so a subtype that implements neither
fails as a `MethodError` naming the hook.
"""
boundary_condition_sides(config, grid, forcing, boundary_config, tracers, boundaries) =
    boundary_condition_sides(config, grid, forcing, boundary_config, tracers)

"""
    AirSeaFluxes()

The air-sea exchange surface: wind stress on `u` and `v`, heat flux on `T` and salt flux on `S`.

Carries no fields — the flux fields are allocated per grid, and everything that writes into them is
NumericalEarth's business. See `air_sea_flux_boundary_conditions` for why the tracer top conditions
carry a `FreshwaterExchange` rather than being bare flux conditions, and why a third tracer gets no
top condition here. A biogeochemical tracer's surface condition belongs in its own
`AbstractBoundaryConditionConfig`.

Split from the bottom drag it used to be fused with, because the two are independent scientific
statements sharing no computation: a setup can change its drag law, or drop drag entirely, without
restating the surface.
"""
struct AirSeaFluxes <: AbstractBoundaryConditionConfig end

boundary_condition_sides(::AirSeaFluxes, grid, forcing, boundary_config, tracers) =
    air_sea_flux_boundary_conditions(grid)

"""
    QuadraticBottomDrag(; coefficient)

Quadratic drag on `u` and `v` at the bottom.

`coefficient` is the dimensionless drag coefficient the discrete-form drag functions take as their
parameter.

Contributes a condition on the *immersed* seabed as well as on the underlying grid's floor, which is
what makes it act at all on a fjord — see `quadratic_bottom_drag_boundary_conditions`.
"""
Base.@kwdef struct QuadraticBottomDrag <: AbstractBoundaryConditionConfig
    coefficient::Float64
end

QuadraticBottomDrag(coefficient::Real) = QuadraticBottomDrag(Float64(coefficient))

boundary_condition_sides(config::QuadraticBottomDrag, grid, forcing, boundary_config, tracers) =
    quadratic_bottom_drag_boundary_conditions(convert(eltype(grid), config.coefficient))

"""
    field_boundary_conditions(config, grid, forcing, boundary_config, tracers, boundaries)

The boundary conditions a whole `AbstractBoundaryConditionSetConfig` describes, materialized into
the `NamedTuple` a model constructor consumes.

The set-level hook, and a separate name from `boundary_condition_sides` on purpose: one function
returning a nested named tuple for one piece and a materialized one for a set of them would be two
shapes behind a single name.

Declared with no method, so a set config that does not implement it fails as a `MethodError` naming
the hook.
"""
function field_boundary_conditions end

"""
    MergedBoundaryConditions(pieces...)

The built-in `AbstractBoundaryConditionSetConfig`: merge every piece's contribution with
`recursive_merge`, then materialize each field's sides into `FieldBoundaryConditions`.

Merging is last-wins per leaf, so argument order is precedence order. Naming no pieces at all is a
valid statement — a model with default boundary conditions everywhere — and yields `(;)`.

A struct rather than the bare `Tuple` this used to dispatch on. A `Tuple` is not a type any subtype
can specialize, so the merge strategy and the `FieldBoundaryConditions` result shape were fixed for
every model FjordSim could ever assemble; and it said nothing about what the tuple *was*. A model
whose boundary objects are something else, or that wants precedence resolved another way, is now a
sibling of this type.

The materialized result may name barotropic `U`/`V` and `η` as well as the prognostic fields, which
`HydrostaticFreeSurfaceModel` accepts: it regularizes them through `assumed_field_location` and hands
`bcs.U`/`bcs.V` to the split-explicit free surface.
"""
struct MergedBoundaryConditions{P<:Tuple} <: AbstractBoundaryConditionSetConfig
    pieces::P
end

MergedBoundaryConditions(pieces::AbstractBoundaryConditionConfig...) =
    MergedBoundaryConditions(pieces)

function field_boundary_conditions(
    config::MergedBoundaryConditions,
    grid,
    forcing,
    boundary_config,
    tracers,
    boundaries = nothing,
)
    merged = mapreduce(
        piece -> boundary_condition_sides(piece, grid, forcing, boundary_config, tracers, boundaries),
        recursive_merge,
        config.pieces;
        init = (;),
    )

    return map(sides -> FieldBoundaryConditions(; sides...), merged)
end

include("open_boundary_conditions.jl")

end  # module BoundaryConditions
