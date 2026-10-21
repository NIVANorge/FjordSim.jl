# The open lateral boundary condition: an `AbstractBoundaryConditionConfig` piece whose exterior
# state is the prepared boundary series `FjordSim.BoundaryData.boundary_series` reads.

using Oceananigans.BoundaryConditions:
    GravityWaveRadiationBoundaryCondition,
    NormalFlowBoundaryCondition,
    NormalRadiation,
    ValueBoundaryCondition
using Oceananigans.Grids: column_depthᶜᶠᵃ, column_depthᶠᶜᵃ
using Oceananigans.Units: Time

using ..Grids: open_edges, LATERAL_EDGES

#####
##### The open lateral boundary
#####
##### Three schemes, one per part of the state, matching what the parent model itself does at its
##### own lateral boundaries (NorKyst names them `zeta: Che`, `ubar/vbar: Shc`, `u/v/temp/salt:
##### RadNud` in its output attributes):
#####
##### - `GravityWaveRadiation` (Flather 1976) on the barotropic transport, given the exterior
#####   transport and elevation. This is the part that actually opens the boundary: it is what lets
#####   the tide and the subtidal barotropic exchange cross it.
##### - `SurfaceWaveRadiation` (Chapman 1985) on `η`. Not stated here at all — Oceananigans pairs it
#####   onto every side where `U` or `V` carries a gravity-wave condition
#####   (`default_free_surface_boundary_conditions`), so a second statement could only disagree.
##### - `NormalRadiation` (Orlanski 1976 with Marchesiello et al. 2001 nudging) on the 3D velocity
#####   and on the tracers: radiate outgoing signals, nudge towards the exterior state on inflow.
#####
##### Naming a `U`/`V` condition also switches the split-explicit solver to `LocalHaloFilling`
##### (`substep_halo_filling`), so the barotropic halo is filled every substep. That is what makes
##### the Flather condition act at all, and it needs no free-surface config change.

"""
    boundary_radiation_scheme(config)

The Orlanski-with-nudging scheme the 3D velocity and the tracers share, from the config's two
timescales.

One scheme object for both because the inflow/outflow distinction is the same physical statement
wherever it is applied: on outflow radiate at the locally diagnosed phase speed and barely nudge, on
inflow abandon radiation and nudge hard towards the data.
"""
boundary_radiation_scheme(config) = NormalRadiation(;
    inflow_timescale = config.inflow_timescale,
    outflow_timescale = config.outflow_timescale,
)

#####
##### Exterior barotropic transport
#####
##### `GravityWaveRadiation` wants a 2-tuple `(transport, elevation)` per boundary node. The
##### prepared boundary file holds a depth-averaged velocity, so the transport is that velocity times
##### the model's *own* column depth — the same `column_depth` operator the scheme's kernel uses, so
##### the boundary transport is expressed in the model's bathymetry rather than the source's.
#####
##### One function per edge rather than one function branching on an index: the boundary node is a
##### different index of a different operator on each side.
#####
##### No land guard: `FjordSim.BoundaryData.fill_boundary_gaps!` has already filled every dry
##### boundary cell from its nearest wet neighbour, so a boundary series is finite everywhere. That is the
##### one place the fill can live — see that function for why a sentinel reaching a scheme is not
##### hypothetical.

@inline boundary_value(series, i, k, time) = @inbounds series[i, 1, k, time]

@inline exterior_state(velocity, elevation, depth) = (velocity * depth, elevation)

"""
    south_boundary_transport(i, k, grid, clock, model_fields, parameters)

Exterior `(transport, elevation)` at the southern boundary face `(i, 1)`: `vbar` times the model's
column depth there, and `eta`, both from the prepared boundary series in `parameters`.
"""
@inline function south_boundary_transport(i, k, grid, clock, model_fields, parameters)
    time = Time(clock.time)
    depth = column_depthᶜᶠᵃ(i, 1, grid.Nz + 1, grid, model_fields.η)
    return exterior_state(
        boundary_value(parameters.vbar, i, 1, time),
        boundary_value(parameters.eta, i, 1, time),
        depth,
    )
end

"""
    north_boundary_transport(i, k, grid, clock, model_fields, parameters)

As `south_boundary_transport`, at the northern boundary face `(i, Ny + 1)`.
"""
@inline function north_boundary_transport(i, k, grid, clock, model_fields, parameters)
    time = Time(clock.time)
    depth = column_depthᶜᶠᵃ(i, grid.Ny + 1, grid.Nz + 1, grid, model_fields.η)
    return exterior_state(
        boundary_value(parameters.vbar, i, 1, time),
        boundary_value(parameters.eta, i, 1, time),
        depth,
    )
end

"""
    west_boundary_transport(j, k, grid, clock, model_fields, parameters)

As `south_boundary_transport`, at the western boundary face `(1, j)`, from `ubar`.
"""
@inline function west_boundary_transport(j, k, grid, clock, model_fields, parameters)
    time = Time(clock.time)
    depth = column_depthᶠᶜᵃ(1, j, grid.Nz + 1, grid, model_fields.η)
    return exterior_state(
        boundary_value(parameters.ubar, j, 1, time),
        boundary_value(parameters.eta, j, 1, time),
        depth,
    )
end

"""
    east_boundary_transport(j, k, grid, clock, model_fields, parameters)

As `west_boundary_transport`, at the eastern boundary face `(Nx + 1, j)`.
"""
@inline function east_boundary_transport(j, k, grid, clock, model_fields, parameters)
    time = Time(clock.time)
    depth = column_depthᶠᶜᵃ(grid.Nx + 1, j, grid.Nz + 1, grid, model_fields.η)
    return exterior_state(
        boundary_value(parameters.ubar, j, 1, time),
        boundary_value(parameters.eta, j, 1, time),
        depth,
    )
end

#####
##### The four groups, one method per edge each
#####
##### Each takes *one* edge's series — the `(; T, S, u, v, eta, ubar, vbar)` group
##### `FjordSim.BoundaryData.boundary_series` returns under that edge's key, not the whole per-edge
##### collection. `open_edge_boundary_conditions` is what unwraps a group per edge and merges the
##### four; these functions never see a second side.

"""
    boundary_series_value(boundaries, name, edge)

The prepared boundary series for `name`, or a stated error naming what to prepare. The exterior
state is data, so a missing variable is a preparation step that has not been run rather than a
default to fall back on.
"""
function boundary_series_value(boundaries, name, edge)
    haskey(boundaries, name) || error(
        "The open :$edge boundary needs `$name` but the prepared boundary file does not carry it. " *
        "Add its source variable to the boundary config's `parameters` and re-run " *
        "`prepare_boundaries`.",
    )
    return getproperty(boundaries, name)
end

"""
    open_normal_velocity_boundary_conditions(::Val{edge}, boundaries, scheme)

The 3D velocity component normal to `edge`, radiating with `scheme` towards its own prepared series.

This is the condition that used to be `NormalFlowBoundaryCondition(nothing)` — a closed wall, since
`getbc` of `nothing` is `zero(grid)`. The series is passed straight through as the condition: a
reduced `FieldTimeSeries` is a boundary condition Oceananigans indexes by the two tangential indices
on its own, so no discrete-form wrapper is needed.
"""
open_normal_velocity_boundary_conditions(::Val{:south}, boundaries, scheme) = (
    v = (south = NormalFlowBoundaryCondition(boundary_series_value(boundaries, :v, :south); scheme),),
)
open_normal_velocity_boundary_conditions(::Val{:north}, boundaries, scheme) = (
    v = (north = NormalFlowBoundaryCondition(boundary_series_value(boundaries, :v, :north); scheme),),
)
open_normal_velocity_boundary_conditions(::Val{:west}, boundaries, scheme) = (
    u = (west = NormalFlowBoundaryCondition(boundary_series_value(boundaries, :u, :west); scheme),),
)
open_normal_velocity_boundary_conditions(::Val{:east}, boundaries, scheme) = (
    u = (east = NormalFlowBoundaryCondition(boundary_series_value(boundaries, :u, :east); scheme),),
)
open_normal_velocity_boundary_conditions(::Val{edge}, boundaries, scheme) where {edge} =
    throw(ArgumentError("open_edge must be one of $LATERAL_EDGES, got :$edge"))

"""
    open_tangential_velocity_boundary_conditions(::Val{edge}, boundaries, scheme)

The 3D velocity component *tangential* to `edge`, as a `Value` condition radiating with `scheme`.

A tangential velocity is Center-located across the boundary, so it takes the same shape of condition
a tracer does: `NormalRadiation` decides inflow from the boundary-normal velocity one cell in and
writes the halo cell, never an interior one.
"""
open_tangential_velocity_boundary_conditions(::Val{:south}, boundaries, scheme) = (
    u = (south = ValueBoundaryCondition(boundary_series_value(boundaries, :u, :south); scheme),),
)
open_tangential_velocity_boundary_conditions(::Val{:north}, boundaries, scheme) = (
    u = (north = ValueBoundaryCondition(boundary_series_value(boundaries, :u, :north); scheme),),
)
open_tangential_velocity_boundary_conditions(::Val{:west}, boundaries, scheme) = (
    v = (west = ValueBoundaryCondition(boundary_series_value(boundaries, :v, :west); scheme),),
)
open_tangential_velocity_boundary_conditions(::Val{:east}, boundaries, scheme) = (
    v = (east = ValueBoundaryCondition(boundary_series_value(boundaries, :v, :east); scheme),),
)
open_tangential_velocity_boundary_conditions(::Val{edge}, boundaries, scheme) where {edge} =
    throw(ArgumentError("open_edge must be one of $LATERAL_EDGES, got :$edge"))

"""
    open_transport_boundary_conditions(::Val{edge}, boundaries)

The barotropic transport normal to `edge`, as a Flather (`GravityWaveRadiation`) condition reading
the exterior transport and elevation through one of the `*_boundary_transport` functions.

`U` for a west or east edge, `V` for a south or north one — the same staggering as the 3D velocity,
one integral up. Both prepared series ride along in `parameters`, which is also what keeps them
paged in: `Oceananigans.OutputReaders.extract_field_time_series` walks a discrete boundary
function's parameters, and reaches this one because barotropic `U` and `V` are in
`Oceananigans.fields(model)`.
"""
open_transport_boundary_conditions(::Val{:south}, boundaries) = (
    V = (
        south = GravityWaveRadiationBoundaryCondition(
            south_boundary_transport;
            discrete_form = true,
            parameters = transport_parameters(boundaries, :vbar, :south),
        ),
    ),
)
open_transport_boundary_conditions(::Val{:north}, boundaries) = (
    V = (
        north = GravityWaveRadiationBoundaryCondition(
            north_boundary_transport;
            discrete_form = true,
            parameters = transport_parameters(boundaries, :vbar, :north),
        ),
    ),
)
open_transport_boundary_conditions(::Val{:west}, boundaries) = (
    U = (
        west = GravityWaveRadiationBoundaryCondition(
            west_boundary_transport;
            discrete_form = true,
            parameters = transport_parameters(boundaries, :ubar, :west),
        ),
    ),
)
open_transport_boundary_conditions(::Val{:east}, boundaries) = (
    U = (
        east = GravityWaveRadiationBoundaryCondition(
            east_boundary_transport;
            discrete_form = true,
            parameters = transport_parameters(boundaries, :ubar, :east),
        ),
    ),
)
open_transport_boundary_conditions(::Val{edge}, boundaries) where {edge} =
    throw(ArgumentError("open_edge must be one of $LATERAL_EDGES, got :$edge"))

"""
    transport_parameters(boundaries, barotropic_name, edge)

The two prepared series a `*_boundary_transport` function reads, under the names it reads them by.

`barotropic_name` is `:ubar` or `:vbar` depending on the edge; both are bound under their own key so
one parameter shape serves all four edges without a rename.
"""
transport_parameters(boundaries, barotropic_name::Symbol, edge) = NamedTuple{(barotropic_name, :eta)}((
    boundary_series_value(boundaries, barotropic_name, edge),
    boundary_series_value(boundaries, :eta, edge),
))

"""
    open_tracer_boundary_conditions(::Val{edge}, boundaries, tracers, scheme)

One radiating `Value` condition per simulated tracer on `edge`, towards that tracer's own prepared
series.

Errors if a tracer has no prepared series: naming an open edge implies the boundary data should
supply every simulated tracer there. A biogeochemical tracer therefore needs its own source variable
in the boundary config's `parameters`.
"""
open_tracer_boundary_conditions(::Val{edge}, boundaries, tracers, scheme) where {edge} = NamedTuple(
    name => (; Symbol(edge) => ValueBoundaryCondition(boundary_series_value(boundaries, name, edge); scheme))
    for name in tracers
)

"""
    OpenLateralBoundaryFromData(; inflow_timescale, outflow_timescale)

The domain's open edge, with every exterior value read from the prepared open-boundary file.

Four groups on **each** edge the setup's `AbstractBoundaryDataConfig` names: a Flather condition on the
barotropic transport (from `ubar`/`vbar` and `eta`), and Orlanski radiation with nudging on the normal
3D velocity, on the tangential 3D velocity, and on every simulated tracer. Oceananigans adds the
Chapman condition on `η` itself, on every side carrying a gravity-wave transport condition — so a
domain in the open ocean naming all four edges needs nothing extra here or there.

The result is a boundary water, heat and salt actually cross — the
first version of this config, `OpenLateralBoundary`, put `NormalFlowBoundaryCondition(nothing)` on the
normal velocity, which is a closed wall, so only the tracers there were ever open.

The values come from the boundary data config a setup names on its `FjordConfig`, read by
`FjordSim.BoundaryData.boundary_series` and handed in by `build_simulation` — hourly NorKyst for the
setups in `examples/`, because a Flather boundary without a tide is not worth prescribing and the interior
forcing is daily means.

Called `OpenLateralBoundaryFromForcing` while the boundary dataset hung off the forcing config. The
values were never forcing — they are a separate hourly file from a separate collection read by a
separate pipeline — and the name is the last place that said otherwise.

# Fields
- `inflow_timescale`: seconds to nudge towards the exterior state on inflow. Marchesiello et al.
  (2001) use about a day.
- `outflow_timescale`: seconds to nudge on outflow, where radiation should do the work instead. About
  a year in the same reference; `Inf` is pure radiation.

Both are fields rather than derivations, unlike the single `relaxation_timescale` this used to take
from the forcing config: a genuinely open boundary treats inflow and outflow differently, and the
outflow rate is information no other config states. Nothing in this stack has a default, so a setup
states both.
"""
struct OpenLateralBoundaryFromData <: AbstractBoundaryConditionConfig
    inflow_timescale::Float64
    outflow_timescale::Float64
end

OpenLateralBoundaryFromData(; inflow_timescale, outflow_timescale) =
    OpenLateralBoundaryFromData(Float64(inflow_timescale), Float64(outflow_timescale))

function boundary_condition_sides(
    config::OpenLateralBoundaryFromData,
    grid,
    forcing,
    boundary_config,
    tracers,
    boundaries,
)
    isnothing(boundary_config) && error(
        "OpenLateralBoundaryFromData puts its schemes on the edges a boundary data config names, " *
        "but this setup's `FjordConfig` names no `boundary_config`. Give it an " *
        "`AbstractBoundaryDataConfig`, or use a boundary condition that carries its own exterior " *
        "state and edges.",
    )

    isnothing(boundaries) && error(
        "OpenLateralBoundaryFromData reads its exterior state from the prepared boundary file, but " *
        "none was read for this setup. Run `julia --project -m FjordSim prepare_boundaries` for it, " *
        "or use a boundary condition that carries its own exterior state.",
    )

    scheme = boundary_radiation_scheme(config)

    return mapreduce(
        edge -> open_edge_boundary_conditions(Val(edge), getproperty(boundaries, edge), tracers, scheme),
        recursive_merge,
        open_edges(boundary_config);
        init = (;),
    )
end

"""
    open_edge_boundary_conditions(edge::Val, series, tracers, scheme)

The four groups one open edge contributes, merged: the Flather condition on the barotropic transport,
and Orlanski radiation with nudging on the normal velocity, the tangential velocity and every tracer.

Split out so `boundary_condition_sides` can `mapreduce` it over several open edges. The merge across
edges is the same `recursive_merge` as within one, and it composes rather than colliding: on a domain
open to the south and the west, `u` picks up a *tangential* condition at `.south` and a *normal* one at
`.west`, and each tracer picks up both sides. A corner cell is the one place two edges' conditions meet,
and which of them fills it is Oceananigans' halo-fill order rather than anything stated here.
"""
open_edge_boundary_conditions(edge::Val, series, tracers, scheme) = recursive_merge(
    open_normal_velocity_boundary_conditions(edge, series, scheme),
    open_tangential_velocity_boundary_conditions(edge, series, scheme),
    open_transport_boundary_conditions(edge, series),
    open_tracer_boundary_conditions(edge, series, tracers, scheme),
)
