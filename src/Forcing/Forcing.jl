module Forcing

using Oceananigans
using Oceananigans.Units: Time
using Oceananigans.OutputReaders: Cyclical, AbstractInMemoryBackend, FlavorOfFTS, time_indices
using Oceananigans.Operators: Ax, Ay, volume
using Oceananigans: fill_halo_regions!, nodes, interior
using Oceananigans.Grids: x_domain, y_domain
using Dates: DateTime, Day
import Dates
using Adapt
using NCDatasets

import Oceananigans: on_architecture
import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend

using ..Configs:
    AbstractForcingConfig,
    FjordConfig,
    bathymetry_path,
    forcing_path,
    forcing_directory,
    coverage_window,
    domain_grid,
    open_edges,
    simulation_grid
using ..Plotting: plot_forcing
using ..Regridding:
    Regridding,
    native_times_to_seconds,
    FORCING_DEFLATE_LEVEL,
    interpolation_architecture,
    ProjectedSourceGrid,
    SourceRecord,
    PreparedVariable,
    uniform_time_steps,
    pad_time_steps,
    water_mask,
    projected_target_nodes,
    source_validity,
    source_fill,
    source_level_coordinate,
    source_field_grid,
    set_source_field!,
    interpolate_to_target!,
    SourceReader,
    blended_slab,
    define_forcing_dimensions!

export forcing_from_file,
    simulation_forcing,
    forcing_date_range,
    prepare_forcing,
    download_forcing,
    forcing_time_steps,
    forcing_source_grid,
    forcing_variable_names,
    forcing_monthly_filename,
    NorKystConfig,
    NorKystHindcastConfig,
    # NorKyst subsetting, shared with the NorKyst adapters in `InitialState` and `Boundaries`
    NORKYST_PROJECTION_VARIABLE,
    NorKystSubset,
    NORKYST_VARIABLE_NAMES,
    NORKYST_CATALOG_URL,
    NORKYST_OPENDAP_URL,
    norkyst_source_grid,
    subset_ranges,
    list_opendap_files,
    list_hindcast_files,
    bounding_range,
    time_length,
    define_output_file,
    write_parameter_chunk!,
    write_time_dependent_coordinates!,
    decoded_attributes,
    # the prepared file's staggering, shared with the initial state's
    data_location,
    forcing_dimension_names

""" Custom backend for FieldTimeSeries """
struct NetCDFBackend <: AbstractInMemoryBackend{Int}
    start::Int
    length::Int
end

NetCDFBackend(length) = NetCDFBackend(1, length)
new_backend(::NetCDFBackend, start, length) = NetCDFBackend(start, length)

Base.length(backend::NetCDFBackend) = backend.length
Base.summary(backend::NetCDFBackend) = string("NetCDFBackend(", backend.start, ", ", backend.length, ")")

const NetCDFFTS = FlavorOfFTS{<:Any,<:Any,<:Any,<:Any,<:NetCDFBackend}

function data_location(variable::Symbol)
    if variable in (:u,)
        return (Face, Center, Center)
    elseif variable in (:v,)
        return (Center, Face, Center)
    else
        return (Center, Center, Center)
    end
end

# Runtime generation of dictionaries based on forcing_variables_names
function oceananigans_field_names(variables)
    Dict(Symbol(variable) => Val{Symbol(variable)}() for variable in variables)
end

function data_locations(variables)
    Dict(Val{Symbol(variable)}() => data_location(Symbol(variable)) for variable in variables)
end


""" A custom forcing callable structure """
struct ForcingFromFile{FTS,V}
    fts_value::FTS
    fts_λ::FTS
    field_name::V
end

Adapt.adapt_structure(to, p::ForcingFromFile) =
    ForcingFromFile(Adapt.adapt(to, p.fts_value), Adapt.adapt(to, p.fts_λ), Adapt.adapt(to, p.field_name))

on_architecture(to, forcing::ForcingFromFile) = ForcingFromFile(
    on_architecture(to, forcing.fts_value),
    on_architecture(to, forcing.fts_λ),
    on_architecture(to, forcing.field_name),
)

"""
    forcing_term_x_flux(λ, flux, i, j, k, grid, args...)

Convert a face-crossing flux density `flux` (units `[C]·length/time`, e.g. a tracer flux in
`units(C)·m/s`) entering cell `(i, j, k)` through its x-face into a volumetric tendency for
whatever field the caller applies it to, mirroring one face-term of Oceananigans' own divergence
operator `∇·f = (1/V)[δx(Ax·fx) + ...]` — but one-sided (a source, not a two-face divergence).

Only dimensionally valid for a tracer-like field: a velocity tendency needs acceleration units,
which nothing here supplies. `ForcingFromFile` is generic in `field_name`, so this is reachable for
any field whose `λ` happens to fall in `λ > 1` — it is the caller's responsibility to only route a
field here when `flux` really is a flux density for that field. There is currently no producer of
this `λ` regime at all — `prepared_variable` writes zeros, and the rivers enter through a file of
their own (`river_forcing`) — so this function is unreachable in practice; it is kept, named
honestly, for the day a producer of `λ > 1` exists.
"""
function forcing_term_x_flux(λ, flux, i, j, k, grid, args...)
    return flux * Ax(i, j, k, grid, Center(), Center(), Center()) / volume(i, j, k, grid, Center(), Center(), Center())
end

"""
    forcing_term_y_flux(λ, flux, i, j, k, grid, args...)

As `forcing_term_x_flux`, for a flux entering through the y-face (`Ay` in place of `Ax`). Reached
when `λ < -1`.
"""
function forcing_term_y_flux(λ, flux, i, j, k, grid, args...)
    return flux * Ay(i, j, k, grid, Center(), Center(), Center()) / volume(i, j, k, grid, Center(), Center(), Center())
end

""" relaxation term """
function forcing_term_relax(λ, value, i, j, k, grid, field)
    return -λ * (field - value)
end

""" 
Return a result to be added to the tendency contributions (a kernel function).
"""
@inline function (p::ForcingFromFile{FTS,V})(i, j, k, grid, clock, fields) where {FTS,V}
    value = @inbounds p.fts_value[i, j, k, Time(clock.time)]
    λ = @inbounds p.fts_λ[i, j, k, Time(clock.time)]
    result = 0.0
    result += @inbounds ifelse(
        λ > 1 && value > -990,
        forcing_term_x_flux(λ, value, i, j, k, grid, fields[i, j, k, p.field_name]),
        0,
    )
    result += @inbounds ifelse(
        λ < -1 && value > -990,
        forcing_term_y_flux(λ, value, i, j, k, grid, fields[i, j, k, p.field_name]),
        0,
    )
    result += @inbounds ifelse(
        -1 < λ < 1 && value > -990,
        forcing_term_relax(λ, value, i, j, k, grid, fields[i, j, k, p.field_name]),
        0,
    )
    return result
end

regularize_forcing(forcing::ForcingFromFile, field, field_name, model_field_names) = forcing

function load_from_netcdf(;
    path::String,
    variable_name::String,
    grid_size::Tuple,
    time_indices_in_memory::Tuple,
    reference_date = nothing,
)
    ds = NCDataset(path)
    variable = ds[variable_name]
    native_times = ds["time"]

    data = zeros(Float64, (grid_size[1:end]..., length(time_indices_in_memory)))
    j = 1
    for i in time_indices_in_memory
        @views data[:, :, :, j] .= coalesce.(variable[:, :, :, i], -999.0)
        j += 1
    end
    start_time = isnothing(reference_date) ? native_times[1] : reference_date
    times = convert.(Int64, native_times_to_seconds(native_times, start_time))

    close(ds)
    return data, times
end

""" Update data in the FieldTimeSeries, e.g. in ForcingFromFile forcing structures. """
function set!(fts::NetCDFFTS, path::String = fts.path, name::String = fts.name)
    ti = time_indices(fts)
    data, _ =
        load_from_netcdf(; path, variable_name = name, grid_size = size(fts)[1:end-1], time_indices_in_memory = ti)

    copyto!(interior(fts, :, :, :, :), data)
    fill_halo_regions!(fts)

    return nothing
end

"""
Return a custom ForcingFromFile forcing, which is a structure keeping 2 FieldTimeSeries
with forcing values and forcing 'lambdas'.
By default both FieldTimeSeries keep only 2 times indices in memory.
"""
function forcing_get_tuple(
    filepath,
    variable_name,
    grid,
    time_indices_in_memory,
    backend,
    field_names,
    locations,
    reference_date,
)
    field_name = field_names[Symbol(variable_name)]
    LX, LY, LZ = locations[field_name]
    grid_size_tupled = size.(nodes(grid, (LX(), LY(), LZ())))
    grid_size = Tuple(x[1] for x in grid_size_tupled)

    data, times =
        load_from_netcdf(; path = filepath, variable_name, grid_size, time_indices_in_memory, reference_date)
    dataλ, timesλ = load_from_netcdf(;
        path = filepath,
        variable_name = variable_name * "_lambda",
        grid_size,
        time_indices_in_memory,
        reference_date,
    )

    fts = FieldTimeSeries{LX,LY,LZ}(
        grid,
        times;
        backend,
        time_indexing = Cyclical(),
        path = filepath,
        name = variable_name,
    )
    copyto!(interior(fts, :, :, :, :), data)
    fill_halo_regions!(fts)

    ftsλ = FieldTimeSeries{LX,LY,LZ}(
        grid,
        timesλ;
        backend,
        time_indexing = Cyclical(),
        path = filepath,
        name = variable_name * "_lambda",
    )
    copyto!(interior(ftsλ, :, :, :, :), dataλ)
    fill_halo_regions!(ftsλ)

    _forcing = ForcingFromFile(fts, ftsλ, field_name)
    result = NamedTuple{(Symbol(variable_name),)}((_forcing,))
    return result
end

"""
Return a named tuple of forcing functions for all available variables in a NetCDF file.

`reference_date` is the instant the time axis is zeroed at, `nothing` meaning the file's own
first record. A coupled run passes its simulation config's `start_date`, so the forcing and the
atmosphere agree on what model time zero stands for.
"""
function forcing_from_file(; grid, filepath, tracers, reference_date = nothing)
    ds = NCDataset(filepath)
    grid.underlying_grid.Nx == ds.dim["Nx"] &&
        grid.underlying_grid.Ny == ds.dim["Ny"] &&
        grid.underlying_grid.Nz == ds.dim["Nz"] ||
        throw(DimensionMismatch("forcing file dimensions not equal to grid dimensions"))
    forcing_variables_names = (map(String, tracers) ∪ ("u", "v")) ∩ keys(ds)
    close(ds)

    field_names = oceananigans_field_names(forcing_variables_names)
    locations = data_locations(forcing_variables_names)

    backend = NetCDFBackend(2)
    time_indices_in_memory = (1, length(backend))
    result = mapreduce(
        variable_name -> forcing_get_tuple(
            filepath,
            variable_name,
            grid,
            time_indices_in_memory,
            backend,
            field_names,
            locations,
            reference_date,
        ),
        merge,
        forcing_variables_names,
    )

    return result
end

"""
    forcing_from_file(config::AbstractForcingConfig; grid, tracers, reference_date = nothing)

Return the forcing named tuple for the prepared forcing file this setup names, resolved by
`forcing_path`. `config` is positional so the method dispatches on it.
"""
forcing_from_file(config::AbstractForcingConfig; grid, tracers, reference_date = nothing) =
    forcing_from_file(; grid, filepath = forcing_path(config), tracers, reference_date)

"""
    simulation_forcing(config::AbstractForcingConfig, grid, filepath, tracers, reference_date)

The forcing term object `coupled_simulation` consumes, read from `filepath` — the file
`build_simulation` resolved via `simulation_forcing_path`, so this hook takes the path instead of
resolving it. The rivers are not in it: `build_simulation` adds them from their own file through
`river_forcing`.

A dispatched hook rather than a hardcoded `forcing_from_file` call, so a forcing source whose
prepared files are not the FjordSim NetCDF layout `forcing_from_file` reads could supply a
completely different reader. Every built-in source shares this one default, since that layout is a
contract fixed by the read side, not a per-source detail — see `AbstractForcingConfig`.
"""
simulation_forcing(config::AbstractForcingConfig, grid, filepath, tracers, reference_date) =
    forcing_from_file(; grid, filepath, tracers, reference_date)

"""
    simulation_forcing(::Nothing, grid, filepath, tracers, reference_date)

No forcing term at all, for a setup whose `forcing_config` is `nothing`: `coupled_simulation`
still receives a splattable `forcing` keyword, just an empty one, which is exactly
`HydrostaticFreeSurfaceModel`'s own default.
"""
simulation_forcing(::Nothing, grid, filepath, tracers, reference_date) = NamedTuple()

"""
    forcing_date_range(config::AbstractForcingConfig, filepath)

First and last date the prepared forcing file at `filepath` covers, as a `(first, last)` tuple, or
`nothing` for a source that cannot report its dates.

What `FjordSim.Simulations.validate_time_coverage` checks the run's window against. Both readers use
`Cyclical()` time indexing, which wraps rather than failing outside its data, so this is what keeps
a run that outlasts its forcing from quietly replaying the beginning.

The supertype default reads `ds["time"]` from the prepared NetCDF, which is the layout every
built-in source writes. A source that overrode `simulation_forcing` to read something else overrides
this too — otherwise its coverage would be validated by a reader it does not use. The exact mirror of
`atmosphere_date_range`, which the forcing side previously had no counterpart to: this was a bare
`forcing_date_range(filepath)` in `Simulations`, dispatched on nothing at all.
"""
forcing_date_range(config::AbstractForcingConfig, filepath) = NCDataset(filepath) do ds
    dates = ds["time"][:]
    (first(dates), last(dates))
end

"""
    forcing_date_range(::Nothing, filepath)

No dates to check, for a setup whose `forcing_config` is `nothing`.
"""
forcing_date_range(::Nothing, filepath) = nothing

# --- Forcing preparation ---

# Longest run of missing days interpolated without a separate warning, matching the default of
# `NumericalEarth.DataWrangling.fill_gaps!`.
const FORCING_MAX_GAP = 6

Regridding.interpolation_architecture(config::AbstractForcingConfig) =
    interpolation_architecture(Val(config.architecture))

# --- Extension hooks ---

"""
    forcing_time_steps(config)

Every source time record available to `config`, as `SourceRecord`s sorted by date with
duplicates dropped. `prepare_forcing` completes them to a gap-free daily axis with
`daily_time_steps`.

A new forcing dataset implements this on its `AbstractForcingConfig` subtype; see
`forcing_time_steps(config::NorKystConfig)` in `src/Forcing/norkyst.jl`.
"""
function forcing_time_steps end

"""
    forcing_source_grid(config, filepath)

Geometry of the source data in `filepath`: the coordinates it lives on and the projection
they are defined in. Return a `ProjectedSourceGrid` for any source on a regular grid in
projected meters — that reuses `source_field_grid` and `projected_target_nodes` unchanged.

A new forcing dataset implements this on its `AbstractForcingConfig` subtype.
"""
function forcing_source_grid end

"""
    forcing_variable_names(config)

Mapping from source variable name to the FjordSim forcing name it becomes, e.g.
`"temperature" => "T"`. Only the intersection with `config.parameters` is prepared, so this
also declares which variables the dataset can supply.

A new forcing dataset implements this on its `AbstractForcingConfig` subtype.
"""
function forcing_variable_names end

"""
    download_forcing(config::FjordConfig)
    download_forcing(target_grid, config::AbstractForcingConfig)

Fetch and subset the source data a later `prepare_forcing` call reads, into
`forcing_directory(config)`.

The `FjordConfig` form is the generic driver: it builds the setup's grid on the CPU and
dispatches on the forcing config, so a new dataset only implements the second method. The
grid is passed rather than the grid config because a dataset needs the domain bounds, which
`x_domain`/`y_domain` provide for any `AbstractGridConfig`.
"""
download_forcing(config::FjordConfig) =
    download_forcing(domain_grid(config.grid_config, CPU()), config.forcing_config)

download_forcing(target_grid, ::Nothing) = nothing

"""
    prepare_forcing(target_grid, config::AbstractForcingConfig)

Regrid the source files already downloaded into `forcing_directory(config)` onto `target_grid`,
and write a forcing NetCDF at `forcing_path(config)` that `forcing_from_file` reads directly.

The dataset enters only through the three hooks `forcing_variable_names`,
`forcing_time_steps` and `forcing_source_grid`; everything after them is the same for every
source.

`target_grid` should be the `ImmersedBoundaryGrid` the simulation runs on, so the output
dimensions and land mask match it exactly — the mask comes from `peripheral_node`, the same
predicate the model uses. Interpolation is one trilinear `Oceananigans.Fields.interpolate` call
per target cell, taken in the source projection and clamped at the vertical ends. Cells the model
treats as land are written as `NaN`, which `forcing_from_file` reads back as the `value > -990`
sentinel the forcing kernel skips. Days missing from the download are interpolated between their
neighbours.

`config.architecture` selects where the interpolation kernel runs, via
`interpolation_architecture`, and is independent of `target_grid`, which must stay on the CPU:
building the masks walks `peripheral_node` cell by cell, which needs scalar access to the
bathymetry. Only the target nodes, mask and output buffer are moved to the interpolation device,
and only the finished slab comes back for each write, so a GPU run keeps the same streaming
memory profile.

Relaxation lambdas are written as **zero everywhere**. The interior sponge band this file used to
carry along the open edge is gone: the open lateral boundary now nudges towards hourly exterior data at
the boundary itself (`NormalRadiation`, `GravityWaveRadiation`), and a band relaxing the same
variables a few cells inside would fight it. So the `*_lambda` variables exist only for
`forcing_from_file`'s own contract, and the prepared values are otherwise read only as initial
conditions (`FromForcing`).

# Why the interpolation is not `NumericalEarth`'s

A source grid in projected meters — NorKyst's is rotated about 59 degrees from east in the
Oslofjord region — cannot be a `LatitudeLongitudeGrid`; treating its 2D `lon`/`lat` as separable
is wrong by roughly 100 grid cells. `NumericalEarth`'s dataset path (`native_grid`,
`set!(field, metadata)`, `DatasetRestoring`) always builds a `LatitudeLongitudeGrid`, and
`Oceananigans`' field-to-field `interpolate!` cannot bridge lon/lat and projected meters, so
neither applies here. The NORA3 atmosphere avoids this only because its file has already been
regridded to a regular lon/lat grid upstream. What does work is projecting the target nodes and
interpolating against a `RectilinearGrid` in projected meters; see `source_field_grid` and
`interpolate_to_target!`.

# Coverage
`coverage` is the `(first, last)` calendar interval the run needs, from `coverage_window`. The time
axis is padded to reach both ends by replicating the nearest step; `nothing` prepares exactly the
downloaded range. See `pad_time_steps`.

# The open edges
`edges` are the lateral boundaries the domain is open on, from the setup's open-boundary data config,
or `nothing` for a setup that names none and is therefore closed all round. They reach `water_mask`,
which unmasks each of those edges' velocity face rows so the prepared values are not dropped exactly
where an open boundary reads them. A keyword rather than a field of `config`, because which edges are
open is a property of the domain and its boundary data, not of the forcing dataset.

# Returns
A named tuple with `output_file`, `times` and `variables` (the written variable names).
"""
function prepare_forcing(target_grid, config::AbstractForcingConfig; coverage = nothing, edges = nothing)
    architecture = interpolation_architecture(config)
    variable_names = forcing_variable_names(config)
    source_names = [name for name in config.parameters if haskey(variable_names, name)]
    isempty(source_names) && error(
        "None of the configured parameters $(config.parameters) map to a forcing variable. " *
        "Known $(nameof(typeof(config))) variables: $(sort(collect(keys(variable_names)))).",
    )

    steps = pad_time_steps(daily_time_steps(forcing_time_steps(config)), coverage)
    reference_file = first(steps).lower.filepath
    source = forcing_source_grid(config, reference_file)
    @info "Preparing forcing from $(length(unique(step -> step.lower.filepath, steps))) file(s), " *
          "$(length(steps)) time steps: $(first(steps).date) to $(last(steps).date)"

    variables =
        [prepared_variable(name, target_grid, source, reference_file, config, edges) for name in source_names]

    output_file = forcing_path(config)
    @info "Writing forcing file to $output_file, interpolating on $(summary(architecture))"
    write_forcing_file(output_file, target_grid, variables, steps, source, architecture)
    @info "Finished preparing forcing"

    return (; output_file, times = [step.date for step in steps], variables = [variable.name for variable in variables])
end

prepare_forcing(target_grid, ::Nothing; coverage = nothing, edges = nothing) = nothing

"""
    prepare_forcing(config::FjordConfig)

Regrid the forcing a whole setup names onto its simulation grid, and write the diagnostic plot.
Returns `nothing` when the setup names no forcing.

This is the setup-level driver, the same shape as `download_forcing(config::FjordConfig)`. The
grid comes from the processed bathymetry rather than the grid config, so the output's land mask
matches the model exactly, which means `prepare_bathymetry` must have run first. It stays on the
CPU regardless of `config.forcing_config.architecture` — that field selects where the
interpolation kernel runs, and building the masks needs scalar access to the bathymetry.

The coverage window comes from the setup's simulation config, so the prepared file spans the run
the setup describes; a setup naming no simulation config gets `nothing` and the downloaded range.
The open edges come from the setup's boundary config, so the mask leaves those edges' velocity faces
where the open boundary needs them; a setup naming no boundary data is closed all round.

# Returns
The `prepare_forcing(target_grid, config)` named tuple with `plot_file` added.
"""
function prepare_forcing(config::FjordConfig)
    isnothing(config.forcing_config) && return nothing

    bathymetry_file = bathymetry_path(config.bathymetry_config)
    isfile(bathymetry_file) || error(
        "Processed bathymetry $bathymetry_file does not exist. " *
        "Run `julia --project -m FjordSim prepare_bathymetry` for this setup first.",
    )

    grid = simulation_grid(config.grid_config, bathymetry_file, CPU())
    result = prepare_forcing(
        grid,
        config.forcing_config;
        coverage = coverage_window(config.simulation_config),
        edges = open_edges(config.boundary_config),
    )
    plot_file = plot_forcing(grid, config.forcing_config)

    @info "Prepared variables: $(join(result.variables, ", "))"
    @info "Time range: $(first(result.times)) to $(last(result.times)) ($(length(result.times)) steps)"
    @info "Forcing file saved to $(result.output_file)"
    @info "Forcing plot saved to $plot_file"

    return (; result..., plot_file)
end

"""
    daily_time_steps(records; max_gap = FORCING_MAX_GAP)

`uniform_time_steps` at a one-day spacing: the cadence of the NorKyst daily-average interior
forcing.
"""
daily_time_steps(records; max_gap = FORCING_MAX_GAP) =
    uniform_time_steps(records, Day(1); max_gap, unit = "day")

"""
    forcing_dimension_names(name)

NetCDF dimension names for a forcing variable, staggered according to its Oceananigans
location. Written in Julia (fastest-first) order, which yields the `(time, Nz, Ny, Nx)` layout
in the file.
"""
function forcing_dimension_names(name)
    LX, LY, _ = data_location(Symbol(name))
    return (LX === Face ? "Nx_faces" : "Nx", LY === Face ? "Ny_faces" : "Ny", "Nz", "time")
end

"""
    prepared_variable(source_name, target_grid, source, filepath, config, edges)

Build the mask, projected target nodes, source fill and lambdas for one source variable at its
target location. `edges` are the open lateral boundaries, or `nothing` for a closed domain.
"""
function prepared_variable(source_name, target_grid, source, filepath, config::AbstractForcingConfig, edges)
    name = forcing_variable_names(config)[source_name]
    LX, LY, LZ = data_location(Symbol(name))
    @info "Preparing target nodes and source fill for $source_name -> $name"

    mask = water_mask(target_grid, LX, LY, edges)
    longitude = Array(λnodes(target_grid, LX()))
    latitude = Array(φnodes(target_grid, LY()))
    size(mask)[1:2] == (length(longitude), length(latitude)) || error(
        "Target node count $(length(longitude))x$(length(latitude)) does not match the mask " *
        "$(size(mask)[1:2]) for $name; the grid topology must be bounded in both directions.",
    )

    x, y = projected_target_nodes(longitude, latitude, source)
    shape = (length(longitude), length(latitude))
    mask_fill = source_fill(source_validity(filepath, source_name))
    # Zero everywhere: the interior relaxation band is gone. See `prepare_forcing`.
    lambda = zeros(Float32, size(mask))

    return PreparedVariable(
        source_name,
        name,
        forcing_dimension_names(name),
        mask,
        reshape(x, shape),
        reshape(y, shape),
        source_level_coordinate(source, znodes(target_grid, LZ())),
        mask_fill,
        lambda,
    )
end

"""
    write_forcing_file(filepath, target_grid, variables, steps, source, architecture)

Write the forcing NetCDF, streaming one time step at a time so peak memory stays at a single
slab rather than a whole variable. Interpolation runs on `architecture`.
"""
function write_forcing_file(filepath, target_grid, variables, steps, source, architecture)
    isdir(dirname(filepath)) || mkpath(dirname(filepath))
    isfile(filepath) && rm(filepath; force = true)

    ds = NCDataset(filepath, "c")
    try
        define_forcing_dimensions!(ds, target_grid, [step.date for step in steps])
        for variable in variables
            # One horizontal level of one time step per chunk. Both this writer and
            # `load_from_netcdf` touch exactly one time index at a time, so a chunk spanning
            # several time steps would be decompressed and recompressed once per step; the
            # netCDF default chunking does span them and makes the write an order of magnitude
            # slower.
            chunksizes = [size(variable.mask, 1), size(variable.mask, 2), 1, 1]
            for name in (variable.name, variable.name * "_lambda")
                defVar(ds, name, Float32, variable.dimensions; chunksizes,
                    deflatelevel = FORCING_DEFLATE_LEVEL, attrib = ["_FillValue" => NaN32])
            end
        end

        # All variables of a source live on the same source grid, so one field is enough.
        source_field = Field{Center,Center,Center}(source_field_grid(source, architecture))

        # Target nodes, mask and output buffer move to the interpolation device once; only the
        # finished slab comes back per step, for the NetCDF write.
        devices = [(
            x = on_architecture(architecture, variable.x),
            y = on_architecture(architecture, variable.y),
            z = on_architecture(architecture, variable.z),
            mask = on_architecture(architecture, variable.mask),
            output = on_architecture(architecture, zeros(Float32, size(variable.mask))),
            host = Array{Float32}(undef, size(variable.mask)),
        ) for variable in variables]

        reader = SourceReader(first(steps).lower.filepath)
        announced = ""

        try
            for (index, step) in enumerate(steps)
                if step.lower.filepath != announced
                    announced = step.lower.filepath
                    @info "Regridding $(basename(announced))"
                end

                for (variable, device) in zip(variables, devices)
                    slab = blended_slab(reader, step, variable.source_name)
                    set_source_field!(source_field, slab, variable.mask_fill)
                    interpolate_to_target!(
                        device.output, source_field, device.x, device.y, device.z, device.mask, architecture,
                    )
                    copyto!(device.host, device.output)
                    ds[variable.name][:, :, :, index] = device.host
                    ds[variable.name*"_lambda"][:, :, :, index] = variable.lambda
                end
            end
        finally
            close(reader)
        end
    finally
        close(ds)
    end

    return filepath
end

include("norkyst.jl")
include("norkyst_hindcast.jl")

end # module
