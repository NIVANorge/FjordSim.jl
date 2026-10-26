module InitialState

# The initial-state pipeline: one source record at the run's `start_date`, regridded onto the
# simulation grid for `FromInitialState` initial conditions. The regridding is `Regridding`'s, the
# same as the forcing's and the open boundaries'; the dataset enters only through the hooks below.

export download_initial_state,
    prepare_initial_state,
    NorKystInitialStateConfig

using Oceananigans
using NCDatasets

using ..Configs:
    AbstractInitialStateConfig,
    AbstractRunConfig,
    AbstractFjordConfig,
    grid_config,
    bathymetry_config,
    initial_state_config,
    boundary_data_config,
    run_config,
    start_date
using ..Grids: domain_grid, open_edges, simulation_grid
using ..Bathymetry: resolve_bathymetry_file
using ..Regridding:
    FORCING_DEFLATE_LEVEL,
    SourceRecord,
    SourceReader,
    water_mask,
    projected_target_nodes,
    source_validity,
    source_fill,
    source_level_coordinate,
    source_field_grid,
    set_source_field!,
    interpolate_to_target!,
    source_slab,
    define_forcing_dimensions!,
    data_location,
    forcing_dimension_names

"""
    initial_state_path(config)
    initial_state_directory(config)

Resolve `config.output_file` and `config.output_directory` against `config.data_root`: the file
`prepare_initial_state` writes, and the directory `download_initial_state` writes into. Defined for
every `AbstractInitialStateConfig`; a dataset without those fields overloads them. An absolute field
is returned unchanged.
"""
initial_state_path(config::AbstractInitialStateConfig) = joinpath(config.data_root, config.output_file)
initial_state_directory(config::AbstractInitialStateConfig) =
    joinpath(config.data_root, config.output_directory)

"""
The global attribute holding the `start_date` a prepared initial-state file was prepared for, which
`FromInitialState` checks against the run's own.
"""
const INITIAL_STATE_DATE_ATTRIBUTE = "start_date"

# --- Extension hooks ---

"""
    initial_state_record(config, date)

The downloaded `SourceRecord` holding the state at `date`. Errors, naming `download_initial_state`,
when it has not been downloaded.

A new initial-state dataset implements this on its `AbstractInitialStateConfig` subtype; see
`src/InitialState/norkyst_initial_state.jl`.
"""
function initial_state_record end

"""
    initial_state_source_grid(config, filepath)

Geometry of the source data in `filepath`, e.g. a `ProjectedSourceGrid`.

A new initial-state dataset implements this on its `AbstractInitialStateConfig` subtype.
"""
function initial_state_source_grid end

"""
    initial_state_variable_names(config)

Mapping from source variable name to the FjordSim state name it becomes, e.g.
`"temperature" => "T"`. Only the intersection with `config.parameters` is prepared.

A new initial-state dataset implements this on its `AbstractInitialStateConfig` subtype.
"""
function initial_state_variable_names end

"""
    initial_state_date(initial_state_config, run_config)

The instant the initial state is for: the run config's `start_date`. `nothing` for a
`nothing` initial state, which needs no date. A setup naming an initial state but no run config
has no such instant, which is an error rather than a no-op.
"""
initial_state_date(::Nothing, run_config) = nothing
initial_state_date(::AbstractInitialStateConfig, config::AbstractRunConfig) = start_date(config)

initial_state_date(::AbstractInitialStateConfig, ::Nothing) = error(
    "The initial state is read at the run config's `start_date`, and this setup names no " *
    "`run_config`.",
)

# --- Drivers ---

"""
    download_initial_state(config)

Fetch the source record a later `prepare_initial_state` reads, into `initial_state_directory`:
builds the bare domain grid on the CPU, takes `date` from the run config's `start_date` and hands
both to the initial-state config's `fetch_initial_state` hook. Returns `nothing` for a `nothing`
initial-state config.
"""
download_initial_state(config::AbstractFjordConfig) = fetch_initial_state(
    domain_grid(grid_config(config), CPU()),
    initial_state_config(config),
    initial_state_date(initial_state_config(config), run_config(config)),
)

"""
    fetch_initial_state(target_grid, initial_state_config, date)

The hook an initial-state dataset that downloads implements. A `nothing` initial-state config
fetches nothing.
"""
fetch_initial_state(target_grid, ::Nothing, date) = nothing

"""
    prepare_initial_state(config)

Regrid the initial state a setup names onto its simulation grid. Returns `nothing` for a `nothing`
initial-state config. Only picks out the five components this step reads and hands them to
`regrid_initial_state`.
"""
prepare_initial_state(config::AbstractFjordConfig) = regrid_initial_state(
    grid_config(config),
    bathymetry_config(config),
    initial_state_config(config),
    boundary_data_config(config),
    run_config(config),
)

"""
    regrid_initial_state(grid_config, bathymetry_config, initial_state_config, boundary_data_config, run_config)
    regrid_initial_state(target_grid, initial_state_config, date; edges = nothing)

Regrid the downloaded record for `date` onto `target_grid` and write it, with `date` as its
`$INITIAL_STATE_DATE_ATTRIBUTE` attribute, at `initial_state_path(config)`. Returns `nothing` for a
`nothing` initial-state config.

The file has the layout of the prepared forcing file with one record and no `_lambda` twins. As
there, `target_grid` is the immersed grid the simulation runs on, cells the model treats as land are
`NaN`, and `edges` are the open lateral boundaries whose velocity face rows are kept. The
component form reads the grid from the processed bathymetry, the date from the run config and
the edges from the boundary data config. Interpolation runs on the CPU: it is one record.

# Returns
A named tuple with `output_file`, `date` (the source record's own) and `variables`.
"""
regrid_initial_state(grid_config, bathymetry_config, ::Nothing, boundary_data_config, run_config) =
    nothing

function regrid_initial_state(
    grid_config,
    bathymetry_config,
    initial_state_config::AbstractInitialStateConfig,
    boundary_data_config,
    run_config,
)
    grid = simulation_grid(grid_config, resolve_bathymetry_file(bathymetry_config), CPU())
    result = regrid_initial_state(
        grid,
        initial_state_config,
        initial_state_date(initial_state_config, run_config);
        edges = open_edges(boundary_data_config),
    )

    @info "Prepared variables: $(join(result.variables, ", "))"
    @info "Initial state saved to $(result.output_file)"

    return result
end

function regrid_initial_state(target_grid, config::AbstractInitialStateConfig, date; edges = nothing)
    variable_names = initial_state_variable_names(config)
    names = [name => variable_names[name] for name in config.parameters if haskey(variable_names, name)]
    isempty(names) && error(
        "None of the configured parameters $(config.parameters) map to a state variable. " *
        "Known $(nameof(typeof(config))) variables: $(sort(collect(keys(variable_names)))).",
    )

    record = initial_state_record(config, date)
    source = initial_state_source_grid(config, record.filepath)
    output_file = initial_state_path(config)
    @info "Preparing the initial state for $date from $(basename(record.filepath)), " *
          "whose record is at $(record.date)"

    write_initial_state_file(output_file, target_grid, names, record, source, date, edges)

    return (; output_file, date = record.date, variables = last.(names))
end

"""
    write_initial_state_file(filepath, target_grid, names, record, source, date, edges)

Interpolate `record` of each `source name => state name` in `names` onto `target_grid` and write it
at `filepath`.
"""
function write_initial_state_file(filepath, target_grid, names, record, source, date, edges)
    mkpath(dirname(filepath))
    source_field = Field{Center,Center,Center}(source_field_grid(source))
    reader = SourceReader(record.filepath)

    try
        NCDataset(filepath, "c") do ds
            define_forcing_dimensions!(ds, target_grid, [record.date])
            ds.attrib[INITIAL_STATE_DATE_ATTRIBUTE] = string(date)

            for (source_name, name) in names
                LX, LY, LZ = data_location(Symbol(name))
                mask = Array{Bool,3}(water_mask(target_grid, LX, LY, edges))
                longitude = Array(λnodes(target_grid, LX()))
                latitude = Array(φnodes(target_grid, LY()))
                shape = (length(longitude), length(latitude))
                x, y = projected_target_nodes(longitude, latitude, source)
                z = source_level_coordinate(source, znodes(target_grid, LZ()))

                mask_fill = source_fill(source_validity(record.filepath, source_name))
                set_source_field!(source_field, source_slab(reader, record, source_name), mask_fill)
                output = zeros(Float32, size(mask))
                interpolate_to_target!(output, source_field, reshape(x, shape), reshape(y, shape), z, mask, CPU())

                variable = defVar(ds, name, Float32, forcing_dimension_names(name);
                    deflatelevel = FORCING_DEFLATE_LEVEL, attrib = ["_FillValue" => NaN32])
                variable[:, :, :, 1] = output
            end
        end
    finally
        close(reader)
    end

    return filepath
end

include("norkyst_initial_state.jl")

end # module
