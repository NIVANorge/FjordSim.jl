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
    AbstractSimulationConfig,
    FjordConfig
using ..Grids: domain_grid, open_edges, simulation_grid
using ..Bathymetry: bathymetry_path
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
    initial_state_date(simulation_config)

The instant the initial state is for: the simulation config's `start_date`. A setup naming no
simulation config has no such instant, which is an error rather than a no-op, since it does name an
initial state.
"""
initial_state_date(config::AbstractSimulationConfig) = config.start_date

initial_state_date(::Nothing) = error(
    "The initial state is read at the simulation config's `start_date`, and this setup names no " *
    "`simulation_config`.",
)

# --- Drivers ---

"""
    download_initial_state(config::FjordConfig)
    download_initial_state(target_grid, config::AbstractInitialStateConfig, date)

Fetch the source record a later `prepare_initial_state` reads, into `initial_state_directory`.
Returns `nothing` when the setup names no `initial_state_config`.

The `FjordConfig` form is the generic driver: it builds the bare domain grid on the CPU and takes
`date` from the simulation config's `start_date`, so a new dataset only implements the second.
"""
function download_initial_state(config::FjordConfig)
    isnothing(config.initial_state_config) && return nothing

    return download_initial_state(
        domain_grid(config.grid_config, CPU()),
        config.initial_state_config,
        initial_state_date(config.simulation_config),
    )
end

"""
    prepare_initial_state(config::FjordConfig)
    prepare_initial_state(target_grid, config::AbstractInitialStateConfig, date; edges = nothing)

Regrid the downloaded record for `date` onto `target_grid` and write it, with `date` as its
`$INITIAL_STATE_DATE_ATTRIBUTE` attribute, at `initial_state_path(config)`. Returns `nothing` when
the setup names no `initial_state_config`.

The file has the layout of the prepared forcing file with one record and no `_lambda` twins. As
there, `target_grid` is the immersed grid the simulation runs on, cells the model treats as land are
`NaN`, and `edges` are the open lateral boundaries whose velocity face rows are kept. The
`FjordConfig` form reads the grid from the processed bathymetry, the date from the simulation config
and the edges from the boundary config. Interpolation runs on the CPU: it is one record.

# Returns
A named tuple with `output_file`, `date` (the source record's own) and `variables`.
"""
function prepare_initial_state(config::FjordConfig)
    isnothing(config.initial_state_config) && return nothing

    bathymetry_file = bathymetry_path(config.bathymetry_config)
    isfile(bathymetry_file) || error(
        "Processed bathymetry $bathymetry_file does not exist. " *
        "Run `julia --project -m FjordSim prepare_bathymetry` for this setup first.",
    )

    grid = simulation_grid(config.grid_config, bathymetry_file, CPU())
    result = prepare_initial_state(
        grid,
        config.initial_state_config,
        initial_state_date(config.simulation_config);
        edges = open_edges(config.boundary_config),
    )

    @info "Prepared variables: $(join(result.variables, ", "))"
    @info "Initial state saved to $(result.output_file)"

    return result
end

function prepare_initial_state(target_grid, config::AbstractInitialStateConfig, date; edges = nothing)
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
