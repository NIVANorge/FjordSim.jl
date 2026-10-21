# NorKyst-800m adapter for the initial state: one daily mean from the collection `NorKystConfig`
# reads, downloaded with the subsetting every NorKyst download shares
# (`src/DataSources/NorKyst.jl`). Copy this file's shape to add a new initial-state dataset — a
# config subtype, its download, and the three hooks.

using Dates: Date, DateTime, @dateformat_str
import Dates
using Oceananigans.Grids: x_domain, y_domain

using ..NorKyst:
    NORKYST_CATALOG_URL,
    NORKYST_OPENDAP_URL,
    NORKYST_VARIABLE_NAMES,
    norkyst_source_grid,
    list_opendap_files,
    open_remote,
    subset_ranges,
    time_length,
    define_output_file,
    write_parameter_chunk!,
    write_time_dependent_coordinates!

"""
    NorKystInitialStateConfig

The initial state from one NorKyst-800m daily mean: the one for the day `start_date` falls in,
downloaded subset to the domain box. The operational archive begins on 2017-02-20, as for
`NorKystConfig`.

A daily mean is stamped 12:00 of its day, so for a run starting at midnight the state is a mean
centred half a day after `start_date`.

`output_directory` and `output_file` are names relative to `data_root`, resolved by
`initial_state_directory` and `initial_state_path`. Setting one to an absolute path overrides
`data_root` for that entry only.

# Fields
- `data_root`: Directory holding this setup's initial-state files. Required.
- `output_directory`: Name of the directory the daily download is written to. Required.
- `output_file`: Name of the prepared NetCDF written by `prepare_initial_state`.
- `catalog_url`, `opendap_url`: The NorKyst dataset's own public endpoints; the only defaulted
  fields besides `output_file`.
- `parameters`: Variable names to extract (e.g. `["temperature", "salinity"]`). Required.
"""
Base.@kwdef mutable struct NorKystInitialStateConfig <: AbstractInitialStateConfig
    data_root::String
    output_directory::String
    output_file::String = "initial_state.nc"
    catalog_url::String = NORKYST_CATALOG_URL
    opendap_url::String = NORKYST_OPENDAP_URL
    parameters::Vector{String}
end

"""
    initial_state_filename(config::NorKystInitialStateConfig, date)

Name of the daily NorKyst-800m NetCDF downloaded for the day `date` falls in.
"""
initial_state_filename(config::NorKystInitialStateConfig, date) =
    "NorKyst-800m_ZDEPTHS_avg_$(Dates.format(date, dateformat"yyyymmdd")).nc"

initial_state_variable_names(config::NorKystInitialStateConfig) = NORKYST_VARIABLE_NAMES

initial_state_source_grid(config::NorKystInitialStateConfig, filepath) = norkyst_source_grid(filepath)

function initial_state_record(config::NorKystInitialStateConfig, date)
    filepath = joinpath(initial_state_directory(config), initial_state_filename(config, date))
    isfile(filepath) || error(
        "No NorKyst initial state for $(Date(date)) at $filepath. " *
        "Run `julia --project -m FjordSim download_initial_state` for this setup first.",
    )

    return NCDataset(filepath) do ds
        SourceRecord(DateTime(first(ds["time"][:])), filepath, 1)
    end
end

"""
    download_initial_state(target_grid, config::NorKystInitialStateConfig, date)

Download the NorKyst-800m daily file for the day `date` falls in, subset to the lon/lat box of
`target_grid`, into `initial_state_directory(config)`. A day already downloaded is skipped.
"""
function download_initial_state(target_grid, config::NorKystInitialStateConfig, date)
    output_path = joinpath(initial_state_directory(config), initial_state_filename(config, date))
    if isfile(output_path)
        @info "Skipping $(Date(date)) (already exists)"
        return output_path
    end
    mkpath(initial_state_directory(config))

    day = Dates.format(date, dateformat"yyyymmdd")
    files = sort([file for file in list_opendap_files(catalog_url = config.catalog_url) if occursin(".$day", file)])
    length(files) == 1 || error(
        "Expected one NorKyst-800m daily file for $(Date(date)) in $(config.catalog_url), found " *
        (isempty(files) ? "none; the archive begins on 2017-02-20." : "$(length(files)): $(join(files, ", ")).")
    )

    url = joinpath(config.opendap_url, only(files))
    @info "Downloading the NorKyst-800m daily mean for $(Date(date)) from $url"

    open_remote(url) do ds
        subset = subset_ranges(ds, x_domain(target_grid), y_domain(target_grid), config.parameters)
        # Written under a temporary name and renamed once complete, so an interrupted download is
        # redone on the next run rather than skipped as already present.
        temporary_path = output_path * ".tmp"
        output = define_output_file(temporary_path, ds, subset, time_length(ds))
        try
            write_time_dependent_coordinates!(output, ds, subset, 1)
            for name in subset.parameters
                write_parameter_chunk!(output, ds, name, subset, 1)
            end
        finally
            close(output)
        end
        mv(temporary_path, output_path; force = true)
    end

    @info "Initial-state download saved to $output_path"
    return output_path
end
