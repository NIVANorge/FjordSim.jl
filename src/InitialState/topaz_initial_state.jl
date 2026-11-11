# TOPAZ5-ECOSMO adapter for the initial state: the daily-mean biogeochemical state of the day
# `start_date` falls in (`DataSources/Topaz.jl`). It supplies no physics, so a setup prepares it
# beside a physical initial state and combines the two (NivaEcosmo's FjordSim extension does).

using Dates: Date, DateTime, @dateformat_str
import Dates
using Oceananigans.Grids: x_domain, y_domain

using ..Topaz:
    TOPAZ_CATALOG_URL,
    TOPAZ_OPENDAP_URL,
    TOPAZ_VARIABLE_NAMES,
    download_topaz,
    topaz_source_grid

"""
    TopazInitialStateConfig

The TOPAZ5-ECOSMO daily mean of the day `start_date` falls in, under the TOPAZ variable names and
units (`TOPAZ_VARIABLES`). The fjords are narrower than its 6.25 km cells, so most of a fjord
domain takes the nearest water TOPAZ has.

# Fields
- `data_root`: Directory holding this setup's initial-state files. Required.
- `output_directory`: Name of the directory the download is written to. Required.
- `output_file`: Name of the prepared NetCDF written by `prepare_initial_state`.
- `margin`: Degrees of longitude/latitude the download reaches past the domain on every side.
- `catalog_url`, `opendap_url`: The collection's public endpoints. Defaulted.
- `parameters`: TOPAZ variable names to extract. Required.
- `start_date`: The instant the state is for, independent of the run's. Required.
"""
Base.@kwdef mutable struct TopazInitialStateConfig <: AbstractInitialStateConfig
    data_root::String
    output_directory::String
    output_file::String = "initial_state_topaz.nc"
    margin::Float64 = 0.2
    catalog_url::String = TOPAZ_CATALOG_URL
    opendap_url::String = TOPAZ_OPENDAP_URL
    parameters::Vector{String}
    start_date::DateTime
end

initial_state_filename(config::TopazInitialStateConfig, date) =
    "TOPAZ5-ECOSMO_$(Dates.format(Date(date), dateformat"yyyymmdd")).nc"

initial_state_variable_names(config::TopazInitialStateConfig) = TOPAZ_VARIABLE_NAMES

initial_state_source_grid(config::TopazInitialStateConfig, filepath) = topaz_source_grid(filepath)

function initial_state_record(config::TopazInitialStateConfig, date)
    filepath = joinpath(initial_state_directory(config), initial_state_filename(config, date))
    isfile(filepath) || error(
        "No TOPAZ5-ECOSMO initial state for $(Date(date)) at $filepath. " *
        "Run `julia --project -m FjordSim download_initial_state` for this setup first.",
    )
    return NCDataset(filepath) do ds
        SourceRecord(DateTime(first(ds["time"][:])), filepath, 1)
    end
end

"""
    fetch_initial_state(target_grid, config::TopazInitialStateConfig, date)

Download the TOPAZ5-ECOSMO daily mean of the day `date` falls in, subset to `target_grid`'s box
widened by `margin`, into `initial_state_directory(config)`. A day already downloaded is skipped.
"""
function fetch_initial_state(target_grid, config::TopazInitialStateConfig, date)
    output_path = joinpath(initial_state_directory(config), initial_state_filename(config, date))
    if isfile(output_path)
        @info "Skipping $(Date(date)) (already exists)"
        return output_path
    end
    mkpath(initial_state_directory(config))

    longitude = (x_domain(target_grid)[1] - config.margin, x_domain(target_grid)[2] + config.margin)
    latitude = (y_domain(target_grid)[1] - config.margin, y_domain(target_grid)[2] + config.margin)
    @info "Downloading the TOPAZ5-ECOSMO daily mean for $(Date(date))"
    result = download_topaz(output_path, [Date(date)], longitude, latitude, config.parameters;
        catalog_url = config.catalog_url, opendap_url = config.opendap_url)
    isnothing(result) && error("TOPAZ5-ECOSMO publishes no file for $(Date(date)).")

    @info "Initial-state download saved to $output_path"
    return output_path
end
