# TOPAZ5-ECOSMO adapter for the open-boundary pipeline: the daily biogeochemical state along the
# open edges (`DataSources/Topaz.jl`). It supplies no physics, so a setup prepares it beside a
# physical boundary source and combines the two (NivaEcosmo's FjordSim extension does).

import Dates
using Dates: Date, Day, Millisecond
using ..Topaz:
    TOPAZ_CATALOG_URL,
    TOPAZ_OPENDAP_URL,
    TOPAZ_VARIABLE_NAMES,
    download_topaz,
    topaz_source_grid

"""
    TopazBoundaryDataConfig

Configuration for downloading and preparing daily TOPAZ5-ECOSMO data (nitrate, phosphate, silicate,
oxygen, DIC, chlorophyll, phytoplankton and zooplankton carbon, ...) along the open lateral
boundaries. The variables keep their TOPAZ names and units (`TOPAZ_VARIABLES`).

The prepared file keeps the daily cadence, so it is read beside an hourly physical boundary file
rather than as one. The download reaches one day past each end of the period, so the records
bracket it.

# Fields
- `data_root`: Directory holding this setup's boundary files. Required.
- `output_directory`: Name of the directory the downloads are written to. Required.
- `output_file`: Name of the prepared boundary NetCDF written by `prepare_boundary_data`.
- `plot_file`: Name of the diagnostic boundary plot.
- `open_edges`: Lateral boundaries the domain is open on. One `Symbol` or a collection. Required.
- `margin`: Degrees of latitude/longitude the download reaches past the boundary row. TOPAZ cells
  are 6.25 km, so it must span at least two of them.
- `architecture`: Where `prepare_boundary_data` interpolates — `:auto`, `:cpu` or `:gpu`.
- `catalog_url`, `opendap_url`: The collection's public endpoints. Defaulted.
- `parameters`: TOPAZ variable names to extract. Required.
- `start_date`, `stop_date`: The period downloaded and prepared, independent of the run's. Required.
"""
mutable struct TopazBoundaryDataConfig <: AbstractBoundaryDataConfig
    data_root::String
    output_directory::String
    output_file::String
    plot_file::String
    open_edges::Vector{Symbol}
    margin::Float64
    architecture::Symbol
    catalog_url::String
    opendap_url::String
    parameters::Vector{String}
    start_date::DateTime
    stop_date::DateTime
end

function TopazBoundaryDataConfig(;
    data_root,
    output_directory,
    output_file = "boundaries_topaz.nc",
    plot_file = "boundaries_topaz.png",
    open_edges,
    margin = 0.2,
    architecture = :auto,
    catalog_url = TOPAZ_CATALOG_URL,
    opendap_url = TOPAZ_OPENDAP_URL,
    parameters,
    start_date,
    stop_date,
)
    edges = lateral_edges(open_edges)
    isempty(edges) && throw(ArgumentError("TopazBoundaryDataConfig names no `open_edges`."))
    unknown = setdiff(parameters, keys(TOPAZ_VARIABLE_NAMES))
    isempty(unknown) || throw(ArgumentError("Not TOPAZ5-ECOSMO variables: $(join(unknown, ", "))."))

    return TopazBoundaryDataConfig(
        data_root,
        output_directory,
        output_file,
        plot_file,
        edges,
        Float64(margin),
        architecture,
        catalog_url,
        opendap_url,
        parameters,
        start_date,
        stop_date,
    )
end

boundary_variable_names(config::TopazBoundaryDataConfig) = TOPAZ_VARIABLE_NAMES

boundary_source_grid(config::TopazBoundaryDataConfig, filepath) = topaz_source_grid(filepath)

# Daily means: no tide to resolve, so the prepared file keeps their cadence.
boundary_record_spacing(config::TopazBoundaryDataConfig) = Day(1)

"""
    topaz_boundary_days(config)

The days downloaded for `config`'s period: one before its start and one after its end, so the
daily records, stamped at 12:00, bracket it.
"""
function topaz_boundary_days(config::TopazBoundaryDataConfig)
    first_day, last_day = Date(config.start_date), Date(config.stop_date - Millisecond(1))
    return collect((first_day - Day(1)):Day(1):(last_day + Day(1)))
end

topaz_boundary_filename(days) =
    "TOPAZ5-ECOSMO_boundary_$(Dates.format(first(days), "yyyymmdd"))-$(Dates.format(last(days), "yyyymmdd")).nc"

"""
    fetch_boundary_data(target_grid, config::TopazBoundaryDataConfig)

Download the period's days, one file per calendar month, subset to the box `boundary_domain`
derives from the open edges. A month whose file exists is skipped, so an interrupted download
resumes.
"""
function fetch_boundary_data(target_grid, config::TopazBoundaryDataConfig)
    directory = boundary_data_directory(config)
    mkpath(directory)
    edges = open_edges(config)
    longitude, latitude = boundary_domain(edges, target_grid, config.margin)
    @info "Downloading TOPAZ5-ECOSMO along the $(join(edges, ", ")) boundary to $directory"
    @info "  Band: longitude $longitude, latitude $latitude"

    days = topaz_boundary_days(config)
    for month_days in [filter(day -> Dates.yearmonth(day) == month, days) for month in unique(Dates.yearmonth.(days))]
        output_path = joinpath(directory, topaz_boundary_filename(month_days))
        if isfile(output_path)
            @info "Skipping $(basename(output_path)) (already exists)"
            continue
        end
        @info "Downloading $(first(month_days)) to $(last(month_days))"
        download_topaz(output_path, month_days, longitude, latitude, config.parameters;
            catalog_url = config.catalog_url, opendap_url = config.opendap_url)
    end

    return directory
end

function boundary_time_steps(config::TopazBoundaryDataConfig)
    directory = boundary_data_directory(config)
    isdir(directory) || error(
        "TOPAZ boundary directory $directory does not exist. " *
        "Run `julia --project -m FjordSim download_boundary_data` for this setup first.",
    )

    days = topaz_boundary_days(config)
    files = filter(name -> startswith(name, "TOPAZ5-ECOSMO_boundary_") && endswith(name, ".nc"), readdir(directory))
    records = SourceRecord[]
    for name in files
        filepath = joinpath(directory, name)
        NCDataset(filepath) do ds
            for (index, date) in enumerate(ds["time"][:])
                Date(date) in days && push!(records, SourceRecord(DateTime(date), filepath, index))
            end
        end
    end

    isempty(records) && error("No TOPAZ5-ECOSMO boundary records for $(first(days)) to $(last(days)) in $directory.")
    sort!(records; by = record -> record.date)
    return unique(record -> record.date, records)
end
