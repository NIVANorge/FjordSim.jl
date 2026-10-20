# NorKyst-800m adapter: the source-specific half of the forcing pipeline. Copy this file's
# shape to add a new forcing dataset — it owns a config subtype, its own download, and the
# three hook methods `forcing_time_steps`, `forcing_source_grid` and `forcing_variable_names`
# (plus `download_forcing` if it fetches data). Everything else comes from the generic core
# in `Forcing.jl` and `Regridding.jl`, and what is NorKyst's own — its endpoints, grid and subset
# I/O — from `src/DataSources/NorKyst.jl`.

using Oceananigans.Grids: x_domain, y_domain

using ..NorKyst:
    NORKYST_CATALOG_URL,
    NORKYST_OPENDAP_URL,
    NORKYST_VARIABLE_NAMES,
    norkyst_source_grid,
    list_opendap_files,
    subset_ranges,
    time_length,
    define_output_file,
    write_parameter_chunk!,
    write_time_dependent_coordinates!

"""
    AbstractNorKystConfig <: AbstractForcingConfig

Supertype of the NorKyst-800m forcing configs. Two collections carry the same variables on the
same kind of polar-stereographic grid — the operational archive (`NorKystConfig`, 2017 onwards)
and the Norkyst-v3 hindcast (`NorKystHindcastConfig`, 2012 onwards) — so `forcing_variable_names`,
`forcing_time_steps` and `forcing_source_grid` are written once here. What differs is where the
files live and how they are named, which is what each subtype states.
"""
abstract type AbstractNorKystConfig <: AbstractForcingConfig end

"""
    NorKystConfig

Configuration for downloading and subsetting NorKyst-800m reanalysis data.

A setup states where its forcing goes, which variables to extract and which years to cover.
Only `catalog_url` and `opendap_url` are defaulted, being the NorKyst dataset's own public
endpoints.

`output_directory`, `output_file` and `plot_file` are names relative to `data_root`, resolved
by `forcing_directory`, `forcing_path` and `plot_path`. Setting one to an absolute path
overrides `data_root` for that entry only.

# Fields
- `data_root`: Directory holding this setup's forcing files. Required.
- `output_directory`: Name of the directory where monthly NetCDF files are written. Required.
- `output_file`: Name of the prepared forcing NetCDF written by `prepare_forcing`.
- `plot_file`: Name of the diagnostic forcing plot.
- `architecture`: Where `prepare_forcing` interpolates — `:auto`, `:cpu` or `:gpu`, resolved by
  `interpolation_architecture`. `:auto` keeps one config usable on a GPU machine and a laptop.
- `catalog_url`: THREDDS catalog URL listing available files.
- `opendap_url`: OPeNDAP base URL for streaming data.
- `parameters`: Variable names to extract (e.g. `["temperature", "salinity"]`). Required.
- `years`: Calendar years to download. Required.

The rivers and the open edges are not here: they belong to the `AbstractRiverConfig` and the
`AbstractBoundaryDataConfig` a setup names on its `FjordConfig`. `prepare_forcing` takes the edges
as a keyword.
"""
Base.@kwdef mutable struct NorKystConfig <: AbstractNorKystConfig
    data_root::String
    output_directory::String
    output_file::String = "forcing.nc"
    plot_file::String = "forcing.png"
    architecture::Symbol = :auto
    catalog_url::String = NORKYST_CATALOG_URL
    opendap_url::String = NORKYST_OPENDAP_URL
    parameters::Vector{String}
    years::Vector{Int}
end

"""
    forcing_monthly_filename(config::NorKystConfig, year, month)

Name of the combined monthly NorKyst-800m NetCDF file written for `year` and `month`.
"""
forcing_monthly_filename(config::NorKystConfig, year, month) =
    "NorKyst-800m_ZDEPTHS_avg_$(year)$(lpad(month, 2, '0')).nc"

"""
    forcing_variable_names(config::AbstractNorKystConfig)

The NorKyst variables this dataset can supply and the FjordSim forcing names they become.
"""
forcing_variable_names(config::AbstractNorKystConfig) = NORKYST_VARIABLE_NAMES

"""
    forcing_time_steps(config::AbstractNorKystConfig)

Every time record of every downloaded monthly file for `config.years`, sorted by date with
duplicates dropped. Errors if the directory or the files are missing.
"""
function forcing_time_steps(config::AbstractNorKystConfig)
    directory = forcing_directory(config)
    isdir(directory) || error(
        "NorKyst directory $directory does not exist. " *
        "Run `julia --project -m FjordSim download_forcing` for this setup first.",
    )

    records = SourceRecord[]
    for year in config.years, month = 1:12
        filepath = joinpath(directory, forcing_monthly_filename(config, year, month))
        isfile(filepath) || continue
        NCDataset(filepath) do ds
            for (index, date) in enumerate(ds["time"][:])
                push!(records, SourceRecord(DateTime(date), filepath, index))
            end
        end
    end

    isempty(records) && error("No NorKyst monthly files for years $(config.years) found in $directory.")
    sort!(records; by = record -> record.date)

    return unique(record -> record.date, records)
end

"""
    forcing_source_grid(config::AbstractNorKystConfig, filepath)

`norkyst_source_grid(filepath)`.
"""
forcing_source_grid(config::AbstractNorKystConfig, filepath) = norkyst_source_grid(filepath)

# --- Download ---

"""
    download_forcing(target_grid, config::NorKystConfig)

Download the NorKyst-800m months covering `config.years`, each combined into one NetCDF file
in `forcing_directory(config)` subset to the lon/lat box of `target_grid`. A month whose file
already exists is skipped, so an interrupted download resumes.
"""
function download_forcing(target_grid, config::NorKystConfig)
    output_directory = forcing_directory(config)
    mkpath(output_directory)

    @info "Downloading NorKyst-800m years $(join(config.years, ", ")) to $output_directory"

    files = list_opendap_files(catalog_url = config.catalog_url)
    for year in config.years, month = 1:12
        process_month(year, month, target_grid, config; files)
    end

    @info "Finished downloading NorKyst-800m forcing"
    return output_directory
end

function process_month(year, month, target_grid, config::NorKystConfig; files = nothing)
    files = isnothing(files) ? list_opendap_files(catalog_url = config.catalog_url) : files
    month_string = ".$(year)$(lpad(month, 2, '0'))"
    year_month = month_string[2:end]
    output_path = joinpath(forcing_directory(config), forcing_monthly_filename(config, year, month))

    if isfile(output_path)
        @info "Skipping $year_month (already exists)"
        return output_path
    end

    @info "Processing $year_month..."
    monthly_files = sort([file for file in files if occursin(month_string, file)])

    if isempty(monthly_files)
        @info "  No files found for $year_month, skipping."
        return nothing
    end

    urls = [joinpath(config.opendap_url, file) for file in monthly_files]
    @info "  Opening $(length(urls)) datasets..."

    datasets = NCDataset[]
    try
        for url in urls
            push!(datasets, NCDataset(url))
            @info "    Opened: $url"
        end

        longitude, latitude = x_domain(target_grid), y_domain(target_grid)
        subset = subset_ranges(first(datasets), longitude, latitude, config.parameters)
        total_time = sum(time_length, datasets)

        @info "  Writing output to: $output_path"
        # Written under a temporary name and renamed once complete, so a month interrupted
        # mid-write is redone on the next run rather than skipped as already present.
        temporary_path = output_path * ".tmp"
        output = define_output_file(temporary_path, first(datasets), subset, total_time)
        try
            time_start = 1
            for ds in datasets
                write_time_dependent_coordinates!(output, ds, subset, time_start)
                for name in subset.parameters
                    write_parameter_chunk!(output, ds, name, subset, time_start)
                end
                time_start += time_length(ds)
            end
        finally
            close(output)
        end
        mv(temporary_path, output_path; force = true)
    finally
        foreach(close, datasets)
    end

    @info "Finished $year_month"
    return output_path
end
