# Norkyst-v3 hindcast adapter: the same NorKyst-800m fields as `norkyst.jl`, from the collection
# that reaches back to 2012 instead of the operational archive that starts 2017-02-20. Everything
# that reads or writes a subset is shared through `AbstractNorKystConfig`; this file states only
# what is different — where the files are, what they are called, and that a day is 24 hourly
# records rather than one daily mean.

const NORKYST_HINDCAST_CATALOG_URL =
    "https://thredds.met.no/thredds/catalog/romshindcast/norkyst_v3/zdepth"
const NORKYST_HINDCAST_OPENDAP_URL =
    "https://thredds.met.no/thredds/dodsC/romshindcast/norkyst_v3/zdepth"

# The record of each daily file to keep. The hindcast publishes 24 hourly instants per day where
# the operational collection publishes one daily mean, and the interior forcing wants a daily
# cadence: `prepare_forcing` writes `λ = 0` everywhere, so the file's values are read only for
# initial conditions and as the carrier `add_rivers` patches. Keeping every hour would multiply
# both the download and the prepared file by 24 for no gain. Hour 12 is chosen rather than a daily
# mean because it reproduces the operational product's own 12:00 stamp exactly, which is what the
# `pad_time_steps` arithmetic in the setups is written against.
const NORKYST_HINDCAST_HOUR = 12

"""
    NorKystHindcastConfig

Configuration for downloading and subsetting the Norkyst-v3 hindcast, MET Norway's continuous
NorKyst-800m free run from 2012 onwards.

Identical in role to [`NorKystConfig`](@ref), and carrying the same variables on the same
polar-stereographic 800 m grid, so the two share every hook except the ones below. Reach for this
one whenever the run window predates 2017-02-20, which is where the operational archive
`NorKystConfig` reads begins.

Two differences are worth knowing before using it:

- **The files are per day, under `YYYY/MM/`**, rather than one flat catalog of monthly-tagged
  names, so `list_hindcast_files` walks a month's sub-catalog instead of filtering a list.
- **Each file holds 24 hourly instants**, not one daily mean. Only hour `NORKYST_HINDCAST_HOUR` is
  kept, which is what makes the prepared file the same shape and cadence as the operational one.

`output_directory`, `output_file` and `plot_file` are names relative to `data_root`, resolved by
`forcing_directory`, `forcing_path` and `plot_path`. Setting one to an absolute path overrides
`data_root` for that entry only.

# Fields
- `data_root`: Directory holding this setup's forcing files. Required.
- `output_directory`: Name of the directory where monthly NetCDF files are written. Required.
- `output_file`: Name of the prepared forcing NetCDF written by `prepare_forcing`.
- `plot_file`: Name of the diagnostic forcing plot.
- `architecture`: Where `prepare_forcing` interpolates — `:auto`, `:cpu` or `:gpu`.
- `catalog_url`: THREDDS catalog base URL; `YYYY/MM/catalog.xml` is appended per month.
- `opendap_url`: OPeNDAP base URL; `YYYY/MM/<file>` is appended per day.
- `parameters`: Variable names to extract. Required.
- `years`: Calendar years to download. Required.
- `rivers`: An `AbstractRiverConfig` for `add_rivers`, or `nothing`.
"""
Base.@kwdef mutable struct NorKystHindcastConfig{R} <: AbstractNorKystConfig
    data_root::String
    output_directory::String
    output_file::String = "forcing.nc"
    plot_file::String = "forcing.png"
    architecture::Symbol = :auto
    catalog_url::String = NORKYST_HINDCAST_CATALOG_URL
    opendap_url::String = NORKYST_HINDCAST_OPENDAP_URL
    parameters::Vector{String}
    years::Vector{Int}
    rivers::R = nothing
end

"""
    forcing_monthly_filename(config::NorKystHindcastConfig, year, month)

Name of the combined monthly NetCDF written for `year` and `month`. Distinct from the operational
collection's name so the two can share one directory without colliding.
"""
forcing_monthly_filename(config::NorKystHindcastConfig, year, month) =
    "Norkyst-v3_ZDEPTHS_$(year)$(lpad(month, 2, '0')).nc"

"""
    list_hindcast_files(catalog_url, year, month)

The daily file names the hindcast publishes for one month, sorted. Walks that month's own THREDDS
sub-catalog, since the hindcast has no flat catalog to filter. An absent month returns empty rather
than raising: the collection's ends are ragged, and `process_month` reports the gap.
"""
function list_hindcast_files(catalog_url, year, month)
    url = join((catalog_url, string(year), lpad(month, 2, '0'), "catalog.xml"), "/")

    catalog = try
        read(Downloads.download(url), String)
    catch exception
        exception isa InterruptException && rethrow()
        @warn "Could not read the Norkyst-v3 catalog for $year-$(lpad(month, 2, '0'))" url exception
        return String[]
    end

    files = String[]
    for match in eachmatch(r"<dataset\b[^>]*\bname=\"([^\"]+\.nc)\""i, catalog)
        push!(files, match.captures[1])
    end

    return sort!(unique!(files))
end

"""
    hindcast_urls(config, year, month)

OPeNDAP URLs for every daily file of `year`/`month`, in date order.
"""
function hindcast_urls(config, year, month)
    directory = join((config.opendap_url, string(year), lpad(month, 2, '0')), "/")
    return [join((directory, file), "/") for file in list_hindcast_files(config.catalog_url, year, month)]
end

"""
    hindcast_time_index(ds)

The record of `ds` stamped `NORKYST_HINDCAST_HOUR`, or `nothing` if the file does not carry that
hour. Selected by the timestamp rather than by position so a short or shifted file is skipped
instead of contributing the wrong hour.
"""
function hindcast_time_index(ds)
    index = findfirst(date -> Dates.hour(DateTime(date)) == NORKYST_HINDCAST_HOUR, ds["time"][:])
    return index
end

# --- Download ---

"""
    download_forcing(target_grid, config::NorKystHindcastConfig)

Download the Norkyst-v3 hindcast months covering `config.years`, each combined into one NetCDF file
in `forcing_directory(config)` subset to the lon/lat box of `target_grid` and to one record a day.
A month whose file already exists is skipped, so an interrupted download resumes.
"""
function download_forcing(target_grid, config::NorKystHindcastConfig)
    output_directory = forcing_directory(config)
    mkpath(output_directory)

    @info "Downloading Norkyst-v3 hindcast years $(join(config.years, ", ")) to $output_directory"

    for year in config.years, month = 1:12
        process_hindcast_month(year, month, target_grid, config)
    end

    @info "Finished downloading Norkyst-v3 hindcast forcing"
    return output_directory
end

"""
    process_hindcast_month(year, month, target_grid, config::NorKystHindcastConfig)

`process_month`'s counterpart for the hindcast: one source file per day instead of per month, and
one record kept from each.

The time subsetting rides on machinery that is already there rather than adding any. `subset.ranges`
is consulted by `variable_indices` for *every* dimension, so putting the wanted time record in it
under the time dimension's own name makes `write_parameter_chunk!` and
`write_time_dependent_coordinates!` read that record alone — from OPeNDAP, so the 23 unused hours
are never transferred. Each day therefore contributes exactly one record, which is why `time_start`
advances by one rather than by `time_length`.
"""
function process_hindcast_month(year, month, target_grid, config::NorKystHindcastConfig)
    year_month = "$(year)-$(lpad(month, 2, '0'))"
    output_path = joinpath(forcing_directory(config), forcing_monthly_filename(config, year, month))

    if isfile(output_path)
        @info "Skipping $year_month (already exists)"
        return output_path
    end

    urls = hindcast_urls(config, year, month)
    if isempty(urls)
        @info "  No Norkyst-v3 files found for $year_month, skipping."
        return nothing
    end

    @info "Processing $year_month ($(length(urls)) daily files)..."

    datasets = NCDataset[]
    indices = Int[]
    try
        for url in urls
            ds = NCDataset(url)
            index = hindcast_time_index(ds)

            if isnothing(index)
                @warn "  No $(NORKYST_HINDCAST_HOUR):00 record, skipping" url
                close(ds)
                continue
            end

            push!(datasets, ds)
            push!(indices, index)
        end

        if isempty(datasets)
            @info "  No usable Norkyst-v3 records for $year_month, skipping."
            return nothing
        end

        spatial_subset = subset_ranges(first(datasets), target_grid, config)
        time_name = time_dimension(first(datasets))

        @info "  Writing output to: $output_path"
        output = define_output_file(output_path, first(datasets), spatial_subset, length(datasets))
        try
            for (time_start, (ds, index)) in enumerate(zip(datasets, indices))
                subset = NorKystSubset(
                    merge(spatial_subset.ranges, Dict(time_name => index:index)),
                    spatial_subset.mask,
                    spatial_subset.spatial_dimensions,
                    spatial_subset.parameters,
                )

                write_time_dependent_coordinates!(output, ds, subset, time_start)
                for name in subset.parameters
                    write_parameter_chunk!(output, ds, name, subset, time_start)
                end
            end
        finally
            close(output)
        end
    finally
        foreach(close, datasets)
    end

    @info "Finished $year_month"
    return output_path
end
