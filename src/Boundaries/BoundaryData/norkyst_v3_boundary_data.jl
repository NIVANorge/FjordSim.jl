# Operational Norkyst-v3 adapter for the open-boundary pipeline: the same exterior state as
# `norkyst_boundary_data.jl`, from MET's operational system since 2024-01-01, where the operational
# NorKyst-800m archive ends on 2025-10-05. Shares every read-side hook with the other NorKyst
# configs through `AbstractNorKystBoundaryDataConfig`; what this file states is where the files are
# and how the barotropic pair is obtained.
#
# One file a day per member, 24 hourly instants on 15 z levels down to 300 m: `temperature`,
# `salinity`, `u_eastward`, `v_northward` (already eastward/northward) and `zeta`, but no
# barotropic velocity at all. The download therefore writes the pair itself, as the depth mean of
# the downloaded `u_eastward`/`v_northward` over each column's `h` (`depth_mean`), under the
# hindcast's names `ubar_eastward`/`vbar_northward`. Reading it back is then the generic
# `blended_slab`, with no `boundary_source_slab` method. The mean is over `h` rather than `h + ζ`,
# which differs by under 1 % where the column is 100 m deep and the tide a metre. Measured against
# the Norkyst-v3 hindcast's own barotropic pair along the Oslofjord's southern edge (58.95-59.05°N,
# 10.2-11.02°E, 2024-10-01, 3156 points), the depth mean of its z levels cut to these 15 correlates
# at 0.986 (east) and 0.991 (north), RMSE 1.3 cm/s against a 7-9 cm/s spread and 1 mm/s bias, no
# different from the mean over all 25 of its levels.

using Dates: Date
using ..Regridding: depth_mean
using ..NorKyst:
    NORKYST_V3_CATALOG_URL,
    NORKYST_V3_OPENDAP_URL,
    list_v3_days,
    norkyst_v3_url,
    open_v3_file,
    subset_ranges,
    time_length,
    define_output_file,
    write_parameter_chunk!,
    write_time_dependent_coordinates!

# The members are published at two resolutions; the boundary reads the 800 m parent of whichever.
const NORKYST_V3_BOUNDARY_RESOLUTION = "800"

# The barotropic pair the download derives, and the z-level component each is the depth mean of.
const NORKYST_V3_DEPTH_MEANS = Dict("ubar_eastward" => "u_eastward", "vbar_northward" => "v_northward")

const NORKYST_V3_BOUNDARY_VARIABLE_NAMES = Dict(
    "temperature" => "T",
    "salinity" => "S",
    "u_eastward" => "u",
    "v_northward" => "v",
    "zeta" => "eta",
    "ubar_eastward" => "ubar",
    "vbar_northward" => "vbar",
)

"""
    NorKystV3BoundaryDataConfig

Configuration for downloading and subsetting hourly operational Norkyst-v3 data along the open
lateral boundaries, from 2024-01-01 onward. The counterpart of [`NorKystBoundaryDataConfig`](@ref),
whose archive ends on 2025-10-05.

`ubar_eastward`/`vbar_northward` are not in the source: naming one in `parameters` has the download
write it as the depth mean of `u_eastward`/`v_northward`, which must then be named too.

# Fields
- `data_root`: Directory holding this setup's boundary files. Required.
- `output_directory`: Name of the directory the monthly downloads are written to. Required.
- `output_file`: Name of the prepared boundary NetCDF written by `prepare_boundary_data`.
- `plot_file`: Name of the diagnostic boundary plot.
- `open_edges`: Lateral boundaries the domain is open on, each one of `:south`, `:north`, `:west` or
  `:east`. One `Symbol` or a collection; the constructor normalizes either to a `Vector{Symbol}`.
  Required, and stated only here.
- `margin`: Degrees of latitude/longitude the download reaches past the boundary row, on both sides.
- `architecture`: Where `prepare_boundary_data` interpolates — `:auto`, `:cpu` or `:gpu`.
- `catalog_url`, `opendap_url`: The collection's public endpoints. Defaulted.
- `member`: The member whose 800 m file is read: `"m00"`, the reference run and the only one
  published from 2024-01-01, or the two-way parent of a 160 m nest, `"m70"` (Sulafjorden) or
  `"m71"` (Oslofjord, from 2024-11-23). Required.
- `parameters`: Source variable names to extract. Required.
- `start_date`, `stop_date`: The period downloaded and prepared, independent of the run's: every
  month it touches is downloaded, its end exclusive (`coverage_months`), and the prepared axis is
  padded to its ends. Required.
"""
mutable struct NorKystV3BoundaryDataConfig <: AbstractNorKystBoundaryDataConfig
    data_root::String
    output_directory::String
    output_file::String
    plot_file::String
    open_edges::Vector{Symbol}
    margin::Float64
    architecture::Symbol
    catalog_url::String
    opendap_url::String
    member::String
    parameters::Vector{String}
    start_date::DateTime
    stop_date::DateTime
end

# Hand-written for the reason `NorKystBoundaryDataConfig`'s is: `open_edges` may be one `Symbol`.
function NorKystV3BoundaryDataConfig(;
    data_root,
    output_directory,
    output_file = "boundaries.nc",
    plot_file = "boundaries.png",
    open_edges,
    margin = 0.05,
    architecture = :auto,
    catalog_url = NORKYST_V3_CATALOG_URL,
    opendap_url = NORKYST_V3_OPENDAP_URL,
    member,
    parameters,
    start_date,
    stop_date,
)
    edges = lateral_edges(open_edges)
    isempty(edges) && throw(
        ArgumentError(
            "NorKystV3BoundaryDataConfig names no `open_edges`. A boundary dataset exists to supply " *
            "the exterior state along at least one open edge; a closed domain names no boundary " *
            "config at all.",
        ),
    )
    for (mean, component) in NORKYST_V3_DEPTH_MEANS
        mean in parameters && component ∉ parameters && throw(
            ArgumentError("`$mean` is the depth mean of `$component`, so `parameters` must name it too."),
        )
    end

    return NorKystV3BoundaryDataConfig(
        data_root,
        output_directory,
        output_file,
        plot_file,
        edges,
        Float64(margin),
        architecture,
        catalog_url,
        opendap_url,
        member,
        parameters,
        start_date,
        stop_date,
    )
end

"""
    boundary_monthly_filename(config::NorKystV3BoundaryDataConfig, year, month)

Name of the combined monthly boundary NetCDF written for `year` and `month`, distinct from the other
NorKyst configs' and per member, so they can share one directory.
"""
boundary_monthly_filename(config::NorKystV3BoundaryDataConfig, year, month) =
    "Norkyst-v3-operational_$(config.member)_boundary_$(year)$(lpad(month, 2, '0')).nc"

boundary_variable_names(config::NorKystV3BoundaryDataConfig) = NORKYST_V3_BOUNDARY_VARIABLE_NAMES

# --- Download ---

"""
    fetch_boundary_data(target_grid, config::NorKystV3BoundaryDataConfig)

Download the operational Norkyst-v3 months the config's period touches, each combined into one NetCDF
file in `boundary_data_directory(config)` subset to the box `boundary_domain` derives from
`open_edges(config)` and `target_grid`. A month whose file already exists is skipped, so an
interrupted download resumes.
"""
function fetch_boundary_data(target_grid, config::NorKystV3BoundaryDataConfig)
    output_directory = boundary_data_directory(config)
    mkpath(output_directory)

    edges = open_edges(config)
    longitude, latitude = boundary_domain(edges, target_grid, config.margin)
    months = coverage_months(coverage_window(config))
    @info "Downloading operational Norkyst-v3 ($(config.member)) $(month_span(months)) " *
          "along the $(join(edges, ", ")) boundary to $output_directory"
    @info "  Band: longitude $longitude, latitude $latitude"

    for (year, month) in months
        process_v3_boundary_month(year, month, longitude, latitude, config)
    end

    @info "Finished downloading operational Norkyst-v3 boundary data"
    return output_directory
end

"""
    process_v3_boundary_month(year, month, longitude, latitude, config)

Download one month of hourly operational Norkyst-v3 into a single NetCDF, subset to the boundary
band, with the depth-mean velocities `config.parameters` names written after each day's records.
"""
function process_v3_boundary_month(year, month, longitude, latitude, config::NorKystV3BoundaryDataConfig)
    year_month = "$(year)-$(lpad(month, 2, '0'))"
    output_path =
        joinpath(boundary_data_directory(config), boundary_monthly_filename(config, year, month))

    if isfile(output_path)
        @info "Skipping $year_month (already exists)"
        return output_path
    end

    days = list_v3_days(config.catalog_url, year, month)
    if isempty(days)
        @info "  No Norkyst-v3 files found for $year_month, skipping."
        return nothing
    end

    means = [name for name in config.parameters if haskey(NORKYST_V3_DEPTH_MEANS, name)]
    downloaded = [name for name in config.parameters if !haskey(NORKYST_V3_DEPTH_MEANS, name)]
    @info "Processing $year_month ($(length(days)) daily files)..."

    datasets = NCDataset[]
    try
        for day in days
            date = Date(year, month, day)
            url = norkyst_v3_url(config.opendap_url, NORKYST_V3_BOUNDARY_RESOLUTION, config.member, date)
            push!(datasets, open_v3_file(url, date))
        end

        subset = subset_ranges(first(datasets), longitude, latitude, downloaded)
        total_time = sum(time_length, datasets)
        @info "  Subset $(join(("$dimension=$(length(range))" for (dimension, range) in sort(collect(subset.ranges))), ", ")), $total_time records"

        @info "  Writing output to: $output_path"
        # Written under a temporary name and renamed once complete, so a month interrupted
        # mid-write is redone on the next run rather than skipped as already present.
        temporary_path = output_path * ".tmp"
        output = define_output_file(temporary_path, first(datasets), subset, total_time)
        try
            for name in means
                define_depth_mean!(output, name, NORKYST_V3_DEPTH_MEANS[name])
            end

            time_start = 1
            for ds in datasets
                write_time_dependent_coordinates!(output, ds, subset, time_start)
                for name in downloaded
                    write_parameter_chunk!(output, ds, name, subset, time_start)
                end
                records = time_start:(time_start + time_length(ds) - 1)
                for name in means
                    write_depth_mean!(output, name, NORKYST_V3_DEPTH_MEANS[name], records)
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

"""
    define_depth_mean!(output, name, component)

Define `name`, the depth mean of the z-level variable `component`, on `component`'s dimensions less
`depth`.
"""
function define_depth_mean!(output, name, component)
    dimensions = Tuple(dimension for dimension in dimnames(output[component]) if dimension != "depth")
    defVar(
        output,
        name,
        Float32,
        dimensions;
        deflatelevel = 5,
        attrib = [
            "units" => "meter second-1",
            "long_name" => "Depth mean of $component over the water column, from its z levels",
        ],
    )
    return output
end

"""
    write_depth_mean!(output, name, component, records)

Write `name` at `records` as the `depth_mean` of `component` already written there, over the
subset's own bottom depth `h`.
"""
function write_depth_mean!(output, name, component, records)
    dimnames(output[component]) == ("X", "Y", "depth", "time") || error(
        "Expected `$component` on (X, Y, depth, time), found $(dimnames(output[component])).",
    )
    values = coalesce.(output[component][:, :, :, records], NaN32)
    bottom_depth = coalesce.(output["h"][:, :], NaN)
    output[name][:, :, records] = depth_mean(values, Float64.(output["depth"][:]), bottom_depth)
    return output
end
