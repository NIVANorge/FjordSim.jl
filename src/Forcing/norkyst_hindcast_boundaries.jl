# Norkyst-v3 hindcast adapter for the open-boundary pipeline: the same exterior state as
# `norkyst_boundaries.jl`, from the collection that reaches back to 2012 instead of the operational
# one that starts 2017-02-20. Shares every read-side hook with it through
# `AbstractNorKystBoundariesConfig`; what this file states is where the files are and how the
# barotropic pair is obtained.
#
# The hindcast splits what the operational collection keeps in one file. `zdepth` carries the
# z-interpolated fields — `temperature`, `salinity`, `u_eastward`, `v_northward`, `zeta` — on 25
# levels, and `sdepth` carries the same run on its native 40 terrain-following levels. The
# barotropic pair lives only in `sdepth`, so a day's download opens one file from each. They share
# the grid exactly (both 2747 x 1148, X and Y both 0:800:… metres), so the two halves land in one
# subset with no regridding, and no s-level transform is needed because the barotropic fields are
# two-dimensional.
#
# The pair is `ubar_eastward`/`vbar_northward` here, already rotated to geographic axes by met.no,
# where the operational collection publishes ROMS' grid-relative `ubar`/`vbar` and
# `NorKystBoundariesConfig` has to derotate them itself. So this config wants no
# `boundary_source_slab` method at all, and the generic `blended_slab` default applies to all seven
# variables.

const NORKYST_HINDCAST_BOUNDARY_CATALOG_URL =
    "https://thredds.met.no/thredds/catalog/romshindcast/norkyst_v3/zdepth"
const NORKYST_HINDCAST_BOUNDARY_OPENDAP_URL =
    "https://thredds.met.no/thredds/dodsC/romshindcast/norkyst_v3/zdepth"
const NORKYST_HINDCAST_BAROTROPIC_CATALOG_URL =
    "https://thredds.met.no/thredds/catalog/romshindcast/norkyst_v3/sdepth"
const NORKYST_HINDCAST_BAROTROPIC_OPENDAP_URL =
    "https://thredds.met.no/thredds/dodsC/romshindcast/norkyst_v3/sdepth"

# The two variables that come from the `sdepth` file rather than the `zdepth` one.
const NORKYST_HINDCAST_BAROTROPIC_VARIABLES = ("ubar_eastward", "vbar_northward")

# Source variable names and the FjordSim boundary names they become. The mapping is what renames
# the derotated barotropic pair to the `ubar`/`vbar` the read side expects, so nothing downstream
# has to know which collection a variable came from.
const NORKYST_HINDCAST_BOUNDARY_VARIABLE_NAMES = Dict(
    "temperature" => "T",
    "salinity" => "S",
    "u_eastward" => "u",
    "v_northward" => "v",
    "zeta" => "eta",
    "ubar_eastward" => "ubar",
    "vbar_northward" => "vbar",
)

"""
    NorKystHindcastBoundariesConfig

Configuration for downloading and subsetting hourly Norkyst-v3 hindcast data along the open lateral
boundaries. The counterpart of [`NorKystBoundariesConfig`](@ref) for a run window predating
2017-02-20.

Every record is kept, unlike the interior forcing's one-a-day: a Flather boundary compares the
model's own elevation against the exterior one, so the boundary row is the one place the tide has to
survive, and a daily sample would alias it away.

# Fields
- `data_root`: Directory holding this setup's boundary files. Required.
- `output_directory`: Name of the directory the monthly downloads are written to. Required.
- `output_file`: Name of the prepared boundary NetCDF written by `prepare_boundaries`.
- `plot_file`: Name of the diagnostic boundary plot.
- `open_edges`: Lateral boundaries the domain is open on, each one of `:south`, `:north`, `:west` or
  `:east`. One `Symbol` or a collection; the constructor normalizes either to a `Vector{Symbol}`.
  Required, and stated only here.
- `margin`: Degrees of latitude/longitude the download reaches past the boundary row, on both sides.
  Must leave room for the source cells the bilinear interpolation reads around the row.
- `architecture`: Where `prepare_boundaries` interpolates — `:auto`, `:cpu` or `:gpu`.
- `catalog_url`, `opendap_url`: the `zdepth` collection, which supplies all five full-field
  variables. Defaulted.
- `barotropic_catalog_url`, `barotropic_opendap_url`: the `sdepth` collection, which supplies
  `ubar_eastward` and `vbar_northward` and nothing else. Defaulted. Stated separately rather than
  derived from a shared parent so that neither URL means something different from the other config's.
- `parameters`: Source variable names to extract. Required.
- `years`: Calendar years to download. Required.
"""
mutable struct NorKystHindcastBoundariesConfig <: AbstractNorKystBoundariesConfig
    data_root::String
    output_directory::String
    output_file::String
    plot_file::String
    open_edges::Vector{Symbol}
    margin::Float64
    architecture::Symbol
    catalog_url::String
    opendap_url::String
    barotropic_catalog_url::String
    barotropic_opendap_url::String
    parameters::Vector{String}
    years::Vector{Int}
end

# A hand-written keyword constructor for the same reason `NorKystBoundariesConfig` has one:
# `open_edges` may be written as one `Symbol` or as a collection, and `@kwdef` would `convert`,
# which has no method from `Symbol` to `Vector{Symbol}`.
function NorKystHindcastBoundariesConfig(;
    data_root,
    output_directory,
    output_file = "boundaries.nc",
    plot_file = "boundaries.png",
    open_edges,
    margin = 0.05,
    architecture = :auto,
    catalog_url = NORKYST_HINDCAST_BOUNDARY_CATALOG_URL,
    opendap_url = NORKYST_HINDCAST_BOUNDARY_OPENDAP_URL,
    barotropic_catalog_url = NORKYST_HINDCAST_BAROTROPIC_CATALOG_URL,
    barotropic_opendap_url = NORKYST_HINDCAST_BAROTROPIC_OPENDAP_URL,
    parameters,
    years,
)
    edges = lateral_edges(open_edges)
    isempty(edges) && throw(
        ArgumentError(
            "NorKystHindcastBoundariesConfig names no `open_edges`. A boundary dataset exists to " *
            "supply the exterior state along at least one open edge; a closed domain names no " *
            "boundary config at all.",
        ),
    )

    return NorKystHindcastBoundariesConfig(
        data_root,
        output_directory,
        output_file,
        plot_file,
        edges,
        Float64(margin),
        architecture,
        catalog_url,
        opendap_url,
        barotropic_catalog_url,
        barotropic_opendap_url,
        parameters,
        years,
    )
end

"""
    boundary_monthly_filename(config::NorKystHindcastBoundariesConfig, year, month)

Name of the combined monthly Norkyst-v3 boundary NetCDF written for `year` and `month`. Distinct
from the operational collection's name so the two can share one directory without colliding.
"""
boundary_monthly_filename(config::NorKystHindcastBoundariesConfig, year, month) =
    "Norkyst-v3_boundary_$(year)$(lpad(month, 2, '0')).nc"

"""
    boundary_variable_names(config::NorKystHindcastBoundariesConfig)

The Norkyst-v3 variables this dataset can supply along a boundary and the FjordSim names they
become. This is where `ubar_eastward`/`vbar_northward` become the `ubar`/`vbar` the read side reads.
"""
boundary_variable_names(config::NorKystHindcastBoundariesConfig) =
    NORKYST_HINDCAST_BOUNDARY_VARIABLE_NAMES

# --- Download ---

"""
    download_boundaries(target_grid, config::NorKystHindcastBoundariesConfig)

Download the Norkyst-v3 hindcast months covering `config.years`, each combined into one NetCDF file
in `boundary_data_directory(config)` subset to the box `boundary_domain` derives from
`open_edges(config)` and `target_grid`. A month whose file already exists is skipped, so an
interrupted download resumes.
"""
function download_boundaries(target_grid, config::NorKystHindcastBoundariesConfig)
    output_directory = boundary_data_directory(config)
    mkpath(output_directory)

    edges = open_edges(config)
    longitude, latitude = boundary_domain(edges, target_grid, config.margin)
    @info "Downloading Norkyst-v3 hindcast years $(join(config.years, ", ")) along the " *
          "$(join(edges, ", ")) boundary to $output_directory"
    @info "  Band: longitude $longitude, latitude $latitude"

    for year in config.years, month = 1:12
        process_hindcast_boundary_month(year, month, longitude, latitude, config)
    end

    @info "Finished downloading Norkyst-v3 hindcast boundary data"
    return output_directory
end

"""
    boundary_subset_ranges(ds, longitude, latitude, config::NorKystHindcastBoundariesConfig)

The Norkyst-v3 index window covering the `longitude`/`latitude` band, as the same `NorKystSubset`
every other download here builds — so `define_output_file`, `write_parameter_chunk!` and
`write_time_dependent_coordinates!` are reused verbatim.

`parameters` is passed in rather than taken from the config, because a day is read from two files
and each contributes only its own share of them.
"""
function boundary_subset_ranges(
    ds,
    longitude,
    latitude,
    config::NorKystHindcastBoundariesConfig,
    parameters,
)
    latitude_variable = variable(ds, "lat")
    longitude_variable = variable(ds, "lon")
    source_latitude = Array(latitude_variable[ntuple(_ -> :, ndims(latitude_variable))...])
    source_longitude = Array(longitude_variable[ntuple(_ -> :, ndims(longitude_variable))...])

    mask = (source_latitude .>= latitude[1]) .&
           (source_latitude .<= latitude[2]) .&
           (source_longitude .>= longitude[1]) .&
           (source_longitude .<= longitude[2])

    ranges = Dict{String,UnitRange{Int}}()
    for (index, dimension) in enumerate(dimnames(latitude_variable))
        ranges[dimension] = bounding_range(mask, index)
    end

    spatial_dimensions = dimnames(latitude_variable)
    subset_mask = mask[(ranges[dimension] for dimension in spatial_dimensions)...]
    return NorKystSubset(ranges, subset_mask, spatial_dimensions, parameters)
end

"""
    define_extra_parameters!(output, template, subset)

Define `subset.parameters` on an already-open output file from a second template.

`define_output_file` builds the file from one template, and a Norkyst-v3 day has two: the `zdepth`
file, which supplies the dimensions and the five full fields, and the `sdepth` file, which supplies
the barotropic pair alone. The pair is two-dimensional over `X`, `Y` and `time`, all of which the
`zdepth` template already defined, so only the variables themselves are missing.
"""
function define_extra_parameters!(output, template, subset::NorKystSubset)
    for name in subset.parameters
        haskey(output, name) && continue
        variable = template[name]
        defVar(
            output,
            name,
            nonmissingtype(eltype(variable)),
            dimnames(variable);
            deflatelevel = 5,
            attrib = decoded_attributes(variable),
        )
    end

    return output
end

"""
    process_hindcast_boundary_month(year, month, longitude, latitude, config)

Download one month of hourly Norkyst-v3 into a single NetCDF, subset to the boundary band.

Two source collections per day, so the driver differs from the single-collection ones in three
places: the output file is defined from the `zdepth` template and then extended with the `sdepth`
variables, each day's records are written from whichever file holds them, and the `sdepth` files are
opened one at a time inside the loop rather than all up front — the `zdepth` handles are already
held open for the whole month to total the time axis, and doubling that costs handles for nothing.
"""
function process_hindcast_boundary_month(
    year,
    month,
    longitude,
    latitude,
    config::NorKystHindcastBoundariesConfig,
)
    year_month = "$(year)-$(lpad(month, 2, '0'))"
    output_path =
        joinpath(boundary_data_directory(config), boundary_monthly_filename(config, year, month))

    if isfile(output_path)
        @info "Skipping $year_month (already exists)"
        return output_path
    end

    files = list_hindcast_files(config.catalog_url, year, month)
    if isempty(files)
        @info "  No Norkyst-v3 files found for $year_month, skipping."
        return nothing
    end

    month_directory = join((string(year), lpad(month, 2, '0')), "/")
    barotropic_names =
        [name for name in config.parameters if name in NORKYST_HINDCAST_BAROTROPIC_VARIABLES]
    full_names = [name for name in config.parameters if name ∉ NORKYST_HINDCAST_BAROTROPIC_VARIABLES]

    @info "Processing $year_month ($(length(files)) daily files)..."

    datasets = NCDataset[]
    try
        for file in files
            push!(datasets, NCDataset(join((config.opendap_url, month_directory, file), "/")))
        end

        full_subset =
            boundary_subset_ranges(first(datasets), longitude, latitude, config, full_names)
        total_time = sum(time_length, datasets)
        @info "  Subset $(join(("$dimension=$(length(range))" for (dimension, range) in sort(collect(full_subset.ranges))), ", ")), $total_time records"

        @info "  Writing output to: $output_path"
        output = define_output_file(output_path, first(datasets), full_subset, total_time)
        try
            barotropic_subset = NorKystSubset(
                full_subset.ranges,
                full_subset.mask,
                full_subset.spatial_dimensions,
                barotropic_names,
            )

            # Define the barotropic pair before any of it is written, rather than on the first day
            # inside the loop. Defining a variable after data has been written to the file only
            # works because the output happens to be NetCDF-4, and relying on that is a trap for
            # whoever changes the format; doing it here costs one extra open of a file the loop
            # opens anyway.
            if !isempty(barotropic_names)
                url = join((config.barotropic_opendap_url, month_directory, first(files)), "/")
                NCDataset(url) do barotropic
                    define_extra_parameters!(output, barotropic, barotropic_subset)
                end
            end

            time_start = 1
            for (file, ds) in zip(files, datasets)
                write_time_dependent_coordinates!(output, ds, full_subset, time_start)
                for name in full_names
                    write_parameter_chunk!(output, ds, name, full_subset, time_start)
                end

                if !isempty(barotropic_names)
                    url = join((config.barotropic_opendap_url, month_directory, file), "/")
                    NCDataset(url) do barotropic
                        time_length(barotropic) == time_length(ds) || error(
                            "Norkyst-v3 sdepth file $url has $(time_length(barotropic)) records " *
                            "against the matching zdepth file's $(time_length(ds)). The two " *
                            "collections are written from one run and are expected to agree.",
                        )

                        for name in barotropic_names
                            write_parameter_chunk!(
                                output,
                                barotropic,
                                name,
                                barotropic_subset,
                                time_start,
                            )
                        end
                    end
                end

                time_start += time_length(ds)
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
