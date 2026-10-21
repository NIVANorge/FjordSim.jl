module NorKyst

# The NorKyst-800m dataset itself, shared by the forcing, initial-state and open-boundary adapters:
# where its collections live, what its variables are called, its source grid, and how a subset of
# it is listed, downloaded and written. Nothing here names a pipeline — a config type and the hooks
# it implements belong to the adapter in that pipeline's directory.

using Downloads
using NCDatasets

using ..Utils: with_retries
using ..Regridding: ProjectedSourceGrid

# The operational daily means (`..._avg.an.*`, one record a file, from 2017-02-20).
const NORKYST_CATALOG_URL = "https://thredds.met.no/thredds/catalog/fou-hi/norkyst800m/catalog.xml"
const NORKYST_OPENDAP_URL = "https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m/"

# The operational hourly instants (`..._his.an.*`, 24 records a file), which keep the tide.
const NORKYST_HOURLY_CATALOG_URL = "https://thredds.met.no/thredds/catalog/fou-hi/norkyst800m-1h/catalog.xml"
const NORKYST_HOURLY_OPENDAP_URL = "https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h/"

# The Norkyst-v3 hindcast (from 2012), one file a day under `YYYY/MM/`. `zdepth` carries the
# z-interpolated fields; `sdepth` the same run on its native levels, and the only one carrying the
# barotropic velocity pair.
const NORKYST_HINDCAST_CATALOG_URL =
    "https://thredds.met.no/thredds/catalog/romshindcast/norkyst_v3/zdepth"
const NORKYST_HINDCAST_OPENDAP_URL =
    "https://thredds.met.no/thredds/dodsC/romshindcast/norkyst_v3/zdepth"
const NORKYST_HINDCAST_BAROTROPIC_CATALOG_URL =
    "https://thredds.met.no/thredds/catalog/romshindcast/norkyst_v3/sdepth"
const NORKYST_HINDCAST_BAROTROPIC_OPENDAP_URL =
    "https://thredds.met.no/thredds/dodsC/romshindcast/norkyst_v3/sdepth"

# The daily-mean NorKyst variables and the FjordSim names they become, for the forcing and the
# initial state. Only the intersection of this mapping with a config's `parameters` is prepared.
const NORKYST_VARIABLE_NAMES = Dict(
    "temperature" => "T",
    "salinity" => "S",
    "u_eastward" => "u",
    "v_northward" => "v",
)
# thredds.met.no answers `429 Too Many Requests` to bursts: opening a month's daily files back to
# back, at three requests an open, tripped it every five or so opens. Requests are spaced by this.
const NORKYST_REQUEST_INTERVAL_SECONDS = 0.5
const NORKYST_OPEN_REQUESTS = 3
const NORKYST_NEXT_REQUEST = Ref(0.0)

# Scalar NorKyst variable whose `proj4` attribute defines the projected X/Y coordinates.
const NORKYST_PROJECTION_VARIABLE = "projection_stere"

"""
    norkyst_source_grid(filepath)

Read the projected coordinates, depth levels and projection of a downloaded NorKyst subset.
Errors unless the projected coordinates are regularly spaced, which `source_field_grid` needs in
order to express them as a `RectilinearGrid`. Every NorKyst collection shares this grid, its
projection variable and its depth axis.
"""
function norkyst_source_grid(filepath)
    return NCDataset(filepath) do ds
        x = Array{Float64}(ds["X"][:])
        y = Array{Float64}(ds["Y"][:])
        depths = Array{Float64}(ds["depth"][:])
        proj4 = NCDatasets.variable(ds, NORKYST_PROJECTION_VARIABLE).attrib["proj4"]

        length(x) >= 2 && length(y) >= 2 ||
            error(
                "NorKyst subset in $filepath is too small to interpolate from: " *
                "$(length(x))x$(length(y)). For an open-boundary download, raise the boundary " *
                "config's `margin`.",
            )
        all(difference -> isapprox(difference, x[2] - x[1]), diff(x)) &&
            all(difference -> isapprox(difference, y[2] - y[1]), diff(y)) ||
            error("NorKyst projected coordinates in $filepath are not regularly spaced.")

        return ProjectedSourceGrid(x, y, depths, proj4)
    end
end

function list_opendap_files(; catalog_url)
    catalog = with_retries(() -> read(Downloads.download(catalog_url), String), catalog_url)

    files = String[]
    for match in eachmatch(r"<dataset\b[^>]*\bname=\"([^\"]+\.nc)\""i, catalog)
        push!(files, match.captures[1])
    end

    return files
end


"""
    list_hindcast_files(catalog_url, year, month)

The daily file names the hindcast publishes for one month, sorted. Walks that month's own THREDDS
sub-catalog, since the hindcast has no flat catalog to filter. An absent month returns empty rather
than raising: the collection's ends are ragged, and the caller reports the gap.
"""
function list_hindcast_files(catalog_url, year, month)
    url = join((catalog_url, string(year), lpad(month, 2, '0'), "catalog.xml"), "/")

    catalog = try
        with_retries(() -> read(Downloads.download(url), String), url)
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
    open_remote(url)
    open_remote(body, url)

Open the OPeNDAP dataset at `url`, retrying as `with_retries` does. With `body`, run it on the
dataset and close it; a failure inside `body` retries it too, from a fresh open.
"""
open_remote(url) = with_retries(() -> paced(() -> NCDataset(url), NORKYST_OPEN_REQUESTS), url)
open_remote(body, url) =
    with_retries(() -> paced(() -> NCDataset(body, url), NORKYST_OPEN_REQUESTS), url)

"""
    read_remote(variable, indices, source)

The `indices` of `variable` from the OPeNDAP dataset `source`, paced and retried.
"""
read_remote(variable, indices, source) =
    with_retries(() -> paced(() -> variable[indices...], 1), "$(name(variable)) in $(path(source))")

"""
    paced(body, requests)

Run `body()`, which makes `requests` HTTP requests, no sooner than the previous paced call allows,
and hold the next one back by `requests` times `NORKYST_REQUEST_INTERVAL_SECONDS`.
"""
function paced(body, requests)
    wait_seconds = NORKYST_NEXT_REQUEST[] - time()
    wait_seconds > 0 && sleep(wait_seconds)
    NORKYST_NEXT_REQUEST[] = time() + requests * NORKYST_REQUEST_INTERVAL_SECONDS

    return body()
end

function bounding_range(mask, dimension)
    reduced_dims = Tuple(i for i in 1:ndims(mask) if i != dimension)
    axis_mask = vec(any(mask; dims = reduced_dims))
    indices = findall(axis_mask)
    isempty(indices) && error("No NorKyst points found inside the requested lon/lat range.")
    return first(indices):last(indices)
end

"""
    NorKystSubset

The spatial window of a NorKyst dataset selected by a setup, plus the variables to extract.
Derived from a longitude/latitude box by `subset_ranges` and passed to every function that writes
the subset out.

# Fields
- `ranges`: Index range per dataset dimension covering the requested lon/lat box.
- `mask`: In-box flag per point of the subset, used to blank out points outside the box.
- `spatial_dimensions`: Names of the dimensions `mask` is indexed by, in order.
- `parameters`: Variable names to extract.
"""
struct NorKystSubset{M<:AbstractArray,S<:Tuple}
    ranges::Dict{String,UnitRange{Int}}
    mask::M
    spatial_dimensions::S
    parameters::Vector{String}
end

"""
    subset_ranges(ds, longitude, latitude, parameters)

The NorKyst index window of `ds` covering the `(minimum, maximum)` `longitude` and `latitude` box,
extracting `parameters`.
"""
function subset_ranges(ds, longitude, latitude, parameters)
    latitude_variable = variable(ds, "lat")
    longitude_variable = variable(ds, "lon")
    source_latitude = Array(latitude_variable[ntuple(_ -> :, ndims(latitude_variable))...])
    source_longitude = Array(longitude_variable[ntuple(_ -> :, ndims(longitude_variable))...])
    mask = (source_latitude .>= latitude[1]) .&
           (source_latitude .<= latitude[2]) .&
           (source_longitude .>= longitude[1]) .&
           (source_longitude .<= longitude[2])

    ranges = Dict{String,UnitRange{Int}}()
    for (i, dimension) in enumerate(dimnames(latitude_variable))
        ranges[dimension] = bounding_range(mask, i)
    end

    spatial_dimensions = dimnames(latitude_variable)
    subset_mask = mask[(ranges[dimension] for dimension in spatial_dimensions)...]
    return NorKystSubset(ranges, subset_mask, spatial_dimensions, parameters)
end

function variable_indices(variable, ranges)
    return ntuple(ndims(variable)) do i
        dimension = dimnames(variable)[i]
        get(ranges, dimension, :)
    end
end

function time_dimension(ds)
    time_dimensions = [name for name in dimnames(ds) if occursin("time", lowercase(name))]
    !isempty(time_dimensions) && return first(time_dimensions)
    haskey(ds, "time") && return first(dimnames(ds["time"]))
    error("Could not find a time dimension in NorKyst dataset.")
end

function time_length(ds)
    return NCDatasets.dim(ds, time_dimension(ds))
end

function copy_attributes!(destination, source)
    for (key, value) in source.attrib
        key == "_FillValue" && continue
        destination.attrib[key] = value
    end
end

filtered_attributes(source) = [(key, value) for (key, value) in source.attrib if key != "_FillValue"]
decoded_attributes(source) = [
    (key, value) for (key, value) in source.attrib
    if key ∉ ("_FillValue", "missing_value", "scale_factor", "add_offset")
]

function concrete_float_data(data)
    T = nonmissingtype(eltype(data))
    T <: AbstractFloat || error("Expected floating-point NorKyst data, got $(eltype(data)).")
    output = Array{T}(undef, size(data))

    for index in eachindex(data)
        value = data[index]
        output[index] = ismissing(value) ? convert(T, NaN) : value
    end

    return output
end

function define_subset_variable!(output, source, name, subset::NorKystSubset; deflatelevel = 5)
    variable = NCDatasets.variable(source, name)
    dimensions = dimnames(variable)
    indices = variable_indices(variable, subset.ranges)
    data = variable[indices...]
    variable_type = data isa AbstractArray ? eltype(data) : typeof(data)

    output_variable = defVar(
        output,
        name,
        variable_type,
        dimensions;
        deflatelevel,
        attrib = filtered_attributes(variable),
    )

    if time_dimension(source) ∉ dimensions
        if ndims(output_variable) == 0
            output_variable[] = data
        else
            output_variable[ntuple(_ -> :, ndims(output_variable))...] = data
        end
    end

    return output_variable
end

function copy_auxiliary_variable(name, variable, time_dimension_name, subset::NorKystSubset)
    name in subset.parameters && return false
    dimensions = dimnames(variable)
    return name == time_dimension_name || time_dimension_name ∉ dimensions
end

function define_output_file(output_path, template, subset::NorKystSubset, total_time)
    time_dimension_name = time_dimension(template)
    isfile(output_path) && rm(output_path; force = true)

    output = NCDataset(output_path, "c")
    try
        for dimension in dimnames(template)
            dimension_length = if dimension == time_dimension_name
                total_time
            elseif haskey(subset.ranges, dimension)
                length(subset.ranges[dimension])
            else
                NCDatasets.dim(template, dimension)
            end
            defDim(output, dimension, dimension_length)
        end

        for name in keys(template)
            variable = NCDatasets.variable(template, name)
            copy_auxiliary_variable(name, variable, time_dimension_name, subset) || continue
            all(dimension -> dimension in dimnames(output), dimnames(variable)) || continue
            define_subset_variable!(output, template, name, subset; deflatelevel = 0)
        end

        for name in subset.parameters
            variable = template[name]
            variable_type = nonmissingtype(eltype(variable))
            defVar(
                output,
                name,
                variable_type,
                dimnames(variable);
                deflatelevel = 5,
                attrib = decoded_attributes(variable),
            )
        end

        copy_attributes!(output, template)
    catch
        close(output)
        rethrow()
    end

    return output
end

function masked_fill_value(data)
    return if eltype(data) <: AbstractFloat
        convert(eltype(data), NaN)
    else
        error("Cannot mask data with element type $(eltype(data)).")
    end
end

function apply_spatial_mask(data, variable, subset::NorKystSubset)
    spatial_indices = [findfirst(==(dimension), dimnames(variable)) for dimension in subset.spatial_dimensions]
    any(isnothing, spatial_indices) && return data

    masked_data = copy(data)
    fill_value = masked_fill_value(masked_data)

    for mask_index in CartesianIndices(subset.mask)
        subset.mask[mask_index] && continue
        data_indices = ntuple(ndims(masked_data)) do i
            mask_dimension = findfirst(==(i), spatial_indices)
            mask_dimension === nothing ? (:) : mask_index[mask_dimension]
        end
        masked_data[data_indices...] .= fill_value
    end

    return masked_data
end

function write_parameter_chunk!(output, source, name, subset::NorKystSubset, time_start)
    variable = source[name]
    input_indices = variable_indices(variable, subset.ranges)
    data = concrete_float_data(read_remote(variable, input_indices, source))
    data = apply_spatial_mask(data, variable, subset)

    output_variable = output[name]
    output_indices = ntuple(ndims(output_variable)) do i
        dimension = dimnames(output_variable)[i]
        if dimension == time_dimension(source)
            time_start:(time_start + size(data, i) - 1)
        else
            (:)
        end
    end

    output_variable[output_indices...] = data
    return size(data, findfirst(==(time_dimension(source)), dimnames(variable)))
end

function write_time_dependent_coordinates!(output, source, subset::NorKystSubset, time_start)
    time_dimension_name = time_dimension(source)

    for name in keys(source)
        name in subset.parameters && continue
        variable = NCDatasets.variable(source, name)
        time_index = findfirst(==(time_dimension_name), dimnames(variable))
        time_index === nothing && continue
        haskey(output, name) || continue

        input_indices = variable_indices(variable, subset.ranges)
        data = read_remote(variable, input_indices, source)
        output_variable = NCDatasets.variable(output, name)
        output_indices = ntuple(ndims(output_variable)) do i
            if i == time_index
                time_start:(time_start + size(data, i) - 1)
            else
                (:)
            end
        end
        output_variable[output_indices...] = data
    end
end

end # module
