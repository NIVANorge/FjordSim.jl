module Topaz

# The TOPAZ5-ECOSMO dataset (Copernicus Marine ARCTIC_ANALYSISFORECAST_BGC_002_004, MET Norway),
# shared by the initial-state and open-boundary adapters: where its daily files live, how a day's
# file is chosen, its source grid, and how a subset of it is downloaded. Nothing here names a
# pipeline.
#
# HYCOM with ECOSMO II on a 6.25 km polar stereographic grid, daily means on 40 fixed depths. The
# product is described as covering north of 63°N, but the files carry the whole model grid: off
# Rogaland (59°N) its water reaches 266 m.

using Dates: Date, DateTime, Day, @dateformat_str
import Dates
using Downloads
using NCDatasets

using ..Utils: with_retries
using ..Regridding: ProjectedSourceGrid

const TOPAZ_CATALOG_URL = "https://thredds.met.no/thredds/catalog/cmems/topaz5_bgc_dm_files"
const TOPAZ_OPENDAP_URL = "https://thredds.met.no/thredds/dodsC/cmems/topaz5_bgc_dm_files"

# The published biogeochemical variables, under their own names. `dissic` is mol m⁻³, the
# concentrations mmol m⁻³ and `chl` mg m⁻³; the units are kept, so a setup converts them.
const TOPAZ_VARIABLES = ["no3", "po4", "si", "o2", "dissic", "chl", "phyc", "zooc", "ph", "kd", "nppv"]
const TOPAZ_VARIABLE_NAMES = Dict(name => name for name in TOPAZ_VARIABLES)

# The projected coordinates are published in units of 100 km.
const TOPAZ_COORDINATE_SCALE = 1.0e5
const TOPAZ_PROJECTION_VARIABLE = "stereographic"

# thredds.met.no answers bursts with `429 Too Many Requests` (`NorKyst.paced`); opens are spaced.
const TOPAZ_OPEN_INTERVAL_SECONDS = 1.0

"""
    topaz_day_files(catalog_url, year, month)

The file for each day of one month, as `Date => file name`. A day is published several times, once
by each daily run that reached it; the latest run's is kept, since it is an analysis for a day it
has passed and a forecast otherwise. An absent month returns an empty `Dict`.
"""
function topaz_day_files(catalog_url, year, month)
    url = join((catalog_url, string(year), lpad(month, 2, '0'), "catalog.xml"), "/")
    catalog = try
        with_retries(() -> read(Downloads.download(url), String), url)
    catch exception
        exception isa InterruptException && rethrow()
        @warn "Could not read the TOPAZ5 catalog for $year-$(lpad(month, 2, '0'))" url exception
        return Dict{Date,String}()
    end
    return latest_bulletins(catalog)
end

"""
    latest_bulletins(catalog)

`Date => file name` for every day a THREDDS catalog lists, keeping the file of the latest run
(`-b<yyyymmdd>-`) for each day.
"""
function latest_bulletins(catalog)
    files = Dict{Date,Tuple{Date,String}}()
    pattern = r"\bname=\"((\d{8})_dm-[^\"]*-b(\d{8})-[^\"]*\.nc)\""
    for match in eachmatch(pattern, catalog)
        name, day, bulletin = match.captures
        day, bulletin = Date(day, dateformat"yyyymmdd"), Date(bulletin, dateformat"yyyymmdd")
        (!haskey(files, day) || bulletin > first(files[day])) && (files[day] = (bulletin, name))
    end
    return Dict(day => last(entry) for (day, entry) in files)
end

topaz_url(opendap_url, day, name) = join((opendap_url, Dates.format(day, dateformat"yyyy/mm"), name), "/")

"""
    topaz_source_grid(filepath)

The projected coordinates (in metres), depths and projection of a downloaded TOPAZ subset.
"""
function topaz_source_grid(filepath)
    return NCDataset(filepath) do ds
        x = Array{Float64}(ds["x"][:])
        y = Array{Float64}(ds["y"][:])
        depths = Array{Float64}(ds["depth"][:])
        proj4 = NCDatasets.variable(ds, TOPAZ_PROJECTION_VARIABLE).attrib["proj4"]
        length(x) >= 2 && length(y) >= 2 || error(
            "TOPAZ subset in $filepath is too small to interpolate from: $(length(x))x$(length(y)). " *
            "Raise the config's `margin`.",
        )
        return ProjectedSourceGrid(x, y, depths, proj4)
    end
end

"""
    topaz_window(ds, longitude, latitude)

The `x` and `y` index ranges of the TOPAZ grid points inside the `(minimum, maximum)` `longitude`
and `latitude` box.
"""
function topaz_window(ds, longitude, latitude)
    λ = ds["longitude"][:, :]
    φ = ds["latitude"][:, :]
    inside = (longitude[1] .<= λ .<= longitude[2]) .& (latitude[1] .<= φ .<= latitude[2])
    any(inside) || error("No TOPAZ points inside longitude $longitude, latitude $latitude.")
    span(dimension) = findall(vec(any(inside; dims = 3 - dimension)))
    return first(span(1)):last(span(1)), first(span(2)):last(span(2))
end

"""
    download_topaz(output_path, days, longitude, latitude, parameters; catalog_url, opendap_url)

Download `parameters` for every day of `days` that TOPAZ publishes, subset to the `longitude` and
`latitude` box, into one NetCDF at `output_path` on `(x, y, depth, time)`: `x` and `y` in metres,
land `NaN`, one record per day stamped at the day's 12:00. Written under a temporary name and
renamed once complete. Returns the path, or `nothing` when TOPAZ has none of the days.
"""
function download_topaz(output_path, days, longitude, latitude, parameters; catalog_url, opendap_url)
    months = unique((Dates.year(day), Dates.month(day)) for day in days)
    published = merge((topaz_day_files(catalog_url, year, month) for (year, month) in months)...)
    available = [day for day in days if haskey(published, day)]
    missing_days = setdiff(days, available)
    isempty(missing_days) || @warn "TOPAZ5 has no file for $(join(missing_days, ", "))"
    isempty(available) && return nothing

    temporary_path = output_path * ".tmp"
    isfile(temporary_path) && rm(temporary_path)
    first_url = topaz_url(opendap_url, first(available), published[first(available)])

    with_retries(first_url) do
        NCDataset(first_url) do template
            x_range, y_range = topaz_window(template, longitude, latitude)
            NCDataset(temporary_path, "c") do output
                define_topaz_file!(output, template, x_range, y_range, length(available), parameters)
                for (index, day) in enumerate(available)
                    sleep(TOPAZ_OPEN_INTERVAL_SECONDS)
                    url = topaz_url(opendap_url, day, published[day])
                    with_retries(url) do
                        NCDataset(url) do source
                            output["time"][index] = DateTime(day) + Dates.Hour(12)
                            for name in parameters
                                values = source[name][x_range, y_range, :, 1]
                                output[name][:, :, :, index] = Float32.(coalesce.(values, NaN32))
                            end
                        end
                    end
                end
            end
        end
    end

    mv(temporary_path, output_path; force = true)
    return output_path
end

function define_topaz_file!(output, template, x_range, y_range, records, parameters)
    defDim(output, "x", length(x_range))
    defDim(output, "y", length(y_range))
    defDim(output, "depth", NCDatasets.dim(template, "depth"))
    defDim(output, "time", records)

    defVar(output, "x", Float64.(template["x"][x_range]) .* TOPAZ_COORDINATE_SCALE, ("x",);
        attrib = ["units" => "m", "standard_name" => "projection_x_coordinate"])
    defVar(output, "y", Float64.(template["y"][y_range]) .* TOPAZ_COORDINATE_SCALE, ("y",);
        attrib = ["units" => "m", "standard_name" => "projection_y_coordinate"])
    defVar(output, "depth", Float64.(template["depth"][:]), ("depth",);
        attrib = ["units" => "m", "positive" => "down"])
    defVar(output, "longitude", Float64.(template["longitude"][x_range, y_range]), ("x", "y"))
    defVar(output, "latitude", Float64.(template["latitude"][x_range, y_range]), ("x", "y"))
    defVar(output, "time", Float64, ("time",); attrib = ["units" => "hours since 1950-01-01 00:00:00"])
    defVar(output, TOPAZ_PROJECTION_VARIABLE, Int32, ();
        attrib = ["proj4" => template[TOPAZ_PROJECTION_VARIABLE].attrib["proj4"]])

    for name in parameters
        attributes = ["units" => get(template[name].attrib, "units", ""), "_FillValue" => NaN32]
        defVar(output, name, Float32, ("x", "y", "depth", "time"); deflatelevel = 5, attrib = attributes)
    end
    return output
end

end # module
