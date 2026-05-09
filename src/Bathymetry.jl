module Bathymetry

export GeonorgeBathymetry, prepare_geonorge_bathymetry, write_bathymetry_file

using Downloads
using ZipFile
using Scratch
using ArchGDAL
using NCDatasets
using NumericalEarth
using Oceananigans
using Oceananigans.Architectures: on_architecture
using Oceananigans.Fields: interior
using Oceananigans.Grids: x_domain, y_domain

using NumericalEarth.DataWrangling: AbstractStaticBathymetry, Metadatum, metadata_path

import NumericalEarth.DataWrangling:
    dataset_variable_name,
    default_download_directory,
    download_dataset,
    latitude_interfaces,
    longitude_interfaces,
    metadata_filename,
    reversed_vertical_axis

const GEONORGE_DYBDEDATA_URL =
    "https://nedlasting.geonorge.no/geonorge/Basisdata/DybdedataKurverGeneraliserte/Shape/" *
    "Basisdata_0000_Norge_25833_DybdedataKurverGeneraliserte_Shape.zip"
const GEONORGE_DYBDEDATA_ZIP = "DybdedataKurverGeneraliserte_Shape.zip"
# NumericalEarth bathymetry regridding constructs a native grid with halo = (10, 10, 1).
# Keep the generated raw dataset comfortably larger than that minimum.
const MIN_NATIVE_BATHYMETRY_SIZE = 24

download_bathymetry_cache::String = ""

function __init__()
    global download_bathymetry_cache = @get_scratch!("Bathymetry")
end

"""
    GeonorgeBathymetry

Static NumericalEarth-compatible metadata wrapper for a regional raw bathymetry
NetCDF derived from Geonorge generalized depth contours.

# Fields
- `metadata_filename`: Filename of the generated raw NetCDF.
- `default_download_directory`: Directory containing the raw NetCDF.
- `longitude_interfaces`: Native longitude bounds of the raw dataset.
- `latitude_interfaces`: Native latitude bounds of the raw dataset.
- `size`: Native `(Nx, Ny, Nz)` size consumed by `NumericalEarth.regrid_bathymetry`.
"""
struct GeonorgeBathymetry <: AbstractStaticBathymetry
    metadata_filename::String
    default_download_directory::String
    longitude_interfaces::NTuple{2, Float64}
    latitude_interfaces::NTuple{2, Float64}
    size::NTuple{3, Int}
end

const GeonorgeBathymetryMetadatum = Metadatum{<:GeonorgeBathymetry}

default_download_directory(dataset::GeonorgeBathymetry) = dataset.default_download_directory
metadata_filename(dataset::GeonorgeBathymetry, args...) = dataset.metadata_filename
longitude_interfaces(dataset::GeonorgeBathymetry) = dataset.longitude_interfaces
latitude_interfaces(dataset::GeonorgeBathymetry) = dataset.latitude_interfaces
reversed_vertical_axis(::GeonorgeBathymetry) = false
Base.size(dataset::GeonorgeBathymetry) = dataset.size

dataset_variable_name(::GeonorgeBathymetryMetadatum) = "z"

"""
    download_dataset(metadata::GeonorgeBathymetryMetadatum)

Return the generated regional raw bathymetry file.
This dataset is materialized locally by `prepare_geonorge_bathymetry`.
"""
function download_dataset(metadata::GeonorgeBathymetryMetadatum)
    filepath = metadata_path(metadata)
    isfile(filepath) || error("Raw bathymetry file $filepath does not exist. Run prepare_geonorge_bathymetry first.")
    return filepath
end

"""
    prepare_geonorge_bathymetry(target_grid; output_path, raw_dir=download_bathymetry_cache,
                                raw_resolution_factor=4, padding_cells=2, regrid_kw...)

Download open Geonorge generalized depth contours, build a regional
NumericalEarth-style raw bathymetry dataset in scratch storage, regrid it onto
`target_grid` with `NumericalEarth.regrid_bathymetry`, and write a processed
NetCDF file compatible with `FjordSim.Grids.ImmersedBoundaryGrid`.

# Keyword arguments
- `output_path`: Destination for the processed FjordSim bathymetry NetCDF.
- `raw_dir`: Scratch directory used for the downloaded archive and regional raw dataset.
- `raw_resolution_factor`: Native raw-grid refinement relative to `target_grid`.
- `padding_cells`: Number of target-grid cell widths added around the requested region.
- `regrid_kw...`: Forwarded directly to `NumericalEarth.regrid_bathymetry`.

# Returns
A named tuple with `dataset`, `raw_path`, `output_path`, and `bottom_height`.
"""
function prepare_geonorge_bathymetry(
    target_grid;
    output_path::String,
    raw_dir::String = download_bathymetry_cache,
    raw_resolution_factor::Int = 4,
    padding_cells::Int = 2,
    regrid_kw...,
)
    raw_resolution_factor >= 1 || throw(ArgumentError("raw_resolution_factor must be >= 1"))
    padding_cells >= 0 || throw(ArgumentError("padding_cells must be >= 0"))

    dataset = geonorge_dataset(
        target_grid;
        raw_dir,
        raw_resolution_factor,
        padding_cells,
    )

    metadata = Metadatum(:bottom_height; dataset)
    bottom_height = NumericalEarth.regrid_bathymetry(target_grid, metadata; regrid_kw...)
    write_bathymetry_file(output_path, target_grid, bottom_height)

    return (; dataset, raw_path = metadata_path(metadata), output_path, bottom_height)
end

"""
    write_bathymetry_file(filepath, target_grid, bottom_height)

Write a processed NetCDF bathymetry file compatible with
`FjordSim.Grids.ImmersedBoundaryGrid`.

# Arguments
- `filepath`: Output NetCDF path.
- `target_grid`: Oceananigans target grid.
- `bottom_height`: Bottom-height field already regridded onto `target_grid`.
"""
function write_bathymetry_file(filepath::String, target_grid, bottom_height)
    isdir(dirname(filepath)) || mkpath(dirname(filepath))

    Nx, Ny, _ = size(target_grid)
    longitude = center_coordinates(x_domain(target_grid), Nx)
    latitude = center_coordinates(y_domain(target_grid), Ny)
    z_faces = vertical_faces(target_grid)

    cpu_bottom_height = on_architecture(CPU(), bottom_height)
    h = Array(interior(cpu_bottom_height, :, :, 1))

    isfile(filepath) && rm(filepath; force = true)

    ds = NCDataset(filepath, "c")
    try
        defDim(ds, "lon", Nx)
        defDim(ds, "lat", Ny)
        defDim(ds, "zf", length(z_faces))

        lon = defVar(ds, "lon", Float64, ("lon",))
        lat = defVar(ds, "lat", Float64, ("lat",))
        zf = defVar(ds, "z_faces", Float64, ("zf",))
        hvar = defVar(ds, "h", Float32, ("lon", "lat"))

        lon[:] = longitude
        lat[:] = latitude
        zf[:] = z_faces
        hvar[:, :] = h
    finally
        close(ds)
    end

    return filepath
end

"""
    geonorge_dataset(target_grid; raw_dir, raw_resolution_factor, padding_cells)

Construct the regional raw bathymetry dataset wrapper used as input to
`NumericalEarth.regrid_bathymetry`.
"""
function geonorge_dataset(target_grid; raw_dir, raw_resolution_factor, padding_cells)
    isdir(raw_dir) || mkpath(raw_dir)

    Nx, Ny, _ = size(target_grid)
    longitude = expand_domain(x_domain(target_grid), Nx, padding_cells)
    latitude = expand_domain(y_domain(target_grid), Ny, padding_cells)
    raw_size = (
        max(raw_resolution_factor * (Nx + 2 * padding_cells), MIN_NATIVE_BATHYMETRY_SIZE),
        max(raw_resolution_factor * (Ny + 2 * padding_cells), MIN_NATIVE_BATHYMETRY_SIZE),
        1,
    )

    raw_filename = geonorge_raw_filename(longitude, latitude, raw_size)
    raw_path = joinpath(raw_dir, raw_filename)

    if !isfile(raw_path)
        shapefiles = ensure_geonorge_contours!(raw_dir)
        write_native_bathymetry(raw_path, shapefiles; longitude, latitude, size = raw_size)
    end

    return GeonorgeBathymetry(raw_filename, raw_dir, longitude, latitude, raw_size)
end

"""
    ensure_geonorge_contours!(raw_dir)

Download and extract the Geonorge generalized contour shapefiles into `raw_dir`
if they are not already present. Returns sorted `.shp` paths.
"""
function ensure_geonorge_contours!(raw_dir)
    zip_path = joinpath(raw_dir, GEONORGE_DYBDEDATA_ZIP)
    if !isfile(zip_path)
        @info "Downloading Geonorge generalized depth contours to $raw_dir..."
        Downloads.download(GEONORGE_DYBDEDATA_URL, zip_path)
    end

    shapefiles = String[]
    reader = ZipFile.Reader(zip_path)
    try
        for file in reader.files
            keep = endswith(file.name, ".shp") || endswith(file.name, ".shx") || endswith(file.name, ".dbf") || endswith(file.name, ".prj")
            keep || continue

            destination = joinpath(raw_dir, basename(file.name))
            if !isfile(destination)
                open(destination, "w") do io
                    write(io, read(file))
                end
            end

            endswith(destination, ".shp") && push!(shapefiles, destination)
        end
    finally
        close(reader)
    end

    isempty(shapefiles) && error("No Geonorge shapefiles were extracted from $zip_path.")

    return sort(shapefiles)
end

"""
    write_native_bathymetry(filepath, shapefiles; longitude, latitude, size)

Create the regional raw bathymetry NetCDF consumed by `GeonorgeBathymetry`.
The native raster is built by sampling contour vertices, transforming them to
WGS84, and gridding them with GDAL's inverse-distance interpolator.
"""
function write_native_bathymetry(filepath, shapefiles; longitude, latitude, size)
    Nx, Ny, _ = size
    longitude_centers = center_coordinates(longitude, Nx)
    latitude_centers = center_coordinates(latitude, Ny)
    z_data = build_native_bathymetry_data(shapefiles; longitude, latitude, Nx, Ny)

    isfile(filepath) && rm(filepath; force = true)

    ds = NCDataset(filepath, "c")
    try
        defDim(ds, "lon", Nx)
        defDim(ds, "lat", Ny)

        lon = defVar(ds, "lon", Float64, ("lon",))
        lat = defVar(ds, "lat", Float64, ("lat",))
        z = defVar(ds, "z", Float32, ("lon", "lat"))

        lon[:] = longitude_centers
        lat[:] = latitude_centers
        z[:, :] = z_data
    finally
        close(ds)
    end

    return filepath
end

"""
    build_native_bathymetry_data(shapefiles; longitude, latitude, Nx, Ny)

Create the raw regional bathymetry array by sampling contour vertices and
gridding them onto a regular WGS84 longitude-latitude raster.
"""
function build_native_bathymetry_data(shapefiles; longitude, latitude, Nx, Ny)
    filter_bounds = transformed_filter_bounds(longitude, latitude)

    return ArchGDAL.importEPSG(25833; order = :trad) do source_srs
        ArchGDAL.importEPSG(4326; order = :trad) do target_srs
            ArchGDAL.createcoordtrans(source_srs, target_srs) do transform
                create_point_dataset(shapefiles, transform, filter_bounds, target_srs) do point_dataset
                    grid_point_dataset(point_dataset; longitude, latitude, Nx, Ny)
                end
            end
        end
    end
end

"""
    create_point_dataset(shapefiles, transform, filter_bounds, target_srs, f)

Create an in-memory point dataset in WGS84 containing sampled contour vertices,
then invoke `f(point_dataset)`.
"""
function create_point_dataset(f::Function, shapefiles, transform, filter_bounds, target_srs)
    ArchGDAL.create(ArchGDAL.getdriver("Memory")) do point_dataset
        ArchGDAL.createlayer(
            name = "bathymetry_points",
            dataset = point_dataset,
            geom = ArchGDAL.wkbPoint,
            spatialref = target_srs,
        ) do point_layer
            ArchGDAL.addfielddefn!(point_layer, "z", ArchGDAL.OFTReal)
            point_count = sample_contours!(point_layer, shapefiles, transform, filter_bounds)
            point_count > 0 || error("No Geonorge depth contours intersect the requested region.")
            return f(point_dataset)
        end
    end
end

"""
    grid_point_dataset(point_dataset; longitude, latitude, Nx, Ny)

Grid the sampled point dataset onto a regular native bathymetry raster.
"""
function grid_point_dataset(point_dataset; longitude, latitude, Nx, Ny)
    options = [
        "-of", "MEM",
        "-a", "invdistnn:power=2:smoothing=0.2:max_points=16:min_points=1:nodata=0",
        "-zfield", "z",
        "-txe", string(longitude[1]), string(longitude[2]),
        "-tye", string(latitude[1]), string(latitude[2]),
        "-outsize", string(Nx), string(Ny),
        "-a_srs", "EPSG:4326",
        "-l", "bathymetry_points",
    ]

    ArchGDAL.gdalgrid(point_dataset, options) do raster_dataset
        band = ArchGDAL.getband(raster_dataset, 1)
        raster_to_bottom_height(band, Nx, Ny)
    end
end

"""
    sample_contours!(point_layer, shapefiles, transform, filter_bounds)

Sample contour vertices from all shapefiles intersecting `filter_bounds`,
transform them to WGS84, and append them as point features with a `z` value.
Returns the number of sampled points.
"""
function sample_contours!(point_layer, shapefiles, transform, filter_bounds)
    xmin, ymin, xmax, ymax = filter_bounds
    point_count = 0

    for shapefile in shapefiles
        ArchGDAL.read(shapefile) do dataset
            layer = ArchGDAL.getlayer(dataset, 0)
            ArchGDAL.setspatialfilter!(layer, xmin, ymin, xmax, ymax)
            depth_index = ArchGDAL.findfieldindex(layer, "DYBDE", false)

            for feature in layer
                depth = ArchGDAL.getfield(feature, depth_index)
                ismissing(depth) && continue

                line = ArchGDAL.getgeom(feature)
                point_count += add_linestring_points!(point_layer, line, transform, -abs(Float64(depth)))
            end
        end
    end

    return point_count
end

"""
    add_linestring_points!(point_layer, line, transform, bottom_height)

Append each vertex of a contour line to `point_layer` with a constant
`bottom_height` attribute.
"""
function add_linestring_points!(point_layer, line, transform, bottom_height)
    npoints = ArchGDAL.ngeom(line)
    added = 0

    for point_index in 0:npoints-1
        x, y, _ = ArchGDAL.getpoint(line, point_index)
        point = ArchGDAL.createpoint(x, y)
        ArchGDAL.transform!(point, transform)
        longitude, latitude, _ = ArchGDAL.getpoint(point, 0)

        ArchGDAL.createfeature(point_layer) do feature
            ArchGDAL.setfield!(feature, ArchGDAL.findfieldindex(feature, "z"), bottom_height)
            ArchGDAL.setgeom!(feature, 0, ArchGDAL.createpoint(longitude, latitude))
            return nothing
        end

        added += 1
    end

    return added
end

"""
    raster_to_bottom_height(band, Nx, Ny)

Read a GDAL raster band into the `(Nx, Ny)` orientation used by FjordSim's
bathymetry NetCDF files.
"""
function raster_to_bottom_height(band, Nx, Ny)
    data = Array(ArchGDAL.read(band))
    if size(data) == (Ny, Nx)
        data = permutedims(data, (2, 1))
    elseif size(data) != (Nx, Ny)
        error("Unexpected raster shape $(size(data)); expected ($Nx, $Ny) or ($Ny, $Nx).")
    end

    return reverse(Float32.(data), dims = 2)
end

"""
    transformed_filter_bounds(longitude, latitude)

Transform a WGS84 longitude/latitude bounding box into EPSG:25833 bounds for
spatial filtering of the raw Geonorge contour shapefiles.
"""
function transformed_filter_bounds(longitude, latitude)
    corners = (
        (longitude[1], latitude[1]),
        (longitude[1], latitude[2]),
        (longitude[2], latitude[1]),
        (longitude[2], latitude[2]),
    )

    xs = Float64[]
    ys = Float64[]

    ArchGDAL.importEPSG(4326; order = :trad) do source_srs
        ArchGDAL.importEPSG(25833; order = :trad) do target_srs
            ArchGDAL.createcoordtrans(source_srs, target_srs) do transform
                for (lon, lat) in corners
                    point = ArchGDAL.createpoint(lon, lat)
                    ArchGDAL.transform!(point, transform)
                    x, y, _ = ArchGDAL.getpoint(point, 0)
                    push!(xs, x)
                    push!(ys, y)
                end
            end
        end
    end

    return minimum(xs), minimum(ys), maximum(xs), maximum(ys)
end

"""
    center_coordinates(domain, N)

Return `N` uniformly spaced cell-center coordinates for the interval
`domain = (lower, upper)`.
"""
function center_coordinates(domain, N)
    delta = domain_step(domain, N)
    return collect(range(domain[1] + delta / 2, step = delta, length = N))
end

domain_step(domain, N) = (domain[2] - domain[1]) / N

"""
    vertical_faces(grid)

Extract the physical vertical face coordinates from a live Oceananigans grid,
excluding halo entries.
"""
function vertical_faces(grid)
    z_faces_with_halo = getproperty(grid.z, Symbol("cᵃᵃᶠ"))
    first_index = grid.Hz + 1
    last_index = first_index + grid.Nz
    return collect(z_faces_with_halo[first_index:last_index])
end

"""
    expand_domain(domain, N, padding_cells)

Expand a uniformly spaced domain by `padding_cells` cell widths on both sides.
"""
function expand_domain(domain, N, padding_cells)
    delta = domain_step(domain, N)
    lower = domain[1] - padding_cells * delta
    upper = domain[2] + padding_cells * delta
    return (lower, upper)
end

"""
    geonorge_raw_filename(longitude, latitude, size)

Build a deterministic filename for a regional raw bathymetry cache file.
"""
function geonorge_raw_filename(longitude, latitude, size)
    Nx, Ny, _ = size
    lon1 = replace(string(round(longitude[1], digits = 3)), '.' => 'p')
    lon2 = replace(string(round(longitude[2], digits = 3)), '.' => 'p')
    lat1 = replace(string(round(latitude[1], digits = 3)), '.' => 'p')
    lat2 = replace(string(round(latitude[2], digits = 3)), '.' => 'p')

    return "geonorge_bathymetry_$(lon1)_$(lon2)_$(lat1)_$(lat2)_$(Nx)x$(Ny).nc"
end

end  # module Bathymetry