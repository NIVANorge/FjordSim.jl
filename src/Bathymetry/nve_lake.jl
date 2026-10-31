# NVE lake depth-map adapter: the bathymetry of one lake from NVE's lake database (Innsjødatabase),
# whose map service serves the surveyed depth points, depth contours and outline of every lake NVE
# has a depth map for. The Geonorge Sjøkart Dybdedata `DybdedataConfig` reads is a sea chart and
# holds no soundings in an inland lake.

using JSON

const NVE_LAKE_MAP_URL = "https://kart.nve.no/enterprise/rest/services/Innsjodatabase2/MapServer"
# The service's layers: `DybdePunkt`, `DybdeKurve` and `Innsjodatabase`, the last the lake polygons
# with their islands as holes.
const NVE_LAKE_LAYERS = (points = 1, contours = 2, outline = 5)
# Depth below the lake surface, in metres, on the point and contour layers.
const NVE_LAKE_DEPTH_FIELD = "dybde_m"
# The layers' own `maxRecordCount`; a full page sets `exceededTransferLimit`.
const NVE_LAKE_PAGE_SIZE = 2000

"""
    NVELakeBathymetryConfig

Configuration for `prepare_bathymetry` from NVE's lake depth maps, for a domain holding one lake.
The lake is named by its NVE lake number, `vatnLnr` in the lake database (Mjøsa is 118), and
`prepare_bathymetry` downloads its depth points, depth contours and outline into
`data_root/download_directory`, reusing them when present.

Depths are below the lake surface, so the lake surface is `z = 0` whatever its height above the sea.
The points and contour vertices are gridded with the outline's vertices at 0 m by the linear
interpolation `DybdedataConfig` uses, and every cell whose centre is outside the outline — the shore
and the lake's islands — is land (`write_lake_bathymetry`).

`data_root`, `output_file`, `plot_file` and `lake_number` are required.

# Fields
- `download_directory`: where the downloaded GeoJSON and the raw gridded file go, under `data_root`.
- `raw_resolution_factor`, `padding_cells`: the raw grid, as on `DybdedataConfig`.
- `interpolation_passes`, `major_basins`, `minimum_depth`, `regrid_cache`: passed to
  `NumericalEarth.regrid_bathymetry`, as on `DybdedataConfig`.
- `max_island_cells`, `close_narrow_passages`, `spike_ratio`, `minimum_cell_fraction`,
  `max_slope_factor`, `preserve_sills`: the `smooth_bathymetry_gaps!` stages, as on
  `DybdedataConfig`. A lake has no open edge, so the two stages acting on one are not offered.

# Fields derived from the target grid by `native_extent!`
- `raw_file`, `longitude`, `latitude`, `size`: the raw gridded file and its region.
"""
Base.@kwdef mutable struct NVELakeBathymetryConfig <: AbstractBathymetryConfig
    data_root::String
    output_file::String
    plot_file::String
    lake_number::Int
    download_directory::String = "nve_lake"
    raw_resolution_factor::Int = 4
    padding_cells::Int = 0
    interpolation_passes::Int = 1
    major_basins::Int = 1
    minimum_depth::Float64 = 0.0
    max_island_cells::Int = 0
    close_narrow_passages::Bool = false
    spike_ratio::Float64 = 0.0
    minimum_cell_fraction::Float64 = 0.0
    max_slope_factor::Float64 = 0.0
    preserve_sills::Bool = false
    regrid_cache::Bool = true
    raw_file::String = ""
    longitude::NTuple{2,Float64} = (0.0, 0.0)
    latitude::NTuple{2,Float64} = (0.0, 0.0)
    size::NTuple{3,Int} = (0, 0, 1)
end

"""
    bathymetry_dataset(target_grid, config::NVELakeBathymetryConfig)

Download the lake's depth map if it is not on disk, grid it onto the raw grid for `target_grid`,
and return the raw file as a `RawBathymetry`.
"""
function bathymetry_dataset(target_grid, config::NVELakeBathymetryConfig)
    native_extent!(config, target_grid)
    config.raw_file = joinpath(lake_directory(config), lake_raw_filename(config))
    fetch_lake_features(config)

    @info "Gridding the NVE depth map of lake $(config.lake_number)"
    write_lake_bathymetry(config)

    return RawBathymetry(config)
end

regrid_options(config::NVELakeBathymetryConfig) = (;
    cache = config.regrid_cache,
    interpolation_passes = config.interpolation_passes,
    major_basins = config.major_basins,
    minimum_depth = config.minimum_depth,
)

smoothing_options(config::NVELakeBathymetryConfig) = (;
    max_island_cells = config.max_island_cells,
    close_narrow_passages = config.close_narrow_passages,
    spike_ratio = config.spike_ratio,
    minimum_cell_fraction = config.minimum_cell_fraction,
    max_slope_factor = config.max_slope_factor,
    minimum_depth = config.minimum_depth,
    preserve_sills = config.preserve_sills,
)

lake_directory(config::NVELakeBathymetryConfig) = joinpath(config.data_root, config.download_directory)

lake_feature_path(config::NVELakeBathymetryConfig, kind) =
    joinpath(lake_directory(config), "lake_$(config.lake_number)_$(kind).geojson")

function lake_raw_filename(config::NVELakeBathymetryConfig)
    Nx, Ny, _ = config.size
    bounds = (config.longitude..., config.latitude...)
    label = join((replace(string(round(bound, digits = 3)), '.' => 'p') for bound in bounds), "_")
    return "nve_lake_$(config.lake_number)_$(label)_$(Nx)x$(Ny)_linear.nc"
end

"""
    fetch_lake_features(config::NVELakeBathymetryConfig)

Download the lake's depth points, depth contours and outline as GeoJSON, one file each, skipping
any already on disk. NVE's map services need no API key.
"""
function fetch_lake_features(config::NVELakeBathymetryConfig)
    mkpath(lake_directory(config))

    for (kind, layer) in pairs(NVE_LAKE_LAYERS)
        path = lake_feature_path(config, kind)
        isfile(path) && continue

        features = lake_layer_features(config, layer, "$kind of lake $(config.lake_number)")
        open(path, "w") do io
            write(io, JSON.json(Dict("type" => "FeatureCollection", "features" => features)))
        end
        @info "Saved $(length(features)) $kind features to $path"
    end

    return nothing
end

"""
    lake_layer_features(config, layer, description)

Every feature of one layer of the lake database belonging to the config's lake, in EPSG:4326,
paging until the service stops setting `exceededTransferLimit`.
"""
function lake_layer_features(config::NVELakeBathymetryConfig, layer, description)
    features = Any[]
    offset = 0

    while true
        url = string(
            NVE_LAKE_MAP_URL, "/", layer, "/query",
            "?where=vatnlnr%3D", config.lake_number,
            "&outFields=*",
            "&returnGeometry=true",
            "&outSR=4326",
            "&f=geojson",
            "&resultOffset=", offset,
            "&resultRecordCount=", NVE_LAKE_PAGE_SIZE,
        )

        page = lake_map_request(url, "$description from record $offset")
        append!(features, get(page, "features", ()))

        get(page, "exceededTransferLimit", false) === true || break
        offset += NVE_LAKE_PAGE_SIZE
    end

    return features
end

"""
    lake_map_request(url, description)

GET one ArcGIS REST query and parse it. The service reports a rejected query as HTTP 200 carrying
an `error` object, so the body is checked as well as the transport.
"""
function lake_map_request(url, description)
    @info "Downloading $description from NVE's lake database"

    filepath = try
        Downloads.download(url)
    catch exception
        exception isa Downloads.RequestError || rethrow()
        error("NVE map service request for $description failed: $(exception.message). Request was $url")
    end

    response = try
        JSON.parsefile(filepath)
    finally
        rm(filepath; force = true)
    end

    failure = get(response, "error", nothing)
    isnothing(failure) || error(
        "NVE's map service rejected the request for $description: " *
        "$(get(failure, "message", failure)). Request was $url",
    )

    return response
end

"""
    write_lake_bathymetry(config::NVELakeBathymetryConfig)

Grid the downloaded depth map onto the raw grid and write it to `config.raw_file`: the depth points
and contour vertices at their depth and the outline's vertices at 0 m, linearly interpolated, then
every cell outside the outline set to land.

A water cell inside a triangle of shore vertices alone interpolates to 0 m, which is land, and on
Mjøsa's 20 m contours that is 11 % of the lake. Its depth is somewhere between the shore and the
shallowest contour, so it takes half the shallowest contour's depth.
"""
function write_lake_bathymetry(config::NVELakeBathymetryConfig)
    Nx, Ny, _ = config.size
    samples = lake_depth_samples(config)
    length(samples) > 0 || error("NVE's lake database has no depth map for lake $(config.lake_number).")

    bottom_height = ArchGDAL.importEPSG(4326; order = :trad) do srs
        ArchGDAL.create(ArchGDAL.getdriver("Memory")) do dataset
            ArchGDAL.createlayer(
                name = "bathymetry_points",
                dataset = dataset,
                geom = ArchGDAL.wkbPoint,
                spatialref = srs,
            ) do layer
                ArchGDAL.addfielddefn!(layer, "z", ArchGDAL.OFTReal)
                add_point_features!(layer, samples)
                grid_point_dataset(dataset, config)
            end
        end
    end

    water = ArchGDAL.read(lake_feature_path(config, :outline)) do outline
        ArchGDAL.nfeature(ArchGDAL.getlayer(outline, 0)) > 0 ||
            error("NVE's lake database has no lake number $(config.lake_number).")
        rasterize_layer(outline, ArchGDAL.getname(ArchGDAL.getlayer(outline, 0)), config)
    end
    bottom_height[.!water] .= 0

    shore = water .& (bottom_height .>= 0)
    shore_depth = minimum(abs, filter(<(0), samples.bottom_heights)) / 2
    bottom_height[shore] .= -shore_depth
    @info "Lake $(config.lake_number): $(count(water)) of $(Nx * Ny) raw cells are water, " *
          "the deepest at $(-minimum(bottom_height)) m; $(count(shore)) spanned by the shore " *
          "alone are $(shore_depth) m deep"

    return write_raw_bathymetry(
        config.raw_file,
        center_coordinates(config.longitude, Nx),
        center_coordinates(config.latitude, Ny),
        bottom_height,
    )
end

"""
    lake_depth_samples(config::NVELakeBathymetryConfig)

The lake's depth points and contour vertices as bottom heights, and every vertex of its outline —
shore and islands — at 0 m, as `DepthSamples` in longitude and latitude.
"""
function lake_depth_samples(config::NVELakeBathymetryConfig)
    samples = DepthSamples()

    for kind in (:points, :contours)
        ArchGDAL.read(lake_feature_path(config, kind)) do dataset
            layer = ArchGDAL.getlayer(dataset, 0)
            ArchGDAL.nfeature(layer) > 0 || return nothing
            depth_index = ArchGDAL.findfieldindex(layer, NVE_LAKE_DEPTH_FIELD, false)
            for feature in layer
                depth = ArchGDAL.getfield(feature, depth_index)
                ismissing(depth) && continue
                collect_vertices!(samples, ArchGDAL.getgeom(feature), -abs(Float64(depth)))
            end
        end
    end
    length(samples) > 0 || return samples

    ArchGDAL.read(lake_feature_path(config, :outline)) do dataset
        for feature in ArchGDAL.getlayer(dataset, 0)
            collect_vertices!(samples, ArchGDAL.getgeom(feature), 0.0)
        end
    end

    return samples
end

"""
    collect_vertices!(samples::DepthSamples, geometry, bottom_height)

Push every vertex of `geometry` — a point, a line, a ring, or any collection of them — onto
`samples` at `bottom_height`.
"""
function collect_vertices!(samples::DepthSamples, geometry, bottom_height)
    vertex_geometry = ArchGDAL.geomname(geometry) in ("POINT", "LINESTRING", "LINEARRING")

    for index = 0:ArchGDAL.ngeom(geometry)-1
        if vertex_geometry
            x, y, _ = ArchGDAL.getpoint(geometry, index)
            push!(samples.xs, x)
            push!(samples.ys, y)
            push!(samples.bottom_heights, bottom_height)
        else
            collect_vertices!(samples, ArchGDAL.getgeom(geometry, index), bottom_height)
        end
    end

    return samples
end
