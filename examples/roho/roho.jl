# Rogaland: the coast from Jæren to the south end of Bømlo (4.46-6.61°E, 58.83-59.65°N), with
# Boknafjorden and the Ryfylke fjords, Karmsundet and the mouth of the Hardangerfjord system. A
# NorKyst-800m initial state, hourly NorKyst-800m data on the open western, southern and northern
# edges, NVE rivers and NORA3 atmosphere with snowfall, run over 2023 with T and S only.
#
# # Preparing and running it
#
#   julia --project -m FjordSim prepare_bathymetry     --config examples/roho/roho.jl
#   julia --project -m FjordSim download_initial_state --config examples/roho/roho.jl
#   julia --project -m FjordSim prepare_initial_state  --config examples/roho/roho.jl
#   julia --project -m FjordSim add_rivers             --config examples/roho/roho.jl
#   julia --project -m FjordSim download_boundary_data --config examples/roho/roho.jl
#   julia --project -m FjordSim prepare_boundary_data  --config examples/roho/roho.jl
#   julia --project -m FjordSim download_atmosphere    --config examples/roho/roho.jl
#   julia --project -m FjordSim prepare_atmosphere     --config examples/roho/roho.jl
#   julia --project -m FjordSim run_simulation         --config examples/roho/roho.jl   # needs a GPU
#   julia --project -m FjordSim validate_simulation    --config examples/roho/roho.jl
#
# `add_rivers` needs `NVE_API_KEY` (free at https://hydapi.nve.no/Users). `prepare_bathymetry`
# downloads the national Geonorge FileGDB (4.5 GB) if it is missing. Data goes to
# `~/FjordSim_data/roho/`, and results to `~/FjordSim_results/roho/`.
#
# # Settle these after `prepare_bathymetry`, before the rest
#
# If either changes, re-run every prepare step after `prepare_bathymetry`.
#
# 1. **`z_faces`.** The deepest face is fitted to the deepest cell in `bathymetry.nc`, 707.3 m. A
#    change to the bathymetry knobs can move it.
# 2. **`open_edges`.** The western edge is open sea along its whole length, the southern one west of
#    about 5.56°E, and the northern one west of Bømlo, across Bømlafjorden (5.23-5.37°E, to 357 m)
#    and at Saudafjorden's head. The eastern edge is land but for the heads of Lysefjorden and
#    Hylsfjorden, which it cuts a few kilometres short. Check them against the land mask in
#    `bathymetry.png`.
#
# # What follows `examples/oslofjorden/oslofjorden.jl`, and what does not
#
# The physics is `oslofjorden.jl`'s, every coefficient tied to the grid spacing scaled from its
# 160 m to 400 m: bounded WENO for salinity, upwinding at river mouths and in coastal columns,
# CATKE's `minimum_tke` at 7e-7, the biharmonic tracer diffusivity κ at 0, 5.06 m levels between 41
# and 132 m, its bathymetry knobs, bottom drag, open-boundary timescales and CFL. Three things do not
# carry over. Its data are the operational Norkyst-v3, published only from 2024; 2023 is in the
# operational NorKyst-800m archive. Its boundary tide is corrected to Viker; this one is NorKyst's
# own, and `scripts/measure_boundary_tide.jl` at Stavanger measures how far it is off. Its basin
# water is replaced by NIVA casts, which this box has none of.

using FjordSim
using Oceananigans:
    CATKEVerticalDiffusivity,
    HydrostaticSphericalCoriolis,
    SeawaterBuoyancy,
    WENO,
    WENOVectorInvariant
using Oceananigans.TurbulenceClosures: HorizontalScalarBiharmonicDiffusivity
using Oceananigans.Units: day, days, hour, hours, minutes, second
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using NumericalEarth: FreezingLimitedOceanTemperature
using Dates: DateTime, Month, Second
using FjordSim: Bathymetry, Plotting
using FjordSim.Configs: AbstractBathymetryConfig
using FjordSim.Grids: lateral_edges
using ArchGDAL: ArchGDAL
using NCDatasets: NCDataset
using NumericalEarth: regrid_bathymetry
using Oceananigans: interior, set!
using Oceananigans.Grids: λnodes, φnodes, Center

# --- The bathymetry -----------------------------------------------------------------------------
#
# Two corrections to the Geonorge bathymetry, both measured against the chart's depth areas
# (`dybdeareal`), whose union is the charted sea.
#
# East of about 6.51°E the FileGDB has no land polygons at most latitudes, and the gridding fills
# what it leaves open from soundings up to 19 km away, inventing water joined to the fjord heads.
# There the chart's depth areas are the shoreline: water outside them is land.
#
# Fjords whose sound is narrower than a cell come out of the 400 m regrid as basins of their own,
# and `major_basins = 1` deletes them as if they were lakes: 683 cells the chart has as sea, four
# fjords of 48-218 cells among them. `ConnectedFjords` joins them to the sea instead, through
# channels of 1-2 cells, 5-48 m deep: five basins of 4-174 cells, Skjoldafjorden the largest. Every
# water body the domain keeps is joined to the sea inside it, so the rest are land: two fjords that
# reach the sea only beyond the northern edge, one of 77 cells whose only short way to the sea runs
# along that edge's boundary row, and a charted 19-cell basin at 5.28°E 59.49°N, 1.6 km from the sea
# on the native grid too.

"""
    ConnectedFjords(source, chart_shoreline_east_of, minimum_basin_cells, maximum_channel_cells)

`source`'s bathymetry, with two changes.

East of `chart_shoreline_east_of` (°E), every wet native cell outside the chart's depth areas
(`dybdeareal`) is land.

On the model grid, the largest basin is the sea. Every other basin of `minimum_basin_cells` or more
is joined to it through the shortest 4-connected path of at most `maximum_channel_cells` land cells,
each as deep as the deepest native water inside it and at least `minimum_depth`. Smaller basins, any
with no such path, and any touching an open edge that the native grid does not join to the sea inside
the domain become land. This replaces `source`'s `major_basins`. Then `source`'s smoothing stages
run as usual, and any water they cut off from the sea becomes land.
"""
struct ConnectedFjords <: AbstractBathymetryConfig
    source::DybdedataConfig
    chart_shoreline_east_of::Float64
    minimum_basin_cells::Int
    maximum_channel_cells::Int
end

Bathymetry.bathymetry_path(config::ConnectedFjords) = Bathymetry.bathymetry_path(config.source)
Plotting.plot_path(config::ConnectedFjords) = Plotting.plot_path(config.source)

"The chart's depth areas rasterized onto `config`'s native grid: `true` where the chart has sea."
function charted_sea(config::DybdedataConfig)
    return ArchGDAL.importEPSG(25833; order = :trad) do source_srs
        ArchGDAL.importEPSG(4326; order = :trad) do target_srs
            ArchGDAL.createcoordtrans(source_srs, target_srs) do transform
                ArchGDAL.read(Bathymetry.geodatabase_path(config)) do geodatabase
                    areas = Bathymetry.find_layer(geodatabase, "dybdeareal")
                    ArchGDAL.setspatialfilter!(areas, config.filter_bounds...)
                    ArchGDAL.create(ArchGDAL.getdriver("Memory")) do memory
                        ArchGDAL.createlayer(
                            name = "sea",
                            dataset = memory,
                            geom = ArchGDAL.wkbUnknown,
                            spatialref = target_srs,
                        ) do sea
                            for area in areas
                                geometry = ArchGDAL.clone(ArchGDAL.getgeom(area))
                                ArchGDAL.transform!(geometry, transform)
                                ArchGDAL.createfeature(sea) do feature
                                    ArchGDAL.setgeom!(feature, 0, geometry)
                                    return nothing
                                end
                            end
                            return Bathymetry.rasterize_layer(memory, "sea", config)
                        end
                    end
                end
            end
        end
    end
end

function Bathymetry.bathymetry_dataset(target_grid, config::ConnectedFjords)
    dataset = Bathymetry.bathymetry_dataset(target_grid, config.source)
    raw_file = joinpath(dataset.default_download_directory, dataset.metadata_filename)
    output = replace(raw_file, r"\.nc$" => "_charted_east_of_$(config.chart_shoreline_east_of).nc")
    sea = charted_sea(config.source)

    cp(raw_file, output; force = true)
    NCDataset(output, "a") do ds
        z = ds["z"][:, :]
        uncharted = (ds["lon"][:] .> config.chart_shoreline_east_of) .& .!sea .& (z .< 0)
        z[uncharted] .= 0
        ds["z"][:, :] = z
        @info "Turned $(count(uncharted)) native cells east of $(config.chart_shoreline_east_of)°E " *
              "outside the chart's depth areas into land"
    end

    return Bathymetry.RawBathymetry(
        basename(output),
        dirname(output),
        dataset.longitude_interfaces,
        dataset.latitude_interfaces,
        dataset.size,
    )
end

"Label the 4-connected groups of `true` cells in `wet`, `0` elsewhere."
function water_bodies(wet)
    Nx, Ny = size(wet)
    labels = zeros(Int, Nx, Ny)
    count = 0
    for j in 1:Ny, i in 1:Nx
        (wet[i, j] && labels[i, j] == 0) || continue
        count += 1
        labels[i, j] = count
        queue = [(i, j)]
        while !isempty(queue)
            a, b = pop!(queue)
            for (c, d) in ((a + 1, b), (a - 1, b), (a, b + 1), (a, b - 1))
                (1 <= c <= Nx && 1 <= d <= Ny && wet[c, d] && labels[c, d] == 0) || continue
                labels[c, d] = count
                push!(queue, (c, d))
            end
        end
    end
    return labels
end

"""
    channel_to_sea(labels, basin, connected, usable, maximum_cells)

The shortest 4-connected path of land cells marked in `usable` from basin `basin` to a `usable`
cell of a basin marked in `connected`, as a vector of cells, or `nothing` if none is `maximum_cells`
long or shorter.
"""
function channel_to_sea(labels, basin, connected, usable, maximum_cells)
    Nx, Ny = size(labels)
    parent = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
    length_to = Dict{Tuple{Int,Int},Int}()
    front = [Tuple(I) for I in findall(==(basin), labels)]
    foreach(cell -> length_to[cell] = 0, front)
    while !isempty(front)
        next = Tuple{Int,Int}[]
        for (a, b) in front, (c, d) in ((a + 1, b), (a - 1, b), (a, b + 1), (a, b - 1))
            (1 <= c <= Nx && 1 <= d <= Ny && !haskey(length_to, (c, d))) || continue
            label = labels[c, d]
            if label > 0 && connected[label] && usable[c, d]
                path = Tuple{Int,Int}[]
                cell = (a, b)
                while labels[cell...] != basin
                    push!(path, cell)
                    cell = parent[cell]
                end
                return path
            end
            (label == 0 && usable[c, d] && length_to[(a, b)] < maximum_cells) || continue
            parent[(c, d)] = (a, b)
            length_to[(c, d)] = length_to[(a, b)] + 1
            push!(next, (c, d))
        end
        front = next
    end
    return nothing
end

"""
    connect_fjords(h, longitude, latitude, native, edges, config)

`h` with every basin `ConnectedFjords` keeps joined to the sea, and the others dry.

A basin that touches an open edge, and that the native grid does not join to the sea inside the
domain, reaches the sea only beyond the edge: a channel would cut through land that is there, so it
is dry too.
"""
function connect_fjords(h, longitude, latitude, native, edges, config::ConnectedFjords)
    h = copy(h)
    Nx, Ny = size(h)
    labels = water_bodies(h .< 0)
    sizes = [count(==(label), labels) for label in 1:maximum(labels)]
    sea = argmax(sizes)
    connected = falses(length(sizes))
    connected[sea] = true

    on_edge(i, j) = (:west in edges && i == 1) || (:east in edges && i == Nx) ||
                    (:south in edges && j == 1) || (:north in edges && j == Ny)
    open = falses(length(sizes))
    for j in 1:Ny, i in 1:Nx
        labels[i, j] > 0 && on_edge(i, j) && (open[labels[i, j]] = true)
    end
    # Neither through nor into an open edge's boundary row: `raise_boundary_pockets` turns a cell
    # there with land inside it into land, cutting the channel.
    usable = [!on_edge(i, j) for i in 1:Nx, j in 1:Ny]

    # The native water body under each model cell centre.
    native_labels = water_bodies(native.h .< 0)
    columns = [argmin(abs.(native.longitude .- λ)) for λ in longitude]
    rows = [argmin(abs.(native.latitude .- φ)) for φ in latitude]
    native_body(i, j) = native_labels[columns[i], rows[j]]
    native_sea = Set(native_body(i, j) for j in 1:Ny, i in 1:Nx if labels[i, j] == sea)
    delete!(native_sea, 0)
    natively_joined(basin) =
        any(native_body(i, j) in native_sea for j in 1:Ny, i in 1:Nx if labels[i, j] == basin)

    Δλ, Δφ = longitude[2] - longitude[1], latitude[2] - latitude[1]
    minimum_depth = config.source.minimum_depth
    function channel_depth(i, j)
        footprint_columns = findall(λ -> abs(λ - longitude[i]) <= Δλ / 2, native.longitude)
        footprint_rows = findall(φ -> abs(φ - latitude[j]) <= Δφ / 2, native.latitude)
        deepest = minimum(native.h[footprint_columns, footprint_rows]; init = 0)
        return min(deepest, -minimum_depth)
    end

    for basin in sortperm(sizes; rev = true)
        connected[basin] && continue
        name = "the $(sizes[basin])-cell basin at " *
               "$(round(longitude[first(Tuple(findfirst(==(basin), labels)))]; digits = 3))°E"
        if sizes[basin] < config.minimum_basin_cells
            h[labels .== basin] .= 0
            continue
        end
        if open[basin] && !natively_joined(basin)
            h[labels .== basin] .= 0
            @info "Removed $name, which reaches the sea only beyond the edge"
            continue
        end
        path = channel_to_sea(labels, basin, connected, usable, config.maximum_channel_cells)
        if isnothing(path)
            h[labels .== basin] .= 0
            @warn "No channel of $(config.maximum_channel_cells) cells or fewer joins $name to the sea; it is land"
            continue
        end
        for (i, j) in path
            h[i, j] = channel_depth(i, j)
            labels[i, j] = basin
        end
        connected[basin] = true
        @info "Joined $name to the sea through $(length(path)) cells, " *
              "$(round(-maximum(h[i, j] for (i, j) in path); digits = 1)) m deep at its shallowest"
    end
    return h
end

"`h` with every water body but the largest turned into land."
function sea_only(h)
    labels = water_bodies(h .< 0)
    sea = argmax([count(==(label), labels) for label in 1:maximum(labels)])
    cut_off = (labels .> 0) .& (labels .!= sea)
    any(cut_off) && @warn "Turned $(count(cut_off)) cells the smoothing cut off from the sea into land"
    return ifelse.(cut_off, zero(eltype(h)), h)
end

function Bathymetry.process_bathymetry(target_grid, config::ConnectedFjords; edges = nothing)
    dataset = Bathymetry.bathymetry_dataset(target_grid, config)
    metadata = Bathymetry.Metadatum(:bottom_height; dataset)
    options = (; Bathymetry.regrid_options(config.source)..., major_basins = Inf)
    bottom_height = regrid_bathymetry(target_grid, metadata; options...)

    native = Bathymetry.native_bathymetry(metadata)
    longitude = collect(λnodes(target_grid, Center()))
    latitude = collect(φnodes(target_grid, Center()))
    h = Array(interior(bottom_height, :, :, 1))
    set!(bottom_height, connect_fjords(h, longitude, latitude, native, lateral_edges(edges), config))

    smoothing = Bathymetry.smoothing_options(config.source)
    Bathymetry.smooth_bathymetry_gaps!(bottom_height; smoothing..., native, edges)
    set!(bottom_height, sea_only(Array(interior(bottom_height, :, :, 1))))
    output_file = Bathymetry.bathymetry_path(config)
    Bathymetry.write_bathymetry_file(output_file, target_grid, bottom_height)

    return (; dataset, raw_file = Bathymetry.metadata_path(metadata), output_file, bottom_height)
end

data_root = fjord_data_root("roho")
FT = Float32
# The calendar year 2023. Every dataset below is downloaded and prepared for the same period.
start_date = DateTime(2023, 1, 1)
stop_date = DateTime(2024, 1, 1)
# The full-field writers start a new file on the 1st of each month.
monthly_splits = [(start_date + Month(month) - start_date) / Second(1) for month in 1:11]
# `oslofjorden.jl`'s levels: a 1 m surface cell stretching by 1.25, 41-132 m split into eighteen
# 5.06 m levels against WENO's implicit vertical diffusion in the sill basins, and 33.5 m layers
# below, here continued to -668 m, with a 39.5 m bottom layer reaching the deepest cell at 707.3 m
# (gate 1 in the header). If the thin levels hold the time step
# down over the steep fjord walls, as on `examples/moro/moro.jl`, 10 m levels between 41 and 131 m
# are that setup's answer.
#
# 46 levels on 305 x 228 columns is 3.2 M cells.
z_faces = vcat(
    -707.5,
    range(-668.0, -165.5; step = 33.5),
    range(-132.0, -41.0; length = 19),
    [-32.0, -25.0, -19.0, -14.5, -10.8, -7.9, -5.5, -3.7, -2.2, -1.0, 0.0],
)

FjordConfig(
    grid_config = EvenGrid(
        # 0.00705° by 0.0036°: 406 m by 400 m at 58.83°N and 397 m by 400 m at 59.65°N, 2.5 times
        # Oslofjord's 160 m. Every coefficient below tied to Δx is scaled by that factor.
        size       = (305, 228, length(z_faces) - 1),
        halo       = (7, 7, 7),
        longitude  = (4.46, 6.61),
        latitude   = (58.83, 59.65),
        # `z_faces` is written into `bathymetry.nc` and `simulation_grid` reads the file, so a
        # change here needs `prepare_bathymetry` re-run.
        z_faces    = z_faces,
        float_type = FT,
    ),
    # `oslofjorden.jl`'s knobs but for `minimum_depth` and `open_boundary_land_cells`, with the
    # chart's shoreline east of 6.5°E and every fjord joined to the sea ("The bathymetry" above).
    bathymetry_config = ConnectedFjords(
        DybdedataConfig(
            data_root             = data_root,
            output_file           = "bathymetry.nc",
            plot_file             = "bathymetry.png",
            # A 100 m native grid.
            raw_resolution_factor = 4,
            padding_cells         = 2,
            include_contours      = false,
            contour_stride        = 10,
            interpolation_passes  = 1,
            # 5 m rather than Oslofjord's 2, as on `examples/lista/lista.jl`. A river enters
            # through the top cell of its columns, and the shallower the column the less salt
            # water there is to take it: in the Oslofjord runs, 2-3 m columns beside a river mouth
            # are where WENO drove salinity below TEOS10's floor, and where the water overheated.
            minimum_depth         = 5.0,
            # Off, unlike Oslofjord's 5: the Jæren coast runs into the southern edge and Bømlo into
            # the northern one, and the stage floods outwards along the whole edge band.
            open_boundary_land_cells = 0,
            max_island_cells      = 6,
            close_narrow_passages = true,
            spike_ratio           = 0.5,
            # Must equal the `minimum_fractional_cell_height` the grid gives `PartialCellBottom`.
            minimum_cell_fraction = 0.2,
            max_slope_factor      = 0.5,
            preserve_sills        = true,
            raise_boundary_pockets = true,
            geonorge_cache        = true,
            regrid_cache          = false,
        ),
        6.5,
        # Basins of 4 cells (0.64 km²) or more, through channels of up to 3 cells.
        4,
        3,
    ),
    # One NorKyst daily mean, the day `start_date` falls in, for `FromInitialState()` below.
    initial_state_config = NorKystInitialStateConfig(
        data_root        = data_root,
        output_directory = "norkyst_initial",
        output_file      = "initial_state.nc",
        parameters       = ["temperature", "salinity", "u_eastward", "v_northward"],
        start_date       = start_date,
    ),
    river_config = NVERiversConfig(
        data_root   = data_root,
        output_file = "rivers_nve.nc",
        plot_file   = "rivers_nve.png",
        start_date  = start_date,
        stop_date   = stop_date,
        # Mouths are discovered from NVE's ELVIS network and REGINE catchments. At 0.5 m³/s that
        # finds 64, Suldalsvassdraget (114 m³/s), Årdalselva (44), Saudavassdraget (35) and Ulla
        # (35) the largest. Figgjo reaches the sea just south of the box.
        minimum_discharge = 0.5,
        # The vertical Courant bound on the 1 m top cell: 1 mm/s is w·Δt/Δz = 0.18 at the 3-minute
        # `max_time_step`.
        maximum_inflow_velocity = 1.0e-3,
        # Overrides, keyed by each mouth's terminal ELVIS `vassdragsnr`. Every series named here has
        # daily values on at least 362 days of 2023.
        #
        # Suldal, Ulla and Årdal keep their catchment normals. Ulla-Førre transfers 300-2000 km² of
        # each to Blåsjø, and that water reaches the sea through Hylen kraftstasjon at the head of
        # Hylsfjorden, which publishes no discharge. Their mouth gauges carry only what is left:
        # in 2023 Lavika (36.6.0) 48 m³/s, Hauge bru (35.2.0) 4.8 and Leirberget (33.8.0) 18.
        # The normals keep the freshwater volume, at natural timing.
        outlets = [
            # Holmen (38.1.0) drains 117 of the unregulated catchment's 118 km², and is used
            # unscaled. Temperature at Vikedalselva utløp (38.2.0).
            NVERiver(
                vassdragsnr = "038.21D", name = "Vikedalselva",
                discharge_station = "38.1.0", temperature_station = "38.2.0",
            ),
            NVERiver(vassdragsnr = "036.A21", name = "Suldalslågen", temperature_station = "36.61.0"),
            NVERiver(vassdragsnr = "033.A0", name = "Årdalselva", temperature_station = "33.8.0"),
            NVERiver(vassdragsnr = "035.A2", name = "Ulla", temperature_station = "35.2.0"),
        ],
    ),
    boundary_data_config = NorKystBoundaryDataConfig(
        data_root        = data_root,
        output_directory = "norkyst_hourly",
        output_file      = "boundaries.nc",
        plot_file        = "boundaries.png",
        # The Norwegian Coastal Current enters through the southern and western edges and leaves
        # through the northern one, which also carries the Hardangerfjord's exchange (gate 2 in
        # the header). The bounding box of three edge bands is the whole domain, so the hourly
        # download covers all of it.
        open_edges       = (:west, :south, :north),
        margin           = 0.05,
        architecture     = :auto,
        parameters       = [
            "temperature", "salinity", "u_eastward", "v_northward", "zeta", "ubar", "vbar",
        ],
        start_date       = start_date,
        stop_date        = stop_date,
    ),
    atmosphere_config = NORA3Config(
        data_root        = data_root,
        output_directory = "nora3",
        output_file      = "atmosphere.nc",
        plot_file        = "atmosphere.png",
        resolution       = 0.02,
        padding          = 0.1,
        start_date       = start_date,
        stop_date        = stop_date,
        # Snow falls into the water, which supplies the heat that melts it.
        snowfall         = true,
    ),
    model_config = CoupledHydrostaticModel(
        buoyancy           = SeawaterBuoyancy(FT, equation_of_state = TEOS10EquationOfState(FT)),
        closure            = BoundarySponge(
            base = (
                # 7e-7, as in `oslofjorden.jl`, where 7e-6 left the surface layer too salty and too
                # weakly stratified. Unmeasured here: check `e` against the floor.
                CATKEVerticalDiffusivity(FT; minimum_tke = 7e-7),
                # Oslofjord's 9.4e3 m⁴ s⁻¹ scaled by (400/160)⁴, which keeps its 72 min damping of
                # the 2Δx mode. The explicit limit Δt ≤ Δx⁴/32ν₄ is 2100 s at the narrowest 397 m.
                #
                # κ is 0, as there: κ = 2e3 overshot at fronts and left wet cells below 0 g/kg.
                HorizontalScalarBiharmonicDiffusivity(FT; ν = 3.7e5, κ = 0),
            ),
            # Oslofjord's 3.2 km.
            width_cells = 8,
            # Its 20 and 10 m² s⁻¹ scaled by (400/160)². Δt ≤ Δx²/4ν is 315 s, above
            # `max_time_step`.
            viscosity   = 125.0,
            diffusivity = 62.5,
        ),
        # As in `oslofjorden.jl`: T and S upwinded horizontally within WENO5's three-cell stencil
        # of every river column and in coastal columns, where WENO across a fresh front left
        # cells colder or saltier than any neighbour; salinity's WENO bounded, against
        # TEOS10's -32 psu floor. A NamedTuple leaves any tracer it omits on `Centered()`, so it
        # names CATKE's `e` as well.
        tracer_advection   = (
            T = RiverMouthAdvection(interior = WENO(FT), width_cells = 3, narrow_cells = true),
            S = RiverMouthAdvection(interior = WENO(FT; bounds = (0, 40)), width_cells = 3, narrow_cells = true),
            e = WENO(FT),
        ),
        momentum_advection = WENOVectorInvariant(FT),
        tracers            = (:T, :S),
        coriolis           = HydrostaticSphericalCoriolis(FT),
        sea_ice            = FreezingLimitedOceanTemperature(FT),
        free_surface       = SplitExplicitFreeSurfaceConfig(cfl = 0.7),
        extra_kwargs       = (;),
    ),
    boundary_condition_configs = (
        AirSeaFluxes(),
        QuadraticBottomDrag(coefficient = 0.003),
        # Marchesiello et al. (2001): nudge hard on inflow, radiate on outflow. Oslofjord's
        # timescales, applied to every open edge.
        OpenLateralBoundaryFromData(
            inflow_timescale  = 3hours,
            outflow_timescale = 360days,
        ),
    ),
    run_config = RunConfig(
        results_root       = fjord_results_root("roho"),
        architecture       = :auto,
        writers = (
            # Daily, a new file every month. One record of these five fields is about 64 MB, so
            # the year is about 23 GB.
            SnapshotWriter(
                name               = :ocean,
                output_file        = "snapshots_ocean.nc",
                variables          = (:T, :S, :u, :v, :e),
                interval           = 1day,
                split_times        = monthly_splits,
                overwrite_existing = true,
            ),
            # The top level of T, S, u and v, and `η`, hourly for the tide, a new file every month:
            # about 1 GB a month.
            SurfaceSnapshotWriter(
                name               = :surface,
                output_file        = "snapshots_surface.jld2",
                variables          = (:η, :T, :S, :u, :v),
                interval           = 1hour,
                split_times        = monthly_splits,
                overwrite_existing = true,
            ),
            # The two permanent Kartverket gauges in the box, from Kartverket's station register
            # (`vannstand.kartverket.no/tideapi.php?tide_request=stationlist&type=perm`). Ten
            # minutes for the harmonic analysis `validate_simulation` does against their record.
            FieldStationWriter(
                name               = :tidegauge,
                output_file        = "stations_tidegauge.jld2",
                variables          = (:η,),
                stations           = [
                    Station(name = "Stavanger", longitude = 5.730121, latitude = 58.974339),
                    Station(name = "Sandnes",   longitude = 5.746613, latitude = 58.868232),
                ],
                interval           = 10minutes,
                overwrite_existing = true,
            ),
            CheckpointWriter(interval = 12hours, cleanup = true),
        ),
        callbacks = (ProgressCallback(name = :progress, interval = 1hour, report = progress),),
        time_stepping = AdaptiveTimeStep(
            initial_time_step    = 1second,
            cfl                  = 0.3,
            max_time_step        = 3minutes,
            max_time_step_change = 1.01,
        ),
        # The NorKyst daily mean of `start_date`: T, S, u and v.
        initial_conditions = FromInitialState(),
        start_date         = start_date,
        stop_date          = stop_date,
        loops              = 1,
        pickup             = false,
    ),
)
