# A video of a run's ocean snapshots (`SnapshotWriter`) along a vertical section that follows the
# deepest water between two points, one frame per record over all the given files in time order —
# e.g. every monthly part of a daily writer for a year:
#
#     julia --project scripts/section_videos.jl \
#         --from 10.45,61.10 --to 11.22,60.41 [--via LON,LAT]... [--variables T,velocity] [--upper 50] \
#         [--out DIR] ~/FjordSim_results/mjosa/snapshots_ocean_<tag>_part*.nc \
#         [~/FjordSim_results/mjosa/snapshots_surface_<tag>_part*.jld2] \
#         [~/FjordSim_results/mjosa/snapshots_sea_ice_<tag>_part*.jld2]
#
# The path joins the wet cells nearest to `--from`, each `--via` and `--to` by the route of least
# Σ step length × (deepest cell / cell depth)², which keeps to the deep channel rather than the
# shortest line. Each variable gets two rows: the top `--upper` metres (default 50) and the full
# depth, with the bottom along the path in black. A map of the path sits beside them. The colour
# range is the 1st–99th percentile of the wet values over the files (symmetric about zero for the
# velocity), fixed for the whole video so the seasons compare.
#
# The variables are those with a value per cell and record, except `u` and `v`, which give
# `velocity`: the component along the path, + in the direction of increasing distance, from `u` and
# `v` averaged onto the cell centres and the path's direction over ±1 km. Without `--variables`, every one of them
# whose value is not the same everywhere (e.g. a lake's zero salinity).
#
# The `.jld2` files are optional and recognised by their contents. A `SurfaceSnapshotWriter`'s
# (holding `η`) adds a strip of the surface elevation along the path above the velocity, as the
# departure from the mean over the water; a `SeaIceSnapshotWriter`'s (`h`, `ℵ`, and `hs` if
# written) adds one of the ice and snow thickness × concentration above the temperature. They must
# hold a record at every ocean snapshot's time.
#
# The ocean files carry everything else: the model's bottom and land mask, the levels, and the start
# date the seconds of model time count from. The video goes beside the first file unless `--out` is
# given.
module SectionVideos

using Dates: DateTime, Millisecond, format
using JLD2: jldopen
using NCDatasets: NCDataset
using Statistics: mean, quantile
import CairoMakie as M

const FRAMES_PER_SECOND = 10
const COLOUR_RANGE_RECORDS = 50
const EARTH_RADIUS_KM = 6371.0
const CONTOUR_LEVELS = 10
# The path's direction at a cell is from the cell this far before it to the one this far after it.
const DIRECTION_HALF_SPAN_KM = 1.0

const LABELS = Dict(
    "T" => (label = "Temperature (°C)", colormap = :thermal, symmetric = false),
    "S" => (label = "Salinity (g kg⁻¹)", colormap = :haline, symmetric = false),
    "velocity" => (label = "Velocity along the path (m s⁻¹), + in the direction of increasing distance",
                   colormap = :balance, symmetric = true),
)
# The strip drawn above a variable's sections when its files are given.
const STRIP_ABOVE = Dict("T" => :sea_ice, "velocity" => :elevation)

# A cell is wet when its top level is active; its depth is the model's bottom below it.
function snapshot_grid(file)
    return NCDataset(file) do ds
        wet = ds["inactive_nodes_ccc"][:, :, end] .== 0
        depth = ifelse.(wet, max.(-Float64.(ds["bottom_height"][:, :]), 0.0), 0.0)
        (; longitude = ds["λ_caa"][:], latitude = ds["φ_aca"][:], depth)
    end
end

kilometres(a, b) = 2 * EARTH_RADIUS_KM * asin(sqrt(sind((b[2] - a[2]) / 2)^2 +
                                                   cosd(a[2]) * cosd(b[2]) * sind((b[1] - a[1]) / 2)^2))

centre(grid, I) = (grid.longitude[I[1]], grid.latitude[I[2]])

function nearest_wet_cell(grid, point)
    wet = findall(>(0), grid.depth)
    return wet[argmin([kilometres(centre(grid, I), point) for I in wet])]
end

# A binary min-heap of (cost, cell) for `deepest_route`.
function heap_push!(heap, item)
    push!(heap, item)
    i = length(heap)
    while i > 1 && heap[i÷2][1] > heap[i][1]
        heap[i÷2], heap[i] = heap[i], heap[i÷2]
        i ÷= 2
    end
    return heap
end

function heap_pop!(heap)
    top = heap[1]
    heap[1] = heap[end]
    pop!(heap)
    i = 1
    while true
        smallest = i
        for child in (2i, 2i + 1)
            child <= length(heap) && heap[child][1] < heap[smallest][1] && (smallest = child)
        end
        smallest == i && return top
        heap[i], heap[smallest] = heap[smallest], heap[i]
        i = smallest
    end
end

"""
    deepest_route(grid, start, stop)

The wet cells from `start` to `stop` (8-connected) of least Σ step length × (deepest cell / cell depth)².
"""
function deepest_route(grid, start, stop)
    deepest = maximum(grid.depth)
    cost = fill(Inf, size(grid.depth))
    previous = fill(CartesianIndex(0, 0), size(grid.depth))
    cost[start] = 0
    heap = [(0.0, start)]
    steps = [CartesianIndex(di, dj) for di in -1:1, dj in -1:1 if (di, dj) != (0, 0)]
    while !isempty(heap)
        here_cost, here = heap_pop!(heap)
        here == stop && break
        here_cost > cost[here] && continue
        for step in steps
            next = here + step
            checkbounds(Bool, grid.depth, next) && grid.depth[next] > 0 || continue
            next_cost = here_cost + kilometres(centre(grid, here), centre(grid, next)) *
                                    (deepest / grid.depth[next])^2
            if next_cost < cost[next]
                cost[next] = next_cost
                previous[next] = here
                heap_push!(heap, (next_cost, next))
            end
        end
    end
    isinf(cost[stop]) && error("No wet route joins $(centre(grid, start)) and $(centre(grid, stop))")
    route = [stop]
    while route[end] != start
        push!(route, previous[route[end]])
    end
    return reverse(route)
end

# The unit vector (east, north) of the path at each cell, with at least one cell on either side.
function path_directions(grid, cells, distance)
    return map(eachindex(cells)) do n
        before = clamp(searchsortedfirst(distance, distance[n] - DIRECTION_HALF_SPAN_KM), 1, max(n - 1, 1))
        after = clamp(searchsortedlast(distance, distance[n] + DIRECTION_HALF_SPAN_KM),
                      min(n + 1, length(cells)), length(cells))
        (λ₁, φ₁), (λ₂, φ₂) = centre(grid, cells[before]), centre(grid, cells[after])
        east, north = (λ₂ - λ₁) * cosd((φ₁ + φ₂) / 2), φ₂ - φ₁
        (east, north) ./ hypot(east, north)
    end
end

function section_path(grid, points)
    cells = [nearest_wet_cell(grid, point) for point in points]
    path = [cells[1]]
    for (start, stop) in zip(cells[1:end-1], cells[2:end])
        append!(path, deepest_route(grid, start, stop)[2:end])
    end
    distance = cumsum([0; [kilometres(centre(grid, a), centre(grid, b)) for (a, b) in zip(path[1:end-1], path[2:end])]])
    return (; cells = path, distance, depth = [grid.depth[I] for I in path], grid_size = size(grid.depth),
            direction = path_directions(grid, path, distance))
end

"""
    read_section(paths, variable, path, grid)

`variable` along the path at every record of the files: cell × level × record, surface level first
and NaN in the model's inactive cells, with the records' model times and the levels' depths.
"""
function read_section(paths, variable, path)
    columns = Matrix{Float32}[]
    times = Float64[]
    faces, inactive = NCDataset(first(paths)) do ds
        reverse(-ds["z_aaf"][:]), ds["inactive_nodes_ccc"][:, :, :]
    end
    centres = (faces[1:end-1] .+ faces[2:end]) ./ 2
    below_bottom = [inactive[I[1], I[2], end+1-k] != 0 for I in path.cells, k in eachindex(centres)]
    reader = get(DERIVED, variable) do
        (ds, record, path, levels) -> stored_section(ds, variable, record, path, levels)
    end
    for file in paths
        NCDataset(file) do ds
            for record in eachindex(ds["time"])
                section = reader(ds, record, path, eachindex(centres))
                section[below_bottom] .= NaN
                push!(columns, section)
            end
            append!(times, ds["time"][:])
        end
    end
    return (; values = cat(columns...; dims = 3), times, faces, centres)
end

# A field on the cell faces of one horizontal direction has one more point along it than the cells;
# the mean of each cell's two faces is on its centre.
function centred(field, grid_size)
    size(field, 1) == grid_size[1] + 1 && (field = (field[1:end-1, :, :] .+ field[2:end, :, :]) ./ 2)
    size(field, 2) == grid_size[2] + 1 && (field = (field[:, 1:end-1, :] .+ field[:, 2:end, :]) ./ 2)
    return field
end

stored_section(ds, variable, record, path, levels) =
    (field = centred(ds[variable][:, :, :, record], path.grid_size);
     Float32[field[I[1], I[2], end+1-k] for I in path.cells, k in levels])

function along_path_velocity(ds, record, path, levels)
    u = stored_section(ds, "u", record, path, levels)
    v = stored_section(ds, "v", record, path, levels)
    return u .* first.(path.direction) .+ v .* last.(path.direction)
end

# Variables computed from the stored ones, and the stored ones they replace.
const DERIVED = Dict("velocity" => along_path_velocity)
const REPLACED = Dict("velocity" => ["u", "v"])

uniform(section) = (values = filter(isfinite, section.values); isempty(values) || all(==(first(values)), values))

function strip_kind(file)
    names = jldopen(f -> keys(f["timeseries"]), file)
    "η" in names && return :elevation
    "h" in names && "ℵ" in names && return :sea_ice
    return error("$file holds neither η nor sea ice (h and ℵ)")
end

at_path(f, name, iteration, path) = (field = f["timeseries/$name/$iteration"];
                                     Float32[field[I[1], I[2], 1] for I in path.cells])

# The departure from the mean over the water and that mean, both in cm.
function strip_record(::Val{:elevation}, f, iteration, path, wet)
    η = f["timeseries/η/$iteration"][:, :, 1]
    water_mean = mean(η[wet])
    return (anomaly = 100 .* (Float32[η[I] for I in path.cells] .- water_mean), mean = 100 * water_mean)
end

# `h` and `hs` are the thicknesses over the ice-covered fraction `ℵ`.
function strip_record(::Val{:sea_ice}, f, iteration, path, wet)
    concentration = at_path(f, "ℵ", iteration, path)
    ice = at_path(f, "h", iteration, path) .* concentration
    snow = "hs" in keys(f["timeseries"]) ? at_path(f, "hs", iteration, path) .* concentration : zero(ice)
    return (; ice, snow)
end

"""
    read_strip(kind, files, times, path, wet)

The strip of `kind` along the path at each of the model `times`, from the JLD2 snapshot `files`.
"""
function read_strip(kind, files, times, path, wet)
    wanted = Dict(round(Int, t) => n for (n, t) in enumerate(times))
    found = Vector{Any}(nothing, length(times))
    for file in files
        jldopen(file) do f
            for iteration in keys(f["timeseries/t"])
                n = get(wanted, round(Int, f["timeseries/t/$iteration"]), 0)
                n > 0 && (found[n] = strip_record(Val(kind), f, iteration, path, wet))
            end
        end
    end
    absent = count(isnothing, found)
    absent > 0 && error("The $kind files lack $absent of the $(length(times)) ocean snapshot times")
    return [found...]
end

"""
    strip!(figure, row, kind, records, path, xlimits)

Draw the strip of `kind` in `figure`'s `row` and return the function that shows a record.
"""
function strip!(figure, row, ::Val{:sea_ice}, records, path, xlimits)
    bottom = zeros(Float32, length(path.cells))
    ice = M.Observable(records[1].ice)
    top = M.Observable(records[1].ice .+ records[1].snow)
    highest = maximum(record -> maximum(filter(isfinite, record.ice .+ record.snow); init = 0), records)
    axis = M.Axis(figure[row, 1]; title = "Ice and snow × concentration", ylabel = "Thickness (m)",
                  yticks = M.WilkinsonTicks(3))
    M.band!(axis, path.distance, bottom, ice; color = :deepskyblue3, label = "ice")
    M.band!(axis, path.distance, ice, top; color = :gray70, label = "snow")
    M.limits!(axis, xlimits..., 0, max(1.1 * highest, 0.05))
    M.Legend(figure[row, 2], axis; framevisible = false)
    return function (record)
        ice[] = records[record].ice
        top[] = records[record].ice .+ records[record].snow
    end
end

elevation_title(record) = "Surface elevation − mean over the water, mean η = $(round(record.mean; digits = 1)) cm"

function strip!(figure, row, ::Val{:elevation}, records, path, xlimits)
    anomaly = M.Observable(records[1].anomaly)
    title = M.Observable(elevation_title(records[1]))
    low, high = quantile(filter(isfinite, reduce(vcat, [record.anomaly for record in records])), [0.01, 0.99])
    axis = M.Axis(figure[row, 1]; title, ylabel = "η − mean (cm)", yticks = M.WilkinsonTicks(3))
    M.hlines!(axis, 0; color = :gray50, linewidth = 0.7)
    M.lines!(axis, path.distance, anomaly; color = :black, linewidth = 1.5)
    M.limits!(axis, xlimits..., (-1, 1) .* max(abs(low), abs(high), 0.1)...)
    return function (record)
        anomaly[] = records[record].anomaly
        title[] = elevation_title(records[record])
    end
end

function colour_range(values)
    records = unique(round.(Int, range(1, size(values, 3); length = COLOUR_RANGE_RECORDS)))
    samples = filter(isfinite, vec(values[:, :, records]))
    low, high = quantile(samples, [0.01, 0.99])
    return low < high ? (low, high) : (low - 1, high + 1)
end

function map_panel!(position, grid, path)
    (; longitude, latitude, depth) = grid
    aspect = (last(longitude) - first(longitude)) * cosd((first(latitude) + last(latitude)) / 2) /
             (last(latitude) - first(latitude))
    axis = M.Axis(position; title = "Section path", xlabel = "Longitude (°E)", ylabel = "Latitude (°N)",
                  aspect = M.AxisAspect(aspect))
    M.heatmap!(axis, longitude, latitude, ifelse.(depth .> 0, depth, NaN); colormap = :deep,
               nan_color = :gray70)
    points = [centre(grid, I) for I in path.cells]
    M.lines!(axis, first.(points), last.(points); color = :red, linewidth = 2)
    M.scatter!(axis, [first(points)]; color = :red, markersize = 12)
    M.text!(axis, first(points)...; text = "  0 km", color = :red, align = (:left, :center))
    M.limits!(axis, extrema(longitude)..., extrema(latitude)...)
    return axis
end

function record_video(file, sections, strips, path, grid, title_prefix, first_date, upper)
    edges = [path.distance[1]; (path.distance[1:end-1] .+ path.distance[2:end]) ./ 2; path.distance[end]]
    times = last(first(sections)).times
    date(record) = first_date + Millisecond(round(Int, 1000 * times[record]))
    title = M.Observable("")
    deepest = maximum(path.depth)

    strip_of(variable) = get(strips, get(STRIP_ABOVE, variable, nothing), nothing)
    height = sum(isnothing(strip_of(variable)) ? 520 : 680 for (variable, _) in sections) + 80
    figure = M.Figure(size = (1900, height), fontsize = 18)
    M.Label(figure[0, 1:3], title; fontsize = 22, font = :bold)
    frames = Dict{String,Any}()
    show_strips = Function[]
    # One figure grid for every row, so the axes and colour bars line up. Row heights: a strip 0.6,
    # the top metres 1.2 and the full depth 1.8.
    row = 0
    for (n, (variable, section)) in enumerate(sections)
        (; label, colormap, symmetric) = get(LABELS, variable,
                                             (label = variable, colormap = :viridis, symmetric = false))
        colorrange = colour_range(section.values)
        symmetric && (colorrange = (-1, 1) .* maximum(abs, colorrange))
        levels = M.Makie.get_tickvalues(M.LinearTicks(CONTOUR_LEVELS), colorrange...)
        frame = M.Observable(section.values[:, :, 1])
        frames[variable] = frame
        strip = strip_of(variable)
        if !isnothing(strip)
            row += 1
            xlimits = (first(edges), last(edges))
            push!(show_strips, strip!(figure, row, Val(STRIP_ABOVE[variable]), strip, path, xlimits))
            M.rowsize!(figure.layout, row, M.Auto(0.6))
        end
        local plot
        for (panel, (bottom, name, weight)) in enumerate(((upper, "top $(round(Int, upper)) m", 1.2),
                                                          (deepest, "full depth", 1.8)))
            axis = M.Axis(figure[row+panel, 1]; title = "$label, $name", yreversed = true, ylabel = "Depth (m)",
                          xlabel = panel == 2 && n == length(sections) ? "Distance along the path (km)" : "")
            plot = M.heatmap!(axis, edges, section.faces, frame; colormap, colorrange, nan_color = :gray80)
            M.contour!(axis, path.distance, section.centres, frame; levels, color = :black, linewidth = 0.7)
            M.lines!(axis, path.distance, path.depth; color = :black, linewidth = 1.5)
            M.limits!(axis, first(edges), last(edges), bottom, 0)
            M.rowsize!(figure.layout, row + panel, M.Auto(weight))
        end
        # The contours are at the colour bar's ticks.
        M.Colorbar(figure[row+1:row+2, 2], plot; label, ticks = levels)
        row += 2
    end
    map_panel!(figure[1:row, 3], grid, path)
    M.colsize!(figure.layout, 3, M.Relative(0.22))

    M.record(figure, file, eachindex(times); framerate = FRAMES_PER_SECOND) do record
        for (variable, section) in sections
            frames[variable][] = section.values[:, :, record]
        end
        foreach(show -> show(record), show_strips)
        title[] = "$title_prefix, " * format(date(record), "yyyy-mm-dd HH:MM")
    end
    println("Wrote $file")
    return file
end

part_number(file) = (found = match(r"_part(\d+)\.nc$", file); isnothing(found) ? 0 : parse(Int, found[1]))

parse_point(text) = Tuple(parse.(Float64, split(text, ",")))

function main(args)
    flags = ("--from", "--to", "--via", "--variables", "--upper", "--out")
    options = Dict{String,Vector{String}}()
    files = String[]
    i = 1
    while i <= length(args)
        if args[i] in flags
            push!(get!(options, args[i], String[]), args[i+1])
            i += 2
        else
            push!(files, args[i])
            i += 1
        end
    end
    option(flag) = only(get(options, flag, [nothing]))
    for flag in ("--from", "--to")
        isnothing(option(flag)) && error("$flag is required")
    end
    strip_files = filter(endswith(".jld2"), files)
    files = sort(filter(endswith(".nc"), files); by = part_number)
    isempty(files) && error("Give the ocean snapshot files")

    grid = snapshot_grid(first(files))
    first_date = NCDataset(ds -> DateTime(ds.attrib["start_date"]), first(files))
    points = [parse_point(option("--from")); parse_point.(get(options, "--via", String[])); parse_point(option("--to"))]
    path = section_path(grid, points)
    println("Path: $(length(path.cells)) cells, $(round(last(path.distance); digits = 1)) km, " *
            "deepest $(round(maximum(path.depth); digits = 1)) m")

    stored = NCDataset(ds -> [name for name in keys(ds) if ndims(ds[name]) == 4], first(files))
    derived = [name for (name, sources) in REPLACED if sources ⊆ stored]
    replaced = reduce(vcat, [REPLACED[name] for name in derived]; init = String[])
    available = [setdiff(stored, replaced); derived]
    requested = option("--variables")
    variables = split(something(requested, join(available, ",")), ",")
    sections = [variable => read_section(files, variable, path) for variable in variables]
    if isnothing(requested)
        for (variable, section) in filter(uniform ∘ last, sections)
            println("Skipping $variable: the same value everywhere")
        end
        filter!(!uniform ∘ last, sections)
    end
    variables = join(first.(sections), ",")

    times = last(first(sections)).times
    wet = grid.depth .> 0
    kinds = strip_kind.(strip_files)
    strips = Dict(kind => read_strip(kind, strip_files[kinds .== kind], times, path, wet) for kind in unique(kinds))

    upper = parse(Float64, something(option("--upper"), "50"))
    out = something(option("--out"), dirname(abspath(first(files))))
    mkpath(out)
    stem = replace(first(splitext(basename(first(files)))), r"_part\d+$" => "")
    title_prefix = "$(basename(dirname(abspath(first(files))))) section"
    record_video(joinpath(out, "$(stem)_section_$(replace(variables, "," => "_")).mp4"), sections, strips, path,
                 grid, title_prefix, first_date, upper)
    return nothing
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    SectionVideos.main(ARGS)
end
