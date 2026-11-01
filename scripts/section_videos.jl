# A video of a run's ocean snapshots (`SnapshotWriter`) along a vertical section that follows the
# deepest water between two points, one frame per record over all the given files in time order —
# e.g. every monthly part of a daily writer for a year:
#
#     julia --project scripts/section_videos.jl \
#         --from 10.45,61.10 --to 11.22,60.41 [--via LON,LAT]... [--variables T,S] [--upper 50] \
#         [--out DIR] ~/FjordSim_results/mjosa/snapshots_ocean_<tag>_part*.nc
#
# The path joins the wet cells nearest to `--from`, each `--via` and `--to` by the route of least
# Σ step length × (deepest cell / cell depth)², which keeps to the deep channel rather than the
# shortest line. Each variable (default: T and S, those in the files) gets two rows: the top
# `--upper` metres (default 50) and the full depth, with the bottom along the path in black. A map
# of the path sits beside them. The colour range is the 1st–99th percentile of the wet values over
# the files, fixed for the whole video so the seasons compare.
#
# The files carry everything else: the model's bottom and land mask, the levels, and the start date
# the seconds of model time count from. The video goes beside the first file unless `--out` is given.
module SectionVideos

using Dates: DateTime, Millisecond, format
using NCDatasets: NCDataset
using Statistics: quantile
import CairoMakie as M

const FRAMES_PER_SECOND = 10
const COLOUR_RANGE_RECORDS = 50
const EARTH_RADIUS_KM = 6371.0
const CONTOUR_LEVELS = 10

const LABELS = Dict("T" => ("Temperature (°C)", :thermal), "S" => ("Salinity (g kg⁻¹)", :haline))

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

function section_path(grid, points)
    cells = [nearest_wet_cell(grid, point) for point in points]
    path = [cells[1]]
    for (start, stop) in zip(cells[1:end-1], cells[2:end])
        append!(path, deepest_route(grid, start, stop)[2:end])
    end
    distance = cumsum([0; [kilometres(centre(grid, a), centre(grid, b)) for (a, b) in zip(path[1:end-1], path[2:end])]])
    return (; cells = path, distance, depth = [grid.depth[I] for I in path])
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
    for file in paths
        NCDataset(file) do ds
            values = ds[variable]
            for record in eachindex(ds["time"])
                field = values[:, :, :, record]
                section = Float32[field[I[1], I[2], end+1-k] for I in path.cells, k in eachindex(centres)]
                section[below_bottom] .= NaN
                push!(columns, section)
            end
            append!(times, ds["time"][:])
        end
    end
    return (; values = cat(columns...; dims = 3), times, faces, centres)
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

function record_video(file, sections, path, grid, title_prefix, first_date, upper)
    edges = [path.distance[1]; (path.distance[1:end-1] .+ path.distance[2:end]) ./ 2; path.distance[end]]
    times = last(first(sections)).times
    date(record) = first_date + Millisecond(round(Int, 1000 * times[record]))
    title = M.Observable("")
    deepest = maximum(path.depth)

    figure = M.Figure(size = (1900, 520 * length(sections) + 80), fontsize = 18)
    M.Label(figure[0, 1:3], title; fontsize = 22, font = :bold)
    frames = Dict{String,Any}()
    for (row, (variable, section)) in enumerate(sections)
        label, colormap = get(LABELS, variable, (variable, :viridis))
        colorrange = colour_range(section.values)
        levels = M.Makie.get_tickvalues(M.LinearTicks(CONTOUR_LEVELS), colorrange...)
        frame = M.Observable(section.values[:, :, 1])
        frames[variable] = frame
        layout = figure[row, 1] = M.GridLayout()
        local plot
        for (panel, (bottom, name)) in enumerate(((upper, "top $(round(Int, upper)) m"), (deepest, "full depth")))
            axis = M.Axis(layout[panel, 1]; title = "$label, $name", yreversed = true, ylabel = "Depth (m)",
                          xlabel = panel == 2 && row == length(sections) ? "Distance along the path (km)" : "")
            plot = M.heatmap!(axis, edges, section.faces, frame; colormap, colorrange, nan_color = :gray80)
            M.contour!(axis, path.distance, section.centres, frame; levels, color = :black, linewidth = 0.7)
            M.lines!(axis, path.distance, path.depth; color = :black, linewidth = 1.5)
            M.limits!(axis, first(edges), last(edges), bottom, 0)
        end
        M.rowsize!(layout, 1, M.Relative(0.4))
        # The contours are at the colour bar's ticks.
        M.Colorbar(figure[row, 2], plot; label, ticks = levels)
    end
    map_panel!(figure[1:length(sections), 3], grid, path)
    M.colsize!(figure.layout, 3, M.Relative(0.22))

    M.record(figure, file, eachindex(times); framerate = FRAMES_PER_SECOND) do record
        for (variable, section) in sections
            frames[variable][] = section.values[:, :, record]
        end
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
    isempty(files) && error("Give the snapshot files")
    files = sort(files; by = part_number)

    grid = snapshot_grid(first(files))
    first_date = NCDataset(ds -> DateTime(ds.attrib["start_date"]), first(files))
    points = [parse_point(option("--from")); parse_point.(get(options, "--via", String[])); parse_point(option("--to"))]
    path = section_path(grid, points)
    println("Path: $(length(path.cells)) cells, $(round(last(path.distance); digits = 1)) km, " *
            "deepest $(round(maximum(path.depth); digits = 1)) m")

    in_files = NCDataset(ds -> keys(ds), first(files))
    variables = something(option("--variables"), join(intersect(["T", "S"], in_files), ","))
    sections = [variable => read_section(files, variable, path) for variable in split(variables, ",")]

    upper = parse(Float64, something(option("--upper"), "50"))
    out = something(option("--out"), dirname(abspath(first(files))))
    mkpath(out)
    stem = replace(first(splitext(basename(first(files)))), r"_part\d+$" => "")
    title_prefix = "$(basename(dirname(abspath(first(files))))) section"
    record_video(joinpath(out, "$(stem)_section_$(replace(variables, "," => "_")).mp4"), sections, path, grid,
                 title_prefix, first_date, upper)
    return nothing
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    SectionVideos.main(ARGS)
end
