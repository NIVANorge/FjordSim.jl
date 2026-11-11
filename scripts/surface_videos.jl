# Videos of a run's surface snapshots (`SurfaceSnapshotWriter`), one MP4 per variable, for one
# file — one monthly part, for a writer with `split_times`:
#
#     julia --project scripts/surface_videos.jl --config examples/mjosa/mjosa.jl \
#         ~/FjordSim_results/mjosa/snapshots_surface_<tag>_part1.jld2 \
#         [--sea-ice ~/FjordSim_results/mjosa/snapshots_sea_ice_<tag>_part1.jld2] [--out DIR]
#
# The config gives the start date (the files store seconds of model time), the bathymetry (the land
# mask and coordinates) and the title. Each video is made when the files hold its variables:
# temperature, salinity, surface elevation (as the departure from the mean over the water, which is
# in the title, so a tide or a regulated lake level does not hide the set-up), current speed with
# arrows, and from `--sea-ice` (`SeaIceSnapshotWriter`, same records) ice thickness, concentration
# and snow depth. Any other variable gets a plain video of its own.
#
# Land is grey and a NaN in the water is magenta. A video's colour range is the 1st–99th percentile
# of its wet values over the file, except the ice videos' fixed ranges, which let months compare.
# Videos go beside the surface file unless `--out` is given.
#
# A module, because `fjord_config` evaluates the setup file in `Main`, whose globals (`start_date`,
# ...) would clash with the names imported here.
module SurfaceVideos

using Dates: DateTime, Millisecond, format
using FjordSim: fjord_config
using FjordSim.Configs: bathymetry_config, run_config, start_date
using FjordSim.Bathymetry: bathymetry_path
using JLD2: jldopen
using NCDatasets: NCDataset
using Statistics: mean, quantile
import CairoMakie as M

const FRAMES_PER_SECOND = 24
const COLOUR_RANGE_RECORDS = 50
# Current arrows along the longer side of the grid. An arrow at the top of the colour range is one
# arrow spacing long.
const ARROWS_ACROSS = 40

surface_field(files, name, iteration) = files.surface["timeseries/$name/$iteration"][:, :, 1]
sea_ice_field(files, name, iteration) = files.sea_ice["timeseries/$name/$iteration"][:, :, 1]

# `u` and `v` sit on the cell faces; their means over each cell's two faces are on its centre.
function centred_velocities(files, iteration)
    u = surface_field(files, "u", iteration)
    v = surface_field(files, "v", iteration)
    return (u[1:end-1, :] .+ u[2:end, :]) ./ 2, (v[:, 1:end-1] .+ v[:, 2:end]) ./ 2
end

function current_speed(files, iteration, wet)
    u, v = centred_velocities(files, iteration)
    return hypot.(u, v)
end

mean_elevation(files, iteration, wet) = mean(surface_field(files, "η", iteration)[wet])

function elevation_anomaly(files, iteration, wet)
    η = surface_field(files, "η", iteration)
    return 100 .* (η .- mean(η[wet]))
end

function elevation_suffix(files, iteration, wet)
    centimetres = round(100 * mean_elevation(files, iteration, wet); digits = 1)
    return ", mean η = $centimetres cm"
end

# `h` and `hs` are the thicknesses over the ice-covered fraction `ℵ`.
ice_thickness(files, iteration, wet) = sea_ice_field(files, "h", iteration) .* sea_ice_field(files, "ℵ", iteration)
snow_depth(files, iteration, wet) = sea_ice_field(files, "hs", iteration) .* sea_ice_field(files, "ℵ", iteration)
ice_concentration(files, iteration, wet) = sea_ice_field(files, "ℵ", iteration)

no_suffix(files, iteration, wet) = ""

surface_video(variable; name = variable, label = variable, colormap = :viridis) =
    (; name, label, colormap, colorrange = nothing, symmetric = false, arrows = false, suffix = no_suffix,
     values = (files, iteration, wet) -> surface_field(files, variable, iteration))

# Each video with the variables it reads from the surface and sea-ice files. `colorrange = nothing`
# takes the percentiles of the file; `symmetric` centres it on zero.
const VIDEOS = (
    (surface = ("T",), sea_ice = (),
     video = surface_video("T"; name = "temperature", label = "Temperature (°C)", colormap = :thermal)),
    (surface = ("S",), sea_ice = (),
     video = surface_video("S"; name = "salinity", label = "Salinity (g kg⁻¹)", colormap = :haline)),
    (surface = ("η",), sea_ice = (),
     video = (name = "elevation", label = "Surface elevation − mean (cm)", colormap = :balance, colorrange = nothing,
              symmetric = true, arrows = false, suffix = elevation_suffix, values = elevation_anomaly)),
    (surface = ("u", "v"), sea_ice = (),
     video = (name = "currents", label = "Current speed (m s⁻¹)", colormap = :speed, colorrange = nothing,
              symmetric = false, arrows = true, suffix = no_suffix, values = current_speed)),
    (surface = (), sea_ice = ("h", "ℵ"),
     video = (name = "ice_thickness", label = "Ice thickness × concentration (m)", colormap = :ice,
              colorrange = (0, 1), symmetric = false, arrows = false, suffix = no_suffix, values = ice_thickness)),
    (surface = (), sea_ice = ("ℵ",),
     video = (name = "ice_concentration", label = "Ice concentration", colormap = :ice, colorrange = (0, 1),
              symmetric = false, arrows = false, suffix = no_suffix, values = ice_concentration)),
    (surface = (), sea_ice = ("hs", "ℵ"),
     video = (name = "snow", label = "Snow depth × concentration (m)", colormap = :ice, colorrange = (0, 0.3),
              symmetric = false, arrows = false, suffix = no_suffix, values = snow_depth)),
)

variable_names(file) = setdiff(keys(file["timeseries"]), ["t"])
variable_names(::Nothing) = String[]

"""
    available_videos(files)

The `VIDEOS` whose variables the files hold, then a plain video of every other surface variable.
"""
function available_videos(files)
    surface_names = variable_names(files.surface)
    sea_ice_names = variable_names(files.sea_ice)
    known = filter(entry -> entry.surface ⊆ surface_names && entry.sea_ice ⊆ sea_ice_names, VIDEOS)
    read_names = union((entry.surface for entry in VIDEOS)...)
    others = [surface_video(name) for name in setdiff(surface_names, read_names)]
    return [[entry.video for entry in known]; others]
end

function setup_grid(config)
    longitude, latitude, depth = NCDataset(bathymetry_path(bathymetry_config(config))) do ds
        ds["lon"][:], ds["lat"][:], ds["h"][:, :]
    end
    return (; longitude, latitude, land = depth .>= 0)
end

function percentile_range(files, iterations, wet, values)
    records = unique(round.(Int, range(1, length(iterations); length = COLOUR_RANGE_RECORDS)))
    samples = Float32[]
    for record in records
        append!(samples, filter(isfinite, values(files, iterations[record], wet)[wet]))
    end
    isempty(samples) && return (0.0f0, 1.0f0)
    low, high = quantile(samples, [0.01, 0.99])
    return low < high ? (low, high) : (low - 1, high + 1)
end

function arrow_layer(grid, top_speed)
    (; longitude, latitude, land) = grid
    stride = max(1, round(Int, maximum(size(land)) / ARROWS_ACROSS))
    cells = [I for I in CartesianIndices(land)[1:stride:end, 1:stride:end] if !land[I]]
    positions = [M.Point2f(longitude[I[1]], latitude[I[2]]) for I in cells]
    # Degrees of latitude per m s⁻¹.
    arrow_length = stride * (last(latitude) - first(latitude)) / (length(latitude) - 1) / top_speed
    # East in degrees of longitude, so the arrow points the right way on the aspect-corrected axis.
    directions(u, v) = [arrow_length * M.Vec2f(u[I] / cosd(latitude[I[2]]), v[I]) for I in cells]
    return positions, directions
end

"""
    record_video(path, files, iterations, times, grid, video, title_prefix, first_date)

Write one MP4 of `video` over every record of the files.
"""
function record_video(path, files, iterations, times, grid, video, title_prefix, first_date)
    (; longitude, latitude, land) = grid
    wet = .!land
    frame_title(record) = title_prefix * ", " *
        format(first_date + Millisecond(round(Int, 1000 * times[record])), "yyyy-mm-dd HH:MM") *
        video.suffix(files, iterations[record], wet)

    values = M.Observable(video.values(files, first(iterations), wet))
    title = M.Observable(frame_title(1))
    colorrange = something(video.colorrange, percentile_range(files, iterations, wet, video.values))
    video.symmetric && (colorrange = (-1, 1) .* maximum(abs, colorrange))

    aspect = (last(longitude) - first(longitude)) * cosd((first(latitude) + last(latitude)) / 2) /
             (last(latitude) - first(latitude))
    figure = M.Figure(size = aspect < 1 ? (800, 1000) : (1100, 800))
    axis = M.Axis(figure[1, 1]; title, xlabel = "Longitude (°E)", ylabel = "Latitude (°N)",
                  aspect = M.AxisAspect(aspect))
    # Fixed to the grid, so an arrow at the edge is clipped rather than widening the axis.
    M.limits!(axis, extrema(longitude)..., extrema(latitude)...)
    field = M.heatmap!(axis, longitude, latitude, values; colormap = video.colormap, colorrange,
                       nan_color = :magenta)
    M.heatmap!(axis, longitude, latitude, ifelse.(land, 1.0f0, NaN32); colormap = [:gray70, :gray70],
               colorrange = (0, 2), nan_color = :transparent)
    M.Colorbar(figure[1, 2], field; label = video.label)

    if video.arrows
        positions, directions = arrow_layer(grid, last(colorrange))
        arrows = M.Observable(directions(centred_velocities(files, first(iterations))...))
        M.arrows2d!(axis, positions, arrows; color = :black, shaftwidth = 1, tipwidth = 5, tiplength = 5)
    end

    M.record(figure, path, eachindex(iterations); framerate = FRAMES_PER_SECOND) do record
        values[] = video.values(files, iterations[record], wet)
        title[] = frame_title(record)
        video.arrows && (arrows[] = directions(centred_velocities(files, iterations[record])...))
    end

    println("Wrote $path")
    return path
end

record_iterations(file) = sort(parse.(Int, keys(file["timeseries/t"])))

function check_records(iterations, sea_ice, path)
    record_iterations(sea_ice) == iterations || error("$path does not hold the surface file's records")
    return nothing
end
check_records(iterations, ::Nothing, path) = nothing

open_sea_ice(f, ::Nothing) = f(nothing)
open_sea_ice(f, path) = jldopen(f, path)

function main(args)
    value(flag, default) = (i = findfirst(==(flag), args); isnothing(i) ? default : args[i+1])
    config_path = value("--config", nothing)
    isnothing(config_path) && error("--config <setup.jl> is required")
    flagged = [args[i+1] for i in eachindex(args) if startswith(args[i], "--")]
    path = only(filter(arg -> !startswith(arg, "--") && arg ∉ flagged, args))
    sea_ice_path = value("--sea-ice", nothing)
    out = value("--out", dirname(abspath(path)))
    mkpath(out)

    config = fjord_config(config_path)
    # `invokelatest`: methods the setup file adds to a hook are newer than this call.
    grid = Base.invokelatest(setup_grid, config)
    first_date = start_date(run_config(config))
    title_prefix = "$(first(splitext(basename(config_path)))) surface"

    jldopen(path) do surface
        open_sea_ice(sea_ice_path) do sea_ice
            files = (; surface, sea_ice)
            iterations = record_iterations(surface)
            check_records(iterations, sea_ice, sea_ice_path)
            times = [surface["timeseries/t/$iteration"] for iteration in iterations]
            stem = first(splitext(basename(path)))
            for video in available_videos(files)
                record_video(joinpath(out, "$(stem)_$(video.name).mp4"), files, iterations, times, grid, video,
                             title_prefix, first_date)
            end
        end
    end
    return nothing
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    SurfaceVideos.main(ARGS)
end
