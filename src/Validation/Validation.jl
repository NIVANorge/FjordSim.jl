"""
    Validation

Model-observation comparison: score the run's station output against the observations the
`Observations` module reads, and draw the figures.

The pipeline modules before this one all run *before* a simulation; this is the only one that runs
after. The observations come from the setup's `observation_configs`, which `download_observations`
fetched; adding a source is subtyping `Configs.AbstractObservationConfig`, not an edit here.

Five files below it:

- `metrics.jl` — the statistics, over bare vectors. Knows nothing about fjords or files.
- `results.jl` — reading a run's station output back, as the inverse of `StationWriter`.
- `snapshots.jl` — the columns and the domain statistics read from a run's 3D snapshots.
- `pooled.jl` — every observed sample beside the model, pooled into tables by depth band and group.
- `plots.jl` — the figures, one per report figure where there is a counterpart.
"""
module Validation

export validate_simulation

using CairoMakie
using Dates: DateTime, Date, Day
import Dates
using LinearAlgebra
using Printf: @sprintf
using JLD2
using NCDatasets
using Statistics: cor, mean, median, std, var

using ..Configs:
    AbstractRunConfig,
    AbstractFjordConfig,
    AbstractRiverConfig,
    bathymetry_config,
    boundary_data_config,
    river_config,
    observation_configs,
    validation_config,
    run_config,
    results_root,
    writers,
    start_date
using ..Utils: coverage_window
using ..Bathymetry: bathymetry_path
using ..Grids: open_edges
using ..Rivers: river_forcing_path
using ..Observations:
    CsvObservations,
    KartverketSeaLevel,
    ObservationSeries,
    Station,
    fetch_observations,
    interpolate_to_depths,
    observation_series,
    observation_stations,
    station_tag,
    VannmiljoObservations
using ..Simulations: SnapshotWriter, station_output, RESULTS_START_DATE_ATTRIBUTE

include("metrics.jl")
include("results.jl")

"""
    ValidationConfig(; spin_up, depth_bands, station_groups, search_radius, plausible_ranges,
                     stratification)

What `validate_simulation` does beyond scoring each station from its station writer's output. Every
field is a scientific choice the setup states; there are no defaults.

- An observation at a station no station writer records is scored from the run's `SnapshotWriter`
  file instead, at the column `snapshot_column` picks within `search_radius` cells.
- Every record of that file is summarised over the whole domain in `domain_extremes.csv`: the
  volume mean, the extremes and where they are, and how many wet cells lie outside
  `plausible_ranges`. An overshoot or a runaway cell shows there first.
- Every observed sample of a profile is paired with the model at its own time and depth
  (`samples.csv`) and pooled into bias and RMSE by depth band — for each of `station_groups`, each
  station and each month — and each cast into the difference `stratification` names. The pooled
  tables leave out the first `spin_up` of the window, the model's adjustment from its initial
  state.

# Fields
- `spin_up`: how long after the window opens the pooled tables start.
- `depth_bands`: `(upper, lower)` depth ranges in metres, positive down; a sample at depth `d` is in
  the band `upper <= d < lower`, so `(50, Inf)` is everything from 50 m down. Each table also has
  the band `"all"`.
- `station_groups`: `name => station names`, in table order. A station can be in several groups or
  none; the group `"all"` holds every station.
- `search_radius`: how far from a station, in grid cells, its snapshot column may be.
- `plausible_ranges`: `(low, high)` per variable name (`"T"`, `"S"`); a wet cell below `low` or above
  `high` is counted.
- `stratification`: `(variable, upper, lower)`: the cast's `variable` at depth `lower` minus at
  depth `upper`, in metres.
"""
Base.@kwdef struct ValidationConfig
    spin_up::Day
    depth_bands::Vector{Tuple{Float64,Float64}}
    station_groups::Vector{Pair{String,Vector{String}}}
    search_radius::Int
    plausible_ranges::Dict{String,Tuple{Float64,Float64}}
    stratification::Tuple{String,Float64,Float64}
end

include("snapshots.jl")
include("pooled.jl")
include("plots.jl")

"""
    validation_directory(config)

Where a run's validation output goes: `validation/` under the run's `results_root`, beside
the station files it scores. The tables go in its `data/` (`data_file`), the figures in its `plots/`
(`plot_file`).
"""
validation_directory(config::AbstractRunConfig) =
    joinpath(results_root(config), "validation")

data_file(directory, name) = joinpath(directory, "data", name)
plot_file(directory, name) = joinpath(directory, "plots", name)

"""
    station_writers(config)

The writers a run config names whose `station_output` is `true` — the built-in `StationWriter`s and
`FieldStationWriter`s — which is what says where a run's station output is and what is in it.
Reading it from the config rather than from the directory means the validation and the run cannot
disagree about which variables exist.
"""
station_writers(config::AbstractRunConfig) = [writer for writer in writers(config) if station_output(writer)]

"""
    model_series(config, writer, station, variable)

The run's output for one station and variable, or `nothing` when the run did not place that station
— which happens for real reasons (a site outside the grid, or with no water in reach) and is not an
error.

`FieldStationWriter` output is JLD2 and carries no attributes, so the calendar comes from the
run config; `StationWriter` output is NetCDF and carries its own.

`variable` is the observation's name for it, which `model_variable` turns into the model field the
writer recorded.
"""
function model_series(
    config::AbstractRunConfig,
    writer,
    station::Station,
    variable,
)
    name = model_variable(variable)
    name in writer.variables || return nothing

    directory = results_root(config)
    stem = first(splitext(writer.output_file))
    tag = station_tag(station)
    paths = [path for (found, path) in station_files(directory, stem) if found == tag]
    isempty(paths) && return nothing

    read_file(path) = endswith(path, ".jld2") ? read_station_jld2(path, String(name), start_date(config)) :
                      read_station_netcdf(path, String(name))
    return join_station_series(map(read_file, paths))
end

# The model field an observed variable is read from, where the two names differ. Observation sources
# spell the free surface `"eta"`, and the model's field is `η`.
const MODEL_VARIABLES = Dict("eta" => :η)

model_variable(variable) = get(MODEL_VARIABLES, String(variable), Symbol(variable))

# Which model field each observed variable is scored against. `eta` is the free surface, which the
# station writers put in JLD2 for the reason `FieldSnapshotWriter` documents.
const VALIDATION_VARIABLES = ("eta", "T", "S", "u", "v")

"""
    default_observations(run_config, bathymetry_config)

The observation sources a setup is validated against without naming them.

Exactly one today: Kartverket water level over the run's `coverage_window`, at every station any
`FieldStationWriter` in the setup writes the free surface at. The rule is by variable rather than
by writer name — a writer that records `η` at a point is recording water level, and Kartverket is
the public source of observed water level — so a setup gets this for free by naming its gauges and
nothing else. `validate_simulation` fetches these itself, since no setup names them.

Every other source — CTD casts, current-meter moorings, fixed thermistors — a setup names in its
`observation_configs`, which `download_observations` fetches and `validate_simulation` scores after
these.

Returns an empty tuple for a setup with no free-surface station writer, which is not an error —
`validate_simulation` will simply find nothing to score and say so.
"""
function default_observations(run::AbstractRunConfig, bathymetry)
    gauges = Station[]
    # Only a `FieldStationWriter` can record `η`: a `StationWriter` refuses it when it is attached.
    for writer in station_writers(run)
        :η in writer.variables || continue
        append!(gauges, writer.stations)
    end

    isempty(gauges) && return ()

    # Beside the processed bathymetry, the one file every setup has, so the cache lands with its
    # other data.
    data_root = dirname(bathymetry_path(bathymetry))
    start_date, stop_date = coverage_window(run)
    return (KartverketSeaLevel(; data_root, stations = unique(g -> g.name, gauges), start_date, stop_date),)
end

"""
    validate_simulation(config; observations, window, depths)

Score a finished run against observations and write the tables and figures. Returns the
validation directory, or `nothing` for a setup naming no `run_config`. Only picks out the run,
bathymetry, river, open-boundary, observation and validation configs and hands them to `score_run`.
"""
validate_simulation(config::AbstractFjordConfig; kwargs...) = score_run(
    run_config(config),
    bathymetry_config(config),
    river_config(config),
    boundary_data_config(config),
    observation_configs(config),
    validation_config(config);
    kwargs...,
)

"""
    score_run(run, bathymetry_config, rivers, boundaries, sources, validation; gauges, observations, window, depths)

`validate_simulation`'s work, and the one pipeline step that runs *after* `run_simulation` rather than
before it. It reads what each observation source has on disk over the run's window — fetching only
the `default_observations` tide gauges, which no setup names — pairs it with the run's station
output, and writes a skill table, a tidal constituent table, the tide-gauge figures and the station
figures (`write_station_figures`) and a map of the stations, river outlets and open edges
(`write_station_map`) into `validation_directory`. With a `ValidationConfig` as `validation`,
it also scores what no station writer recorded from the 3D snapshots, summarises the domain, and
pools every sample into the depth-band tables (`ValidationConfig`).

Every stage is permissive about missing data on purpose. A validation sweep asks every source for
every station and variable it might hold, and most combinations are empty: an ADCP that was in the
water for ten weeks of a twenty-one month run, a CTD programme that visits a station seven times a
year, a station the grid could not place. Each of those is a blank row, not a failure, and a run
should be scorable on what exists rather than on everything existing.

**Stations are matched between a source and the run by `name`.** A source's station list and the
setup's station writers have to spell a station the same way — `"OF-1"`, not `"OF1"` — because
nothing else identifies them to each other. `default_observations` sidesteps this by building its
station list *from* the writers; a source you supply has to agree by hand, and a mismatch shows up
as "no model output for ..." rather than as an error. Only an observation scored from the 3D
snapshots (`ValidationConfig`) is placed by its station's position instead.

# Arguments
- `rivers`: the setup's `river_config`, whose outlets the map shows, or `nothing`.
- `boundaries`: the setup's `boundary_data_config`, whose `open_edges` the map shows, or `nothing`.
- `sources`: the setup's own observation sources, its `observation_configs`.
- `validation`: the setup's `validation_config`, a `ValidationConfig` or `nothing`.
- `gauges`: the sources `score_run` fetches itself; defaults to `default_observations`.
- `observations`: the sources to score against, as `AbstractObservationConfig`s. Defaults to
  `gauges` followed by `sources`.
- `window`: `(start, stop)`; defaults to the run's own `coverage_window`.
- `depths`: which depths to score profiles at, in metres. Defaults to a spread through the upper
  ocean and the basin; a surface series ignores it.
"""
# A setup that is never simulated has nothing to validate, which is an opt-out rather than an error —
# the same `nothing` `run_simulation` returns for it, reported as a no-op by `CLI.main`.
score_run(::Nothing, bathymetry, rivers, boundaries, sources, validation; kwargs...) = nothing

function score_run(
    run::AbstractRunConfig,
    bathymetry,
    rivers,
    boundaries,
    sources,
    validation;
    gauges = default_observations(run, bathymetry),
    observations = (gauges..., sources...),
    window = coverage_window(run),
    depths = (0.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0),
)
    recorders = station_writers(run)
    isempty(recorders) && isnothing(validation) && throw(
        ArgumentError(
            "This setup's simulation writes no station output, so there is nothing to validate " *
            "against. Add a `StationWriter` or `FieldStationWriter` to its `writers`, or a " *
            "`validation_config` to score its `SnapshotWriter` output.",
        ),
    )

    directory = validation_directory(run)
    mkpath(dirname(data_file(directory, "")))
    mkpath(dirname(plot_file(directory, "")))
    @info "Validating $(basename(results_root(run))) over $(window[1]) .. $(window[2])"

    # The gauges are named, so they need no grid.
    for gauge in gauges
        fetch_observations(nothing, gauge)
    end

    scores = Pair{String,SkillMetrics}[]
    tidal_fits = Pair{String,HarmonicFit}[]
    samples = SamplePair[]
    unrecorded = Tuple{Station,ObservationSeries}[]
    paired = Tuple{ObservationSeries,ModelSeries,Dates.Millisecond}[]

    for source in observations, station in observation_stations(source)
        for variable in VALIDATION_VARIABLES
            observed = observation_series(source, station, variable; window)
            isnothing(observed) && continue

            modelled = nothing
            for writer in recorders
                any(s -> s.name == station.name, writer.stations) || continue
                modelled = model_series(run, writer, station, variable)
                isnothing(modelled) || break
            end

            if isnothing(modelled)
                push!(unrecorded, (station, observed))
                continue
            end

            score_pair!(scores, tidal_fits, directory, observed, modelled, depths)
            append!(samples, sample_pairs(observed, modelled))
            push!(paired, (observed, modelled, Dates.Millisecond(Dates.Hour(1))))
        end
    end

    for (observed, modelled, tolerance) in snapshot_series(run, validation, unrecorded, directory)
        score_pair!(scores, tidal_fits, directory, observed, modelled, depths; tolerance)
        append!(samples, sample_pairs(observed, modelled; tolerance))
        push!(paired, (observed, modelled, tolerance))
    end

    write_validation_summary(directory, scores, tidal_fits)
    write_pooled_tables(directory, samples, window, validation)
    write_station_figures(directory, paired, window, spin_up_end(window, validation))
    write_station_map(directory, bathymetry, rivers, open_edges(boundaries), observations, Set(pair[1].station for pair in paired))
    @info "Validation written to $directory"
    return directory
end

"""
    score_pair!(scores, tidal_fits, directory, observed, modelled, depths; tolerance = Dates.Hour(1))

Score one observation series against one model series at every depth the observations reach, append
the results, and, for a tide gauge, draw its figures.

A surface series is scored once; a profile is scored at each of `depths` that lies inside both the
instrument's range and the model column, so a shallow station contributes the levels it has rather
than a column of blanks. `tolerance` is how far from an observation the nearest model record may
be, as for `match_series`.
"""
function score_pair!(
    scores,
    tidal_fits,
    directory,
    observed,
    modelled,
    depths;
    tolerance = Dates.Hour(1),
)
    label = "$(observed.station) $(observed.variable)"
    surface = length(observed.depths) == 1 && length(modelled.depths) == 1

    targets = if surface
        (nothing,)
    else
        [d for d in depths if d <= min(maximum(observed.depths), maximum(modelled.depths))]
    end

    for depth in targets
        times, o, m = match_series(observed, modelled; depth, tolerance)
        metrics = skill_metrics(o, m)
        metrics.count == 0 && continue

        name = isnothing(depth) ? label : "$label @ $(Int(round(depth))) m"
        push!(scores, name => metrics)

        # A tide gauge gets the harmonic treatment; nothing else does, because a constituent fit of
        # a temperature record is meaningless and would fill the table with noise.
        if observed.variable == "eta" && elapsed_hours(times) >= 30 * 24
            usable = [c for c in TIDAL_CONSTITUENTS if c.period_hours <= elapsed_hours(times)]
            reference = first(times)
            push!(
                tidal_fits,
                "$(observed.station) observed" =>
                    harmonic_fit(times, o; constituents = usable, reference),
            )
            push!(
                tidal_fits,
                "$(observed.station) modelled" =>
                    harmonic_fit(times, m; constituents = usable, reference),
            )

            plot_timeseries_comparison(
                plot_file(directory, "timeseries_$(observed.station).png"),
                observed,
                modelled;
                title = observed.station,
            )
        end
    end

    return scores
end

"""
    spin_up_end(window, validation)

When the pooled tables start, `nothing` without a `ValidationConfig`.
"""
spin_up_end(window, ::Nothing) = nothing
spin_up_end(window, validation::ValidationConfig) = window[1] + validation.spin_up

# The depth of the surface drift figure, in metres: the shallowest every cast samples.
const SURFACE_DRIFT_DEPTH = 1.0
# Stations per drift figure.
const DRIFT_STATIONS_PER_FIGURE = 8

"""
    write_station_figures(directory, paired, window, spin_up_end)

The figures that follow each station through the run, from every `(observed, modelled, tolerance)` of
a temperature or salinity profile in `paired`:

- `profiles_<station>.png`: every cast beside the model (`plot_cast_profiles`).
- `surface_drift.png`: every station at 1 m through `window` (`plot_drift`).
- `basin_drift.png`: every station at its own deep level, `drift_depth`; a station shallower than
  5 m is left out.

A drift figure holds eight stations; more continue in `surface_drift_2.png`, and so on.
"""
function write_station_figures(directory, paired, window, spin_up_end)
    profiles = [pair for pair in paired if pair[1].variable in ("S", "T") && length(pair[2].depths) > 1]
    stations = [
        station => sort([pair for pair in profiles if pair[1].station == station]; by = pair -> pair[1].variable)
        for station in unique(pair[1].station for pair in profiles)
    ]
    isempty(stations) && return nothing

    for (station, pairs) in stations
        title = @sprintf("%s: deepest cast %.0f m, model column %.0f m", station,
            maximum(deepest_depth(pair[1]) for pair in pairs), maximum(deepest_depth(pair[2]) for pair in pairs))
        plot_cast_profiles(plot_file(directory, "profiles_$(file_stem(station)).png"), pairs; title, spin_up_end)
    end

    surface = [(station, SURFACE_DRIFT_DEPTH, pairs) for (station, pairs) in stations]
    basin = [(station, drift_depth(pairs), pairs) for (station, pairs) in stations]
    write_drift_figures(directory, "surface_drift", surface, window)
    write_drift_figures(directory, "basin_drift", [entry for entry in basin if entry[2] >= 5], window)
    return directory
end

"""
    drift_depth(pairs)

A station's deep level for `basin_drift.png`: the median depth its casts reach, no deeper than its
model column, rounded down to 5 m.
"""
function drift_depth(pairs)
    depth = Inf
    for (observed, modelled, _) in pairs
        reached = [
            maximum((observed.depths[k] for k in axes(observed.values, 2) if isfinite(observed.values[n, k])); init = 0.0)
            for n in axes(observed.values, 1)
        ]
        depth = min(depth, median(reached), deepest_depth(modelled))
    end
    return 5 * fld(depth, 5)
end

function write_drift_figures(directory, name, entries, window)
    pages = collect(Iterators.partition(entries, DRIFT_STATIONS_PER_FIGURE))
    for (page, chunk) in enumerate(pages)
        file = page == 1 ? "$name.png" : "$(name)_$page.png"
        plot_drift(plot_file(directory, file), collect(chunk), window...)
    end
    return directory
end

file_stem(name) = replace(name, r"[^A-Za-z0-9]+" => "_")

"""
    write_station_map(directory, bathymetry, rivers, edges, sources, scored)

`map.png`: every station of `sources` on the model depth of the processed bathymetry, filled where
its name is in `scored`, with the outlets of `rivers` (`river_outlets`) and the open `edges`. Skipped, with a warning,
when the bathymetry file is missing.
"""
function write_station_map(directory, bathymetry, rivers, edges, sources, scored)
    path = bathymetry_path(bathymetry)
    isfile(path) || (@warn "No bathymetry file at $path, so no station map."; return nothing)
    longitude, latitude, depth = NCDataset(path) do ds
        bottom = Float64.(coalesce.(ds["h"][:, :], NaN))
        Float64.(ds["lon"][:]), Float64.(ds["lat"][:]), map(h -> h < 0 ? -h : NaN, bottom)
    end
    stations = [(map_label(source), observation_stations(source), scored) for source in sources]
    plot_station_map(plot_file(directory, "map.png"), longitude, latitude, depth, stations,
        river_outlets(rivers, longitude, latitude), edges)
    return directory
end

"""
    river_outlets(rivers, longitude, latitude)

`(name, longitude, latitude, mean discharge)` per outlet in the file `add_rivers` wrote, at the
centre of its column; empty without a river config or its file.
"""
river_outlets(::Nothing, longitude, latitude) = Tuple{String,Float64,Float64,Float64}[]

function river_outlets(rivers::AbstractRiverConfig, longitude, latitude)
    path = river_forcing_path(rivers)
    isfile(path) || return river_outlets(nothing, longitude, latitude)
    return NCDataset(path) do ds
        discharge = Float64.(coalesce.(ds["discharge"][:, :], NaN))
        [
            (String(name), longitude[i], latitude[j], finite_mean(view(discharge, n, :)))
            for (n, (name, i, j)) in enumerate(zip(ds["river_name"][:], ds["river_i"][:], ds["river_j"][:]))
        ]
    end
end

finite_mean(values) = (finite = filter(isfinite, values); isempty(finite) ? 0.0 : mean(finite))

# The legend label of an observation source on the map.
map_label(source) = string(nameof(typeof(source)))
map_label(::KartverketSeaLevel) = "Kartverket tide gauge"
map_label(source::Union{VannmiljoObservations,CsvObservations}) = source.source

"""
    write_validation_summary(directory, scores, tidal_fits)

Write the skill table, the tidal constituent table and the two summary diagrams.

The tables are plain text rather than a figure because they hold numbers a reader will want to quote
beside a published evaluation's own, such as METreport 11/2017's.
"""
function write_validation_summary(directory, scores, tidal_fits)
    if isempty(scores)
        @warn "Nothing was scored: no observation series overlapped the run's station output."
        return directory
    end

    open(data_file(directory, "skill.txt"), "w") do io
        write(io, skill_table(scores))
    end

    plot_taylor(plot_file(directory, "taylor.png"), scores)
    plot_target(plot_file(directory, "target.png"), scores)

    if !isempty(tidal_fits)
        open(data_file(directory, "tides.txt"), "w") do io
            write(io, tidal_constituent_table(tidal_fits))
        end
        # Two fits per station, observed then modelled; plot each station's pair together.
        for index = 1:2:(length(tidal_fits) - 1)
            station = first(split(tidal_fits[index].first, " observed"))
            plot_tidal_constituents(
                plot_file(directory, "constituents_$(station).png"),
                tidal_fits[index:(index + 1)];
                title = station,
            )
        end
    end

    return directory
end

end # module
