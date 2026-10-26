"""
    Validation

Model-observation comparison: fetch observations, score the run's station output against them, and
draw the figures.

The pipeline modules before this one all run *before* a simulation; this is the only one that runs
after. It is deliberately the same shape all the same — an abstract config supertype per source, one
subtype per dataset, and hook methods dispatching on it, exactly as `docs/adding-a-source.md`
describes for bathymetry, forcing, rivers, boundaries and atmosphere. Adding an observation source is
therefore the same two-step it is everywhere else: subtype `AbstractObservationConfig` and overload
its hooks.

Four files below it:

- `observations.jl` — the sources, and one common record type they all return.
- `metrics.jl` — the statistics, over bare vectors. Knows nothing about fjords or files.
- `results.jl` — reading a run's station output back, as the inverse of `StationWriter`.
- `plots.jl` — the figures, one per report figure where there is a counterpart.
"""
module Validation

export validate_simulation

using CairoMakie
using Dates: DateTime, Date
import Dates
using Downloads
using LinearAlgebra
using Printf: @sprintf
using JLD2
using NCDatasets
using Statistics: cor, mean, std, var

using ..Configs:
    AbstractRunConfig,
    AbstractFjordConfig,
    bathymetry_config,
    run_config,
    results_root,
    writers,
    start_date
using ..Utils: coverage_window
using ..Bathymetry: bathymetry_path
using ..Simulations:
    Station,
    station_output,
    station_tag,
    RESULTS_START_DATE_ATTRIBUTE

include("observations.jl")
include("metrics.jl")
include("results.jl")
include("plots.jl")

"""
    validation_directory(config)

Where a run's validation output goes: `validation/` under the run's `results_root`, beside
the station files it scores.
"""
validation_directory(config::AbstractRunConfig) =
    joinpath(results_root(config), "validation")

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
    candidates = [path for (found, path) in station_files(directory, stem) if found == tag]
    isempty(candidates) && return nothing

    path = first(candidates)
    return if endswith(path, ".jld2")
        read_station_jld2(path, String(name), start_date(config))
    else
        read_station_netcdf(path, String(name))
    end
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

The observation sources a setup can be validated against without anyone supplying a file.

Exactly one today: Kartverket water level, at every station any `FieldStationWriter` in the setup
writes the free surface at. The rule is by variable rather than by writer name — a writer that
records `η` at a point is recording water level, and Kartverket is the public source of observed
water level — so a setup gets this for free by naming its gauges and nothing else.

Held programmes — current-meter moorings, monitoring CTD casts, fixed thermistors — are **not** here
and cannot be: they are held by their owners rather than published, so a source for them is a
reader over files someone supplies, such as a `CsvObservations`. Pass those to `validate_simulation`
explicitly once their exports exist.

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
    return (KartverketSeaLevel(; data_root, stations = unique(g -> g.name, gauges)),)
end

"""
    validate_simulation(config; observations, window, depths)

Score a finished run against observations and write the tables and figures. Returns the
validation directory, or `nothing` for a setup naming no `run_config`. Only picks out the run and
bathymetry configs and hands them to `score_run`.
"""
validate_simulation(config::AbstractFjordConfig; kwargs...) =
    score_run(run_config(config), bathymetry_config(config); kwargs...)

"""
    score_run(run, bathymetry_config; observations, window, depths)

`validate_simulation`'s work, and the one pipeline step that runs *after* `run_simulation` rather than
before it. It fetches whatever each observation source can supply over the run's window, pairs it
with the run's station output, and writes a skill table, a tidal constituent table, and one figure
per station-variable pair into `validation_directory`.

Every stage is permissive about missing data on purpose. A validation sweep asks every source for
every station and variable it might hold, and most combinations are empty: an ADCP that was in the
water for ten weeks of a twenty-one month run, a CTD programme that visits a station seven times a
year, a station the grid could not place. Each of those is a blank row, not a failure, and a run
should be scorable on what exists rather than on everything existing.

**Stations are matched between a source and the run by `name`.** A source's station list and the
setup's station writers have to spell a station the same way — `"OF-1"`, not `"OF1"` — because
nothing else identifies them to each other. `default_observations` sidesteps this by building its
station list *from* the writers; a source you supply has to agree by hand, and a mismatch shows up
as "no model output for ..." rather than as an error.

# Arguments
- `observations`: the sources to score against, as `AbstractObservationConfig`s.
- `window`: `(start, stop)`; defaults to the run's own `coverage_window`.
- `depths`: which depths to score profiles at, in metres. Defaults to a spread through the upper
  ocean and the basin; a surface series ignores it.
"""
# A setup that is never simulated has nothing to validate, which is an opt-out rather than an error —
# the same `nothing` `run_simulation` returns for it, reported as a no-op by `CLI.main`.
score_run(::Nothing, bathymetry; kwargs...) = nothing

function score_run(
    run::AbstractRunConfig,
    bathymetry;
    observations = default_observations(run, bathymetry),
    window = coverage_window(run),
    depths = (0.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0),
)
    recorders = station_writers(run)
    isempty(recorders) && throw(
        ArgumentError(
            "This setup's simulation writes no station output, so there is nothing to validate " *
            "against. Add a `StationWriter` or `FieldStationWriter` to its `writers`.",
        ),
    )

    directory = validation_directory(run)
    mkpath(directory)
    @info "Validating $(basename(results_root(run))) over $(window[1]) .. $(window[2])"

    for source in observations
        download_observations(source, window)
    end

    scores = Pair{String,SkillMetrics}[]
    tidal_fits = Pair{String,HarmonicFit}[]

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
                @info "  no model output for $(station.name) $variable"
                continue
            end

            score_pair!(scores, tidal_fits, directory, observed, modelled, depths)
        end
    end

    write_validation_summary(directory, scores, tidal_fits)
    @info "Validation written to $directory"
    return directory
end

"""
    score_pair!(scores, tidal_fits, directory, observed, modelled, depths)

Score one observation series against one model series at every depth the observations reach, append
the results, and draw the figures that pair supports.

A surface series is scored once; a profile is scored at each of `depths` that lies inside both the
instrument's range and the model column, so a shallow station contributes the levels it has rather
than a column of blanks.
"""
function score_pair!(scores, tidal_fits, directory, observed, modelled, depths)
    label = "$(observed.station) $(observed.variable)"
    surface = length(observed.depths) == 1 && length(modelled.depths) == 1

    targets = if surface
        (nothing,)
    else
        [d for d in depths if d <= min(maximum(observed.depths), maximum(modelled.depths))]
    end

    for depth in targets
        times, o, m = match_series(observed, modelled; depth)
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
                joinpath(directory, "timeseries_$(observed.station).png"),
                observed,
                modelled;
                title = observed.station,
            )
        end

        plot_qq_scatter(
            joinpath(directory, "qq_$(replace(name, r"[^A-Za-z0-9]+" => "_")).png"),
            o,
            m;
            title = name,
            unit = observed.units,
        )
    end

    return scores
end

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

    open(joinpath(directory, "skill.txt"), "w") do io
        write(io, skill_table(scores))
    end

    plot_taylor(joinpath(directory, "taylor.png"), scores)
    plot_target(joinpath(directory, "target.png"), scores)

    if !isempty(tidal_fits)
        open(joinpath(directory, "tides.txt"), "w") do io
            write(io, tidal_constituent_table(tidal_fits))
        end
        # Two fits per station, observed then modelled; plot each station's pair together.
        for index = 1:2:(length(tidal_fits) - 1)
            station = first(split(tidal_fits[index].first, " observed"))
            plot_tidal_constituents(
                joinpath(directory, "constituents_$(station).png"),
                tidal_fits[index:(index + 1)];
                title = station,
            )
        end
    end

    return directory
end

end # module
