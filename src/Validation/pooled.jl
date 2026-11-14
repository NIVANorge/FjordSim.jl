# Every observed sample of a profile beside the model at its own time and depth, pooled into tables.
#
# The per-station table scores each station at a few fixed depths. A monitoring programme visits
# dozens of stations a few times a year each, and the question it answers is where in the water
# column and where in the fjord the model is wrong, so these tables pool all of them by depth band,
# by station group, by station and by month.

const SamplePair = @NamedTuple{
    station::String,
    variable::String,
    time::DateTime,
    depth::Float64,
    observed::Float64,
    modelled::Float64,
}

"""
    sample_pairs(observation, model; tolerance = Dates.Hour(1))

Every finite observed value of a profile beside the model at the same depth, from the model record
nearest in time (within `tolerance`, as for `match_series`). Empty for a model series one level deep,
such as the free surface: the pooled tables are about the water column.
"""
function sample_pairs(observation::ObservationSeries, model::ModelSeries; tolerance = Dates.Hour(1))
    samples = SamplePair[]
    length(model.depths) > 1 || return samples

    for (index, time) in enumerate(observation.times)
        record = nearest_record(model.times, time, tolerance)
        isnothing(record) && continue
        modelled = profile_at(model.depths, view(model.values, record, :), observation.depths)
        for (level, depth) in enumerate(observation.depths)
            value = observation.values[index, level]
            isfinite(value) && isfinite(modelled[level]) || continue
            push!(samples, (observation.station, observation.variable, time, depth, value, modelled[level]))
        end
    end

    return samples
end

"""
    profile_at(depths, profile, targets)

A model profile at `targets` depths: linear between level centres, `NaN` below the deepest wet one,
and the top level's value above the top level's centre, which is inside the top cell.
"""
function profile_at(depths, profile, targets)
    values = interpolate_to_depths(depths, profile, targets)
    for (index, target) in enumerate(targets)
        target < first(depths) && (values[index] = profile[1])
    end
    return values
end

"""
    write_pooled_tables(directory, samples, window, validation)

`samples.csv`, every pair; the bias and RMSE tables `skill_bands.csv` (by station group),
`skill_stations.csv` and `skill_months.csv`, all by depth band; and `stratification.csv` with its
summary `stratification_skill.csv`. Everything but `samples.csv` leaves out the first
`validation.spin_up` of `window`. Nothing to write without a `ValidationConfig`.
"""
write_pooled_tables(directory, samples, window, ::Nothing) = nothing

function write_pooled_tables(directory, samples, window, validation::ValidationConfig)
    write_samples(joinpath(directory, "samples.csv"), samples)

    start = window[1] + validation.spin_up
    scored = [s for s in samples if s.time >= start]
    if isempty(scored)
        @warn "No profile sample to pool from $start on."
        return directory
    end

    bands = [band_label(upper, lower) => (upper, lower) for (upper, lower) in validation.depth_bands]
    push!(bands, "all" => (-Inf, Inf))

    stations = sort(unique(s.station for s in scored))
    groups = ["all" => stations; validation.station_groups]
    write_pooled_table(joinpath(directory, "skill_bands.csv"), "group", scored, bands, groups) do sample, members
        sample.station in members
    end
    write_pooled_table(joinpath(directory, "skill_stations.csv"), "station", scored, bands, [s => s for s in stations]) do sample, station
        sample.station == station
    end
    months = [Dates.format(m, "yyyy-mm") => m for m in sort(unique(Dates.firstdayofmonth(Date(s.time)) for s in scored))]
    write_pooled_table(joinpath(directory, "skill_months.csv"), "month", scored, bands, months) do sample, month
        Dates.firstdayofmonth(Date(sample.time)) == month
    end

    write_stratification(directory, samples, start, groups, validation.stratification)
    return directory
end

# `"0-5 m"`, or `">50 m"` for a band with no lower bound.
band_label(upper, lower) = isinf(lower) ? @sprintf(">%g m", upper) : @sprintf("%g-%g m", upper, lower)

"""
    write_samples(path, samples)

`samples.csv`: `station,variable,time,depth,observed,modelled`, one row per pair.
"""
function write_samples(path, samples)
    open(path, "w") do io
        println(io, "station,variable,time,depth,observed,modelled")
        for s in samples
            println(
                io,
                s.station, ",", s.variable, ",", Dates.format(s.time, "yyyy-mm-ddTHH:MM:SS"), ",",
                @sprintf("%.2f,%.4f,%.4f", s.depth, s.observed, s.modelled),
            )
        end
    end
    return path
end

"""
    write_pooled_table(selected, path, column, samples, bands, selections)

One row per selection, variable and depth band: the count, the observed and modelled means, the
bias and the RMSE of the samples `selected(sample, value)` keeps for that selection's `value`.
`selections` are `label => value` pairs, in row order; `column` names the label column. A
combination with fewer than two samples has no row.
"""
function write_pooled_table(selected, path, column, samples, bands, selections)
    variables = sort(unique(s.variable for s in samples))
    open(path, "w") do io
        println(io, column, ",variable,band,count,observed_mean,modelled_mean,bias,rmse")
        for (label, value) in selections
            chosen = [s for s in samples if selected(s, value)]
            for variable in variables, (band, (upper, lower)) in bands
                inside = [s for s in chosen if s.variable == variable && upper <= s.depth < lower]
                metrics = skill_metrics([s.observed for s in inside], [s.modelled for s in inside])
                metrics.count >= 2 || continue
                println(io, label, ",", variable, ",", band, ",", pooled_columns(metrics))
            end
        end
    end
    return path
end

pooled_columns(metrics) = @sprintf(
    "%d,%.3f,%.3f,%.3f,%.3f",
    metrics.count, metrics.observed_mean, metrics.modelled_mean, metrics.bias, metrics.rmse,
)

"""
    write_stratification(directory, samples, start, groups, stratification)

The difference `stratification = (variable, upper, lower)` names — the variable at `lower` minus at
`upper` — for every cast that reaches both depths, observed and modelled, in `stratification.csv`;
and its bias and RMSE per station group from `start` on, in `stratification_skill.csv`.
"""
function write_stratification(directory, samples, start, groups, stratification)
    variable, upper, lower = stratification
    casts = Dict{Tuple{String,DateTime},Vector{SamplePair}}()
    for s in samples
        s.variable == variable && push!(get!(casts, (s.station, s.time), SamplePair[]), s)
    end

    differences = @NamedTuple{station::String, time::DateTime, observed::Float64, modelled::Float64}[]
    for ((station, time), cast) in casts
        sort!(cast; by = s -> s.depth)
        depths = [s.depth for s in cast]
        observed = interpolate_to_depths(depths, [s.observed for s in cast], [upper, lower])
        modelled = interpolate_to_depths(depths, [s.modelled for s in cast], [upper, lower])
        all(isfinite, observed) && all(isfinite, modelled) || continue
        push!(differences, (station, time, observed[2] - observed[1], modelled[2] - modelled[1]))
    end
    sort!(differences; by = d -> (d.station, d.time))

    open(joinpath(directory, "stratification.csv"), "w") do io
        println(io, "station,time,observed,modelled")
        for d in differences
            println(io, d.station, ",", Dates.format(d.time, "yyyy-mm-ddTHH:MM:SS"), @sprintf(",%.3f,%.3f", d.observed, d.modelled))
        end
    end

    open(joinpath(directory, "stratification_skill.csv"), "w") do io
        println(io, "group,count,observed_mean,modelled_mean,bias,rmse")
        for (group, members) in groups
            inside = [d for d in differences if d.time >= start && d.station in members]
            metrics = skill_metrics([d.observed for d in inside], [d.modelled for d in inside])
            metrics.count >= 2 && println(io, group, ",", pooled_columns(metrics))
        end
    end

    return directory
end
