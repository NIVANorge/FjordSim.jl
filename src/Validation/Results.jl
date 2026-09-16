# Reading a run's station output back. The write side is `StationWriter`/`FieldStationWriter` in
# `Simulations.jl`; this is its inverse, and the two are tested against each other rather than
# against a recorded file.
#
# A station file's `time` is seconds from its own run's zero and carries no calendar, exactly as a
# snapshot's is, so the calendar comes from the `start_date` global attribute the writer records —
# the same convention `FromResults` uses to warm-start from a previous run.

"""
    ModelSeries

A run's output at one station for one variable, shaped like an `ObservationSeries` so the two can be
scored against each other without either knowing about the other.

# Fields
- `station`, `variable`: as named by the writer.
- `times`: calendar times, resolved through the file's `start_date` attribute.
- `depths`: metres below the surface, positive down, **shallowest first** — the reverse of the
  model's own bottom-to-top index order, so that it matches how an instrument reports a profile and
  how `ObservationSeries` stores one.
- `values`: `(length(times), length(depths))`, `NaN` in dry cells.
- `metadata`: the file's global attributes, which carry where the station asked to be, where it
  ended up, and how far apart those are.
"""
struct ModelSeries
    station::String
    variable::String
    times::Vector{DateTime}
    depths::Vector{Float64}
    values::Matrix{Float64}
    metadata::Dict{String,Any}
end

"""
    station_files(directory, writer_name)

Every station file a writer produced, as `tag => path`, newest run first.

Files are named `<stem>_<run tag>_<station tag>.<ext>` and the run tag is a timestamp, so sorting
the paths in reverse puts the most recent run first — which is what a validation wants by default
when a results directory holds several.
"""
function station_files(directory, writer_name::AbstractString)
    isdir(directory) || return Pair{String,String}[]

    pattern = Regex("^" * writer_name * "_(\\d{8}T\\d{6}(?:_loop\\d+)?)_(.+)\\.(nc|jld2)\$")
    found = Pair{String,String}[]
    for file in sort(readdir(directory); rev = true)
        match_result = match(pattern, file)
        isnothing(match_result) && continue
        push!(found, match_result.captures[2] => joinpath(directory, file))
    end

    return found
end

"""
    read_station_netcdf(path, variable)

One variable of one `StationWriter` file, as a `ModelSeries`.

The spatial dimensions are singletons — the writer wrote `indices = (i, j, :)` — so they are dropped
here. Which dimension is which is read from the variable's own `dimnames` rather than assumed:
Oceananigans names them by field location, so a `Face`-located velocity and a `Center`-located
tracer do not share a vertical dimension name.
"""
function read_station_netcdf(path, variable::AbstractString)
    return NCDataset(path) do ds
        haskey(ds, variable) || throw(
            ArgumentError(
                "Station file $path has no variable $variable. It carries " *
                "$(join(sort(collect(keys(ds))), ", ")).",
            ),
        )

        metadata = Dict{String,Any}(key => value for (key, value) in ds.attrib)
        reference = station_reference_date(metadata, path)

        raw = ds[variable]
        dimensions = collect(NCDatasets.dimnames(raw))
        time_index = findfirst(name -> occursin("time", lowercase(name)), dimensions)
        isnothing(time_index) &&
            throw(ArgumentError("Station variable $variable in $path has no time dimension."))

        seconds = Array{Float64}(ds["time"][:])
        times = [reference + Dates.Millisecond(round(Int, 1000 * second)) for second in seconds]

        # The vertical dimension is the one that is neither time nor a singleton; a surface field
        # has none, and then the series is one column deep.
        sizes = [NCDatasets.dim(ds, name) for name in dimensions]
        vertical_index = findfirst(
            index -> index != time_index && sizes[index] > 1,
            eachindex(dimensions),
        )

        data = Array(raw[ntuple(_ -> :, length(dimensions))...])
        values = station_value_matrix(data, time_index, vertical_index)

        depths = if isnothing(vertical_index)
            [0.0]
        else
            name = dimensions[vertical_index]
            haskey(ds, name) ? -Array{Float64}(ds[name][:]) : collect(0.0:(size(values, 2) - 1))
        end

        # Shallowest first, so a model profile and an instrument profile are indexed the same way.
        order = sortperm(depths)
        return ModelSeries(
            get(metadata, "station_name", station_name_from_path(path)),
            String(variable),
            times,
            depths[order],
            values[:, order],
            metadata,
        )
    end
end

"""
    station_value_matrix(data, time_index, vertical_index)

Reshape a station variable's array to `(time, depth)`, dropping the singleton spatial dimensions and
turning `missing` into `NaN`.
"""
function station_value_matrix(data, time_index, vertical_index)
    times = size(data, time_index)
    levels = isnothing(vertical_index) ? 1 : size(data, vertical_index)
    values = Matrix{Float64}(undef, times, levels)

    for index in CartesianIndices(data)
        t = index[time_index]
        k = isnothing(vertical_index) ? 1 : index[vertical_index]
        value = data[index]
        values[t, k] = ismissing(value) ? NaN : Float64(value)
    end

    return values
end

"""
    station_reference_date(metadata, path)

The calendar instant a station file's model-time zero stands for, from the `start_date` global
attribute the writer records. Raises rather than guessing: a series placed on the wrong calendar
would be scored against the wrong observations and nothing downstream could tell.
"""
function station_reference_date(metadata, path)
    haskey(metadata, RESULTS_START_DATE_ATTRIBUTE) || throw(
        ArgumentError(
            "Station file $path carries no `$(RESULTS_START_DATE_ATTRIBUTE)` attribute, so its " *
            "time axis cannot be placed on a calendar.",
        ),
    )

    return DateTime(string(metadata[RESULTS_START_DATE_ATTRIBUTE]))
end

"""
    station_name_from_path(path)

The station tag embedded in a station filename, used only when a file carries no `station_name`
attribute — which is the JLD2 case, Oceananigans' layout having nowhere to put one.
"""
function station_name_from_path(path)
    stem = first(splitext(basename(path)))
    match_result = match(r"_\d{8}T\d{6}(?:_loop\d+)?_(.+)$", stem)
    return isnothing(match_result) ? stem : match_result.captures[1]
end

"""
    read_station_jld2(path, variable, reference)

One variable of one `FieldStationWriter` file, as a `ModelSeries`.

Oceananigans' JLD2 layout stores `timeseries/<name>/<iteration>` per record and
`timeseries/t/<iteration>` the model time in seconds, and stores **no grid and no attributes** — so
unlike the NetCDF path, `reference` has to be supplied by the caller from the simulation config.
"""
function read_station_jld2(path, variable::AbstractString, reference::DateTime)
    return JLD2.jldopen(path, "r") do file
        group = "timeseries/$variable"
        haskey(file, group) ||
            throw(ArgumentError("Station file $path has no timeseries for $variable."))

        iterations = sort(parse.(Int, keys(file[group])))
        times = DateTime[]
        columns = Vector{Float64}[]

        for iteration in iterations
            seconds = file["timeseries/t/$iteration"]
            push!(times, reference + Dates.Millisecond(round(Int, 1000 * seconds)))
            record = file["$group/$iteration"]
            push!(columns, [ismissing(v) ? NaN : Float64(v) for v in vec(record)])
        end

        levels = isempty(columns) ? 1 : length(first(columns))
        values = Matrix{Float64}(undef, length(columns), levels)
        for (index, column) in enumerate(columns)
            values[index, :] = column
        end

        return ModelSeries(
            station_name_from_path(path),
            String(variable),
            times,
            levels == 1 ? [0.0] : collect(0.0:(levels - 1)),
            values,
            Dict{String,Any}(),
        )
    end
end

"""
    match_series(observation, model; depth = nothing, tolerance = Dates.Hour(1))

Pair an observation series with a model series at a common depth and on a common time axis, as
`(times, observed, modelled)` ready for `skill_metrics`.

Each observation time is matched to the *nearest* model record, and dropped when the nearest is
further away than `tolerance`. Nearest-in-time rather than interpolated-in-time on purpose: a CTD
cast is an instant, and interpolating the model across a tidal cycle to meet it would smooth away
exactly the variability being scored.

`depth` picks the level to compare at; `nothing` takes the shallowest of each. The model's profile
is interpolated onto the requested depth, so the comparison is not silently at whichever level the
vertical grid happens to put there.
"""
function match_series(
    observation::ObservationSeries,
    model::ModelSeries;
    depth = nothing,
    tolerance = Dates.Hour(1),
)
    target = isnothing(depth) ? min(first(observation.depths), first(model.depths)) : Float64(depth)

    observed_column = depth_column(observation.depths, observation.values, target)
    modelled_column = depth_column(model.depths, model.values, target)

    times = DateTime[]
    observed = Float64[]
    modelled = Float64[]
    isempty(model.times) && return (times, observed, modelled)

    limit = Dates.Millisecond(tolerance)
    for (index, time) in enumerate(observation.times)
        nearest = searchsortedfirst(model.times, time)
        best = nothing
        for candidate in (nearest - 1, nearest)
            checkbounds(Bool, model.times, candidate) || continue
            gap = abs(Dates.Millisecond(model.times[candidate] - time))
            if gap <= limit && (isnothing(best) || gap < best[2])
                best = (candidate, gap)
            end
        end
        isnothing(best) && continue

        push!(times, time)
        push!(observed, observed_column[index])
        push!(modelled, modelled_column[best[1]])
    end

    return (times, observed, modelled)
end

"""
    depth_column(depths, values, target)

The column of `values` at `target` depth, interpolated between the two bracketing levels and `NaN`
outside the range. Shared by both series types, so a model profile and an instrument profile are
sampled by the same rule.
"""
function depth_column(depths, values, target)
    length(depths) == 1 && return values[:, 1]

    upper = searchsortedfirst(depths, target)
    upper > length(depths) && return fill(NaN, size(values, 1))
    depths[upper] == target && return values[:, upper]
    upper == 1 && return fill(NaN, size(values, 1))

    lower = upper - 1
    weight = (target - depths[lower]) / (depths[upper] - depths[lower])
    return (1 - weight) .* values[:, lower] .+ weight .* values[:, upper]
end
