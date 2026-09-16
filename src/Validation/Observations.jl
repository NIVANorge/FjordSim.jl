# Observation sources. The same adapter shape every other dataset in FjordSim uses: an abstract
# config supertype, one subtype per source, and hooks dispatching on it — so adding a programme is a
# new subtype and its methods, never an edit to the scoring code.
#
# Every source returns the same `ObservationSeries`, which is what lets `Metrics.jl` stay ignorant of
# where a number came from.

"""
    AbstractObservationConfig

Supertype of the observation sources a validation scores against.

# Required of a subtype
- `data_root`: directory the cached downloads live under.
- `output_directory`: name of this source's directory inside it.

# Hooks
| Hook | Required | Default |
|---|---|---|
| `observation_stations(config)` | yes | — |
| `observation_series(config, station, variable; window)` | yes | — |
| `download_observations(config, window)` | no | no-op, for a source read from local files |
| `observations_directory(config)` | no | `joinpath(data_root, output_directory)` |
"""
abstract type AbstractObservationConfig end

observations_directory(config::AbstractObservationConfig) =
    joinpath(config.data_root, config.output_directory)

"""
    download_observations(config, window)

Fetch and cache everything `config` can supply over `window`, a `(start, stop)` of `DateTime`s.

Defaults to a no-op: a source read from files someone put on disk has nothing to fetch, and should
not have to say so.
"""
download_observations(::AbstractObservationConfig, window) = nothing

"""
    observation_stations(config)

The `Station`s this source can supply, as positions rather than names, so the same lon/lat that
placed the model's output places the observation beside it.
"""
function observation_stations end

"""
    observation_series(config, station, variable; window)

One station's record of one variable over `window`, as an `ObservationSeries`, or `nothing` when the
source has nothing for that pair.

Returning `nothing` rather than raising is deliberate. A validation sweep asks every source for every
station and variable it might hold; most combinations are empty, and that is not an error.
"""
function observation_series end

"""
    ObservationSeries

One station's observations of one variable: a time axis, a depth axis, and a value per pair.

`values` is `(length(times), length(depths))`. A surface or point record has `depths == [0.0]` and
one column. A profile record — a CTD programme, an ADCP — has one row per cast or sample time.

**Profiles are interpolated onto one depth axis at read time.** Instruments do not sample the same
depths twice: a CTD cast returns whatever the winch gave, and an ADCP's bins depend on its range.
Carrying that through would make every comparison a special case, and a Hovmøller diagram needs a
rectangular array anyway. The interpolation is linear and does not extrapolate — a depth outside a
given cast's range is `NaN` for that cast, not the nearest value — so a shallow cast cannot invent
deep water.

# Fields
- `station`: the station name, as the source calls it.
- `variable`: the FjordSim name, matching the model field it is scored against (`"T"`, `"S"`, `"u"`,
  `"v"`, `"eta"`).
- `times`: one `DateTime` per row.
- `depths`: metres below the surface, positive down, **ascending**, one per column.
  `depth_column` and `match_series` both bisect it, so an unsorted axis would silently mis-sample.
- `values`: the observations, `NaN` where missing.
- `units`, `source`: provenance, carried so a figure can be labelled without the config.
"""
struct ObservationSeries
    station::String
    variable::String
    times::Vector{DateTime}
    depths::Vector{Float64}
    values::Matrix{Float64}
    units::String
    source::String
end

"""
    interpolate_to_depths(sample_depths, sample_values, depths)

Linear interpolation of one cast onto a common depth axis, with `NaN` outside the cast's own range.
Used by every profile source so they all present the same shape.
"""
function interpolate_to_depths(sample_depths, sample_values, depths)
    order = sortperm(sample_depths)
    d = Float64[sample_depths[i] for i in order if isfinite(sample_values[i])]
    v = Float64[sample_values[i] for i in order if isfinite(sample_values[i])]

    output = fill(NaN, length(depths))
    length(d) >= 2 || return output

    for (index, depth) in enumerate(depths)
        (depth < first(d) || depth > last(d)) && continue
        upper = searchsortedfirst(d, depth)
        if upper == 1 || d[upper] == depth
            output[index] = v[upper]
        else
            lower = upper - 1
            weight = (depth - d[lower]) / (d[upper] - d[lower])
            output[index] = (1 - weight) * v[lower] + weight * v[upper]
        end
    end

    return output
end

# --- Kartverket sea level ---

# Kartverket's own tide API. `api.sehavniva.no`, which both METreports and much of the literature
# cite, no longer resolves; this is the endpoint serving the same data today.
const KARTVERKET_TIDE_URL = "https://vannstand.kartverket.no/tideapi.php"

"""
    KartverketSeaLevel(; data_root, output_directory = "observations/kartverket",
                       stations, reference = "msl", interval_minutes = 10)

Observed water level from the Norwegian Mapping Authority's permanent tide gauges.

The one source in this module that is public and needs no credential, and the one that scores the
part of the model the FjordOs reports were most satisfied with. Verified to serve 10-minute,
quality-flagged data at Viker, Oscarsborg and Oslo across the whole 2014-2015 window.

`reference = "msl"` puts the record on mean sea level, which is what the model's free surface `η` is
about. Kartverket's own default is chart datum, roughly 70 cm lower in this fjord, and a comparison
made against it would report that offset as a model bias.

Downloads are cached one month per file under `observations_directory(config)`, so an interrupted
fetch resumes and re-running the analysis costs nothing — the same policy the NVE river client uses.

# Fields
- `data_root`, `output_directory`: where the cache goes.
- `stations`: the gauges to fetch. Positions come from Kartverket's own station register.
- `reference`: vertical datum code, `"msl"` or `"cd"`.
- `interval_minutes`: sampling interval to request.
"""
Base.@kwdef struct KartverketSeaLevel <: AbstractObservationConfig
    data_root::String
    output_directory::String = joinpath("observations", "kartverket")
    stations::Vector{Station}
    reference::String = "msl"
    interval_minutes::Int = 10
end

observation_stations(config::KartverketSeaLevel) = config.stations

"""
    kartverket_cache_path(config, station, year, month)

Where one station-month of the response is cached.
"""
kartverket_cache_path(config::KartverketSeaLevel, station::Station, year, month) = joinpath(
    observations_directory(config),
    "$(station_tag(station))_$(year)$(lpad(month, 2, '0')).xml",
)

"""
    download_observations(config::KartverketSeaLevel, window)

Fetch every station-month `window` touches that is not already cached.
"""
function download_observations(config::KartverketSeaLevel, window)
    directory = observations_directory(config)
    mkpath(directory)

    start, stop = window
    @info "Downloading Kartverket water level for $(length(config.stations)) stations to $directory"

    for station in config.stations, (year, month) in months_in(start, stop)
        path = kartverket_cache_path(config, station, year, month)
        isfile(path) && continue

        from = DateTime(year, month, 1)
        to = from + Dates.Month(1)
        url = string(
            KARTVERKET_TIDE_URL,
            "?lat=", station.latitude,
            "&lon=", station.longitude,
            "&fromtime=", Dates.format(from, "yyyy-mm-ddTHH:MM"),
            "&totime=", Dates.format(to, "yyyy-mm-ddTHH:MM"),
            "&datatype=obs",
            "&refcode=", config.reference,
            "&lang=en&dst=0&tzone=0",
            "&interval=", config.interval_minutes,
            "&tide_request=locationdata",
        )

        try
            Downloads.download(url, path)
            @info "  $(station.name) $year-$(lpad(month, 2, '0'))"
        catch exception
            exception isa InterruptException && rethrow()
            @warn "Could not fetch Kartverket water level" station = station.name year month exception
            isfile(path) && rm(path; force = true)
        end
    end

    return directory
end

"""
    observation_series(config::KartverketSeaLevel, station, variable; window)

The station's water level over `window`, as metres on the configured datum, under the FjordSim name
`"eta"`. Any other `variable` returns `nothing`.

Kartverket reports centimetres; the model's `η` is metres, so the conversion happens here rather
than in the scoring code, where a unit error would be invisible.
"""
function observation_series(config::KartverketSeaLevel, station::Station, variable; window)
    string(variable) == "eta" || return nothing

    times = DateTime[]
    values = Float64[]
    start, stop = window

    for (year, month) in months_in(start, stop)
        path = kartverket_cache_path(config, station, year, month)
        isfile(path) || continue

        for (time, value) in parse_kartverket(read(path, String))
            start <= time <= stop || continue
            push!(times, time)
            push!(values, value / 100)
        end
    end

    isempty(times) && return nothing
    order = sortperm(times)

    return ObservationSeries(
        station.name,
        "eta",
        times[order],
        [0.0],
        reshape(values[order], :, 1),
        "m",
        "Kartverket ($(uppercase(config.reference)))",
    )
end

"""
    parse_kartverket(xml)

`(time, value)` pairs from a Kartverket `locationdata` response.

By regex rather than an XML parser: the response is a flat list of `<waterlevel .../>` elements with
no nesting, no namespaces and no character data, and adding an XML dependency to read it would be
the larger change. A malformed or error response simply matches nothing and yields an empty list,
which the caller reads as a missing month.
"""
function parse_kartverket(xml::AbstractString)
    records = Tuple{DateTime,Float64}[]
    pattern = r"<waterlevel\s+value=\"([^\"]+)\"\s+time=\"([^\"]+)\""

    for match in eachmatch(pattern, xml)
        value = tryparse(Float64, match.captures[1])
        isnothing(value) && continue
        # `2015-01-01T00:00:00+00:00`; the request pins `tzone=0` and `dst=0`, so the offset is
        # always UTC and the first 19 characters are the timestamp.
        time = tryparse(DateTime, first(match.captures[2], 19))
        isnothing(time) && continue
        push!(records, (time, value))
    end

    return records
end

"""
    months_in(start, stop)

`(year, month)` for every calendar month `[start, stop]` touches.
"""
function months_in(start::DateTime, stop::DateTime)
    months = Tuple{Int,Int}[]
    current = Date(Dates.year(start), Dates.month(start), 1)
    last = Date(Dates.year(stop), Dates.month(stop), 1)

    while current <= last
        push!(months, (Dates.year(current), Dates.month(current)))
        current += Dates.Month(1)
    end

    return months
end

# --- Held records, read from files ---

"""
    CsvObservations(; data_root, output_directory, variables, stations,
                    units = Dict(), source = "", delimiter = ',', time_format = nothing)

Observations read from CSV files someone supplies, for the programmes that are not public.

Every observation source in the FjordOs evaluation other than Kartverket water level is held by its
owner rather than published — the 2014 Statnett ADCP moorings, the NIVA Ytre Oslofjord CTD
programme, the Fagrådet Inner Oslofjord series, the Scanmar mooring and the beach thermistors. They
arrive as spreadsheet exports in whatever shape their owner keeps them, so guessing at a parser per
programme would be guessing; this reads **one stated schema** instead, and converting an export to
it is a few lines of whatever the exporter speaks.

Point series and profiles are one type rather than two because every one of these records is the
same thing — a value at a station, at a time, at a depth — and a point instrument is a profile with
one level.

# The schema

One file per station and variable, named `<station tag>_<variable>.csv` in
`observations_directory(config)`, where `<station tag>` is `station_tag(station)` — so CTD station
`TØ-1`'s salinity is `TO_1_S.csv`. A header row is required and names the columns; order does not
matter and extra columns are ignored.

| Column | Required | Meaning |
|---|---|---|
| `time` | yes | ISO 8601, e.g. `2015-06-14T12:00:00`, or matching `time_format` |
| `value` | yes | the observation, in the units `units[variable]` claims |
| `depth` | no | metres below the surface, positive down; absent means a surface record |

Rows whose `value` does not parse become `NaN` rather than being dropped, so a gap stays a gap and
the record keeps its time axis. Rows whose `time` does not parse are dropped with a warning, since a
sample with no time cannot be placed at all.

Profiles are grouped by timestamp into casts and interpolated onto the union of the depths seen, by
`interpolate_to_depths` — which does not extrapolate, so a shallow cast contributes `NaN` below its
deepest sample rather than its deepest value.

# Fields
- `data_root`, `output_directory`: where the files are.
- `variables`: the FjordSim names this source supplies — `"T"`, `"S"`, `"u"`, `"v"`, `"eta"`.
- `stations`: the positions, which the source does not carry and the caller must state.
- `units`: per variable, for figure labels. Defaults to empty.
- `source`: what to call the programme in a caption, e.g. `"Statnett/NIVA ADCP"`.
- `delimiter`: `','` by default; `';'` and `'\\t'` are the usual alternatives from a Norwegian Excel.
- `time_format`: a `Dates.DateFormat` when the timestamps are not ISO 8601. `nothing` parses ISO.
"""
Base.@kwdef struct CsvObservations <: AbstractObservationConfig
    data_root::String
    output_directory::String = "observations"
    variables::Vector{String}
    stations::Vector{Station}
    units::Dict{String,String} = Dict{String,String}()
    source::String = "CSV"
    delimiter::Char = ','
    time_format::Union{Nothing,Dates.DateFormat} = nothing
end

observation_stations(config::CsvObservations) = config.stations

"""
    csv_observation_path(config, station, variable)

Where one station-variable file is expected: `<station tag>_<variable>.csv`.
"""
csv_observation_path(config::CsvObservations, station::Station, variable) =
    joinpath(observations_directory(config), "$(station_tag(station))_$(variable).csv")

"""
    observation_series(config::CsvObservations, station, variable; window)

Read one station-variable file, or `nothing` when this source does not claim that variable or the
file is absent. An absent file is the normal case — a programme visits some stations and not others
— and is not warned about.
"""
function observation_series(config::CsvObservations, station::Station, variable; window)
    name = String(variable)
    name in config.variables || return nothing

    path = csv_observation_path(config, station, name)
    isfile(path) || return nothing

    rows = read_observation_csv(path, config.delimiter, config.time_format)
    isempty(rows) && return nothing

    start, stop = window
    rows = [row for row in rows if start <= row.time <= stop]
    isempty(rows) && return nothing

    times, depths, values = gather_casts(rows)

    return ObservationSeries(
        station.name,
        name,
        times,
        depths,
        values,
        get(config.units, name, ""),
        config.source,
    )
end

# One parsed row of an observation CSV. Concretely typed so `gather_casts` can key a
# dictionary on it without boxing.
const ObservationRow = @NamedTuple{time::DateTime, depth::Float64, value::Float64}

"""
    read_observation_csv(path, delimiter, time_format)

`(time, depth, value)` rows from one file, by the schema `CsvObservations` documents.

Hand-parsed rather than through a CSV package: the schema is three columns of numbers and a
timestamp, and taking on a dependency to split on a delimiter would be the larger change. Quoted
fields are not supported, which the schema does not need.
"""
function read_observation_csv(path, delimiter::Char, time_format)
    lines = readlines(path)
    isempty(lines) && return ObservationRow[]

    header = [strip(lowercase(field)) for field in split(first(lines), delimiter)]
    time_column = findfirst(==("time"), header)
    value_column = findfirst(==("value"), header)
    depth_column = findfirst(==("depth"), header)

    isnothing(time_column) && throw(ArgumentError("$path has no `time` column; header is $header."))
    isnothing(value_column) && throw(ArgumentError("$path has no `value` column; header is $header."))

    rows = ObservationRow[]
    dropped = 0
    for line in Iterators.drop(lines, 1)
        isempty(strip(line)) && continue
        fields = split(line, delimiter)
        length(fields) >= max(time_column, value_column) || (dropped += 1; continue)

        time = isnothing(time_format) ? tryparse(DateTime, strip(fields[time_column])) :
               tryparse(DateTime, strip(fields[time_column]), time_format)
        if isnothing(time)
            dropped += 1
            continue
        end

        # A value that does not parse is a gap, not a missing row: the sample was taken and the
        # instrument had nothing to say, and dropping it would shorten the record silently.
        value = something(tryparse(Float64, strip(fields[value_column])), NaN)
        depth = if isnothing(depth_column) || length(fields) < depth_column
            0.0
        else
            something(tryparse(Float64, strip(fields[depth_column])), NaN)
        end

        isfinite(depth) || (dropped += 1; continue)
        push!(rows, (; time, depth, value))
    end

    dropped > 0 && @warn "Dropped $dropped unparseable rows" path
    return rows
end

"""
    gather_casts(rows)

Group `(time, depth, value)` rows into `(times, depths, values)`: one row of `values` per distinct
timestamp, on the union of the depths seen, with each cast interpolated onto it.

The union rather than the intersection, so a deep cast is not truncated to the shallowest one in the
file; `interpolate_to_depths` fills the levels a given cast does not reach with `NaN` rather than
extrapolating.
"""
function gather_casts(rows)
    times = sort!(unique(row.time for row in rows))
    depths = sort!(unique(row.depth for row in rows))
    values = fill(NaN, length(times), length(depths))

    by_time = Dict{DateTime,Vector{ObservationRow}}()
    for row in rows
        push!(get!(by_time, row.time, ObservationRow[]), row)
    end

    for (index, time) in enumerate(times)
        cast = by_time[time]
        if length(depths) == 1
            values[index, 1] = first(cast).value
        elseif length(cast) == 1
            # A single sample cannot be interpolated; place it at its own depth and leave the rest.
            level = searchsortedfirst(depths, first(cast).depth)
            values[index, level] = first(cast).value
        else
            values[index, :] =
                interpolate_to_depths([r.depth for r in cast], [r.value for r in cast], depths)
        end
    end

    return times, depths, values
end
