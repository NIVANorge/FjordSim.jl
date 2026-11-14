# Stations, and the files every source writes and reads: one directory per source, holding
# `stations.csv` (where the stations are) and one `<station tag>_<variable>.csv` per station and
# variable, with a `time,depth,value` row per sample. A source that downloads writes this layout,
# and one reader serves them all, as it serves files downloaded some other way (`CsvObservations`).

"""
    Station(; name, longitude, latitude)

One observation site, named by the position the observation was actually taken at rather than by
a grid index. A run's station writers write a time series at it; a source supplies observations at it.

Snapping to the grid is deliberately left to attach time, where the grid exists, so that the offset
between the two can be logged: a comparison that does not record how far a station moved cannot be
read honestly.

# Fields
- `name`: what the station is called in the source the observations come from, e.g. `"OF-1"` or
  `"Km1"`. Used in the output filename and the `output_writers` key, via `station_tag`.
- `longitude`, `latitude`: degrees east and north.
"""
struct Station
    name::String
    longitude::Float64
    latitude::Float64
end

Station(; name, longitude, latitude) = Station(String(name), Float64(longitude), Float64(latitude))

# Norwegian letters mapped to their conventional ASCII transliterations before anything else is
# stripped, so `TØ-1` and `Ø-1` stay distinguishable (`TO_1`, `O_1`) where a blanket
# non-alphanumeric substitution would collapse them towards each other.
const STATION_TAG_SUBSTITUTIONS = ('Æ' => "AE", 'æ' => "ae", 'Ø' => "O", 'ø' => "o", 'Å' => "A", 'å' => "a")

"""
    station_tag(station)

`station.name` reduced to a filename- and `Symbol`-safe ASCII token.
"""
function station_tag(station::Station)
    name = station.name
    for (character, replacement) in STATION_TAG_SUBSTITUTIONS
        name = replace(name, character => replacement)
    end

    tag = replace(name, r"[^A-Za-z0-9]+" => "_")
    tag = strip(tag, '_')
    isempty(tag) && throw(ArgumentError("Station name $(repr(station.name)) reduces to an empty tag."))

    return tag
end

"""
    haversine_distance(longitude1, latitude1, longitude2, latitude2)

Great-circle distance in metres between two points in degrees: how far a station moved when it was
snapped to the grid, and how far a source's location is from a station.
"""
function haversine_distance(longitude1, latitude1, longitude2, latitude2)
    radius = 6371008.8
    φ1, φ2 = deg2rad(latitude1), deg2rad(latitude2)
    Δφ = φ2 - φ1
    Δλ = deg2rad(longitude2 - longitude1)
    a = sin(Δφ / 2)^2 + cos(φ1) * cos(φ2) * sin(Δλ / 2)^2
    return 2 * radius * asin(min(one(a), sqrt(a)))
end

"""
    observation_file(directory, station, variable)

Where one station's record of one variable is: `<station tag>_<variable>.csv` in `directory`.
"""
observation_file(directory, station::Station, variable) =
    joinpath(directory, "$(station_tag(station))_$(variable).csv")

"""
    split_csv_line(line, delimiter = ',')

The fields of one CSV line, a field in double quotes allowed to hold the delimiter (`""` in it is a
quote). The quotes are removed.
"""
function split_csv_line(line, delimiter = ',')
    fields = String[]
    field = IOBuffer()
    quoted = false
    characters = collect(line)
    index = 1
    while index <= length(characters)
        character = characters[index]
        if character == '"' && quoted && index < length(characters) && characters[index+1] == '"'
            write(field, '"')
            index += 1
        elseif character == '"'
            quoted = !quoted
        elseif character == delimiter && !quoted
            push!(fields, String(take!(field)))
        else
            write(field, character)
        end
        index += 1
    end
    push!(fields, String(take!(field)))
    return fields
end

# A field written so `split_csv_line` reads it back: quoted when it holds a comma or a quote.
csv_field(text) = occursin(r"[\",]", text) ? string('"', replace(text, '"' => "\"\""), '"') : text

"""
    read_stations(path)

The `Station`s in a `stations.csv`: the columns `name` (or `station_code`), `longitude` and
`latitude`, found by the header; other columns are ignored. Empty when the file is absent.
"""
function read_stations(path)
    isfile(path) || return Station[]
    lines = filter(!isempty ∘ strip, readlines(path))
    isempty(lines) && return Station[]

    header = lowercase.(strip.(split_csv_line(first(lines))))
    name = something(findfirst(==("name"), header), findfirst(==("station_code"), header), 0)
    longitude = something(findfirst(==("longitude"), header), 0)
    latitude = something(findfirst(==("latitude"), header), 0)
    min(name, longitude, latitude) == 0 &&
        throw(ArgumentError("$path needs `name`, `longitude` and `latitude` columns; header is $header."))

    stations = Station[]
    for line in Iterators.drop(lines, 1)
        fields = split_csv_line(line)
        push!(stations, Station(
            strip(fields[name]),
            parse(Float64, fields[longitude]),
            parse(Float64, fields[latitude]),
        ))
    end
    return stations
end

"""
    write_stations(path, stations, locations)

`stations.csv`: `name,longitude,latitude,location`, one row per station, `location` saying what
in the source the station's data comes from.
"""
function write_stations(path, stations, locations)
    open(path, "w") do io
        println(io, "name,longitude,latitude,location")
        for (station, location) in zip(stations, locations)
            println(
                io,
                csv_field(station.name), ",",
                @sprintf("%.6f,%.6f", station.longitude, station.latitude), ",",
                csv_field(location),
            )
        end
    end
    return path
end

# One sample: what a source writes and the reader reads back, one per row of a station's file.
const ObservationRow = @NamedTuple{time::DateTime, depth::Float64, value::Float64}

"""
    write_observation_file(path, rows)

One station's record of one variable: `time,depth,value`, one row per sample, by time then depth.
"""
function write_observation_file(path, rows)
    open(path, "w") do io
        println(io, "time,depth,value")
        for row in sort(rows; by = row -> (row.time, row.depth))
            println(io, Dates.format(row.time, "yyyy-mm-ddTHH:MM:SS"), ",", row.depth, ",", row.value)
        end
    end
    return path
end

"""
    read_observation_series(path, station, variable, window; units, source, delimiter = ',',
                            time_format = nothing)

One station's file read as an `ObservationSeries` over `window`, or `nothing` when the file is
absent or holds nothing in `window`. An absent file is the normal case — a programme visits some
stations and not others — and is not warned about.
"""
function read_observation_series(
    path,
    station::Station,
    variable,
    window;
    units,
    source,
    delimiter = ',',
    time_format = nothing,
)
    isfile(path) || return nothing

    rows = read_observation_csv(path, delimiter, time_format)
    start, stop = window
    rows = [row for row in rows if start <= row.time <= stop]
    isempty(rows) && return nothing

    times, depths, values = gather_casts(rows)

    return ObservationSeries(station.name, String(variable), times, depths, values, units, source)
end

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
