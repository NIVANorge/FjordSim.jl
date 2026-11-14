# Observations someone else downloaded, read from files in the layout `stations.jl` describes.

"""
    CsvObservations(; data_root, output_directory, variables, stations = Station[],
                    units = Dict(), source = "", delimiter = ',', time_format = nothing)

Observations read from CSV files someone supplies, for the programmes that are not public.

Most observation programmes other than tide gauges are held by their owner rather than published —
current-meter moorings, monitoring-programme CTD casts, fixed thermistors. They arrive as spreadsheet
exports in whatever shape their owner keeps them, so guessing at a parser per programme would be
guessing; this reads **one stated schema** instead, and converting an export to it is a few lines of
whatever the exporter speaks.

Point series and profiles are one type rather than two because every one of these records is the
same thing — a value at a station, at a time, at a depth — and a point instrument is a profile with
one level.

# The schema

One file per station and variable, named `<station tag>_<variable>.csv` in
`observations_directory(config)`, where `<station tag>` is `station_tag(station)` — so CTD station
`TØ-1`'s salinity is `TO_1_S.csv`. A header row is required and names the columns; order does not
matter and extra columns are ignored. The stations can be listed in a `stations.csv` beside them
(`read_stations`).

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
- `stations`: the stations to read, at these positions. Empty reads every station in the
  directory's `stations.csv`.
- `units`: per variable, for figure labels. Defaults to empty.
- `source`: what to call the programme in a caption, e.g. `"Statnett/NIVA ADCP"`.
- `delimiter`: `','` by default; `';'` and `'\\t'` are the usual alternatives from a Norwegian Excel.
- `time_format`: a `Dates.DateFormat` when the timestamps are not ISO 8601. `nothing` parses ISO.
"""
Base.@kwdef struct CsvObservations <: AbstractObservationConfig
    data_root::String
    output_directory::String = "observations"
    variables::Vector{String}
    stations::Vector{Station} = Station[]
    units::Dict{String,String} = Dict{String,String}()
    source::String = "CSV"
    delimiter::Char = ','
    time_format::Union{Nothing,Dates.DateFormat} = nothing
end

observation_stations(config::CsvObservations) = isempty(config.stations) ?
    read_stations(joinpath(observations_directory(config), "stations.csv")) : config.stations

"""
    csv_observation_path(config, station, variable)

Where one station-variable file is expected: `<station tag>_<variable>.csv`.
"""
csv_observation_path(config::CsvObservations, station::Station, variable) =
    observation_file(observations_directory(config), station, variable)

"""
    observation_series(config::CsvObservations, station, variable; window)

Read one station-variable file, or `nothing` when this source does not claim that variable or the
file is absent.
"""
function observation_series(config::CsvObservations, station::Station, variable; window)
    name = String(variable)
    name in config.variables || return nothing

    return read_observation_series(
        csv_observation_path(config, station, name),
        station,
        name,
        window;
        units = get(config.units, name, ""),
        source = config.source,
        delimiter = config.delimiter,
        time_format = config.time_format,
    )
end
