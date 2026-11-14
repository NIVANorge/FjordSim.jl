# Kartverket's tide gauges: observed water level, the one source that is public and needs no key.

# Kartverket's own tide API. `api.sehavniva.no`, which both METreports and much of the literature
# cite, no longer resolves; this is the endpoint serving the same data today.
const KARTVERKET_TIDE_URL = "https://vannstand.kartverket.no/tideapi.php"

"""
    KartverketSeaLevel(; data_root, output_directory = "observations/kartverket",
                       stations, start_date, stop_date, reference = "msl", interval_minutes = 10)

Observed water level from the Norwegian Mapping Authority's permanent tide gauges.

The one source in this module that is public and needs no credential, serving quality-flagged
data at up to 10-minute sampling.

`reference = "msl"` puts the record on mean sea level, which is what the model's free surface `η` is
about. Kartverket's own default is chart datum, which lies well below it — roughly 70 cm in the
Oslofjord — and a comparison made against it would report that offset as a model bias.

Downloads are cached one month per file under `observations_directory(config)`, so an interrupted
fetch resumes and re-running the analysis costs nothing — the same policy the NVE river client uses.

# Fields
- `data_root`, `output_directory`: where the cache goes.
- `stations`: the gauges to fetch. Positions come from Kartverket's own station register.
- `start_date`, `stop_date`: the period to fetch, as `DateTime`s; whole months are cached.
- `reference`: vertical datum code, `"msl"` or `"cd"`.
- `interval_minutes`: sampling interval to request.
"""
Base.@kwdef struct KartverketSeaLevel <: AbstractObservationConfig
    data_root::String
    output_directory::String = joinpath("observations", "kartverket")
    stations::Vector{Station}
    start_date::DateTime
    stop_date::DateTime
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
    fetch_observations(target_grid, config::KartverketSeaLevel)

Fetch every station-month the config's period touches that is not already cached. The gauges are
named, so `target_grid` is not needed.
"""
function fetch_observations(target_grid, config::KartverketSeaLevel)
    directory = observations_directory(config)
    mkpath(directory)

    start, stop = config.start_date, config.stop_date
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
