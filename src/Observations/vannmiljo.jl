# The Norwegian Environment Agency's Vannmiljø database: the water measurements every monitoring
# programme in Norway registers, whoever made them — NIVA's and IMR's CTD casts among them.

const VANNMILJO_URL = "https://vannmiljoapi.miljodirektoratet.no/api/Public/"
const VANNMILJO_API_KEY_VARIABLE = "VANNMILJO_API_KEY"

# Vannmiljø's parameter ids for the FjordSim variable names, and the units its values are in.
const VANNMILJO_PARAMETERS = Dict("T" => "TEMP", "S" => "SALIN")
const VANNMILJO_UNITS = Dict("T" => "°C", "S" => "psu")

# How long after its first value a value at one location still belongs to the same cast
# (`cast_starts`). A monitoring cruise visits a station once a day.
const VANNMILJO_CAST_SPAN = Dates.Hour(2)

"""
    VannmiljoObservations(; data_root, output_directory = "observations/vannmiljo", start_date,
                          stop_date, variables = ["T", "S"], contractors = String[],
                          sampling_methods = ["Marint vann: Målesonde in situ"],
                          water_categories = ["C"], stations = Station[], match_distance = 500.0,
                          margin = 0.0, source = "Vannmiljø")

Temperature and salinity from Vannmiljø (`vannmiljo.miljodirektoratet.no`), the database every
Norwegian monitoring programme registers its water measurements in. Its API needs a key, free from
the Norwegian Environment Agency, read from the `$VANNMILJO_API_KEY_VARIABLE` environment variable.

`fetch_observations` asks for every registration inside the grid's box over the period, one
calendar year per request, keeps the ones `contractors` and `sampling_methods` name, and writes
them in the layout `stations.jl` describes. A registration's depth is the middle of its sampled
interval: a 1 m CTD bin registered as 70-71 m is 70.5 m, a point sample 5-5 m is 5 m. Values at one
location less than two hours after a cast's first value are that cast, at its first value's time
(`cast_starts`). A value registered twice at the same time and depth is kept once.

The same cast can be registered by more than one programme: each contractor registers its own, so
`contractors` is how a setup keeps one copy. The same station code can name more than one location
(a FerryBox point and a CTD station), and a location's code can be spelt differently from the
programme's own (`"Ep 1"` for Ep1), or be missing; so a station a setup names is matched by
position, not by name.

# Fields
- `data_root`, `output_directory`: where the files go.
- `start_date`, `stop_date`: the period to fetch, as `DateTime`s.
- `variables`: the FjordSim names to fetch, of `"T"` and `"S"`.
- `contractors`: the contractors whose registrations to keep, as Vannmiljø names them (`"NIVA"`,
  `"Havforskningsinstituttet"`). Empty keeps every one.
- `sampling_methods`: the sampling methods to keep, as Vannmiljø names them. The default is CTD
  casts; `"FerryBox sensormålinger (Petersen et al. 2003)"` is the FerryBox. Empty keeps every one.
- `water_categories`: Vannmiljø's water categories to ask for: `"C"` coastal, `"L"` lake, `"R"`
  river.
- `stations`: the stations to write, each from the location nearest it within `match_distance`
  metres, under its own name and position. Empty writes every location found, named by its code
  (`Name_Alternative`), else its name.
- `match_distance`: how far, in metres, a station's location may be from it.
- `margin`: degrees added to every side of the grid's box.
- `source`: what to call the source in a caption.
"""
Base.@kwdef struct VannmiljoObservations <: AbstractObservationConfig
    data_root::String
    output_directory::String = joinpath("observations", "vannmiljo")
    start_date::DateTime
    stop_date::DateTime
    variables::Vector{String} = ["T", "S"]
    contractors::Vector{String} = String[]
    sampling_methods::Vector{String} = ["Marint vann: Målesonde in situ"]
    water_categories::Vector{String} = ["C"]
    stations::Vector{Station} = Station[]
    match_distance::Float64 = 500.0
    margin::Float64 = 0.0
    source::String = "Vannmiljø"
end

observation_stations(config::VannmiljoObservations) = isempty(config.stations) ?
    read_stations(joinpath(observations_directory(config), "stations.csv")) : config.stations

function observation_series(config::VannmiljoObservations, station::Station, variable; window)
    name = String(variable)
    name in config.variables || return nothing

    return read_observation_series(
        observation_file(observations_directory(config), station, name),
        station,
        name,
        window;
        units = get(VANNMILJO_UNITS, name, ""),
        source = config.source,
    )
end

"""
    fetch_observations(target_grid, config::VannmiljoObservations)

Ask Vannmiljø for every registration of `config.variables` inside `target_grid`'s box, widened by
`config.margin`, over the config's period, and write the ones the config keeps
(`write_vannmiljo_observations`). Every file is rewritten: a fetch is the whole period.
"""
function fetch_observations(target_grid, config::VannmiljoObservations)
    key = vannmiljo_api_key()
    longitude, latitude = download_box(target_grid, config.margin)
    box = utm33_box(longitude, latitude)
    parameters = vannmiljo_parameter.(config.variables)

    @info "Asking Vannmiljø for $(join(config.variables, ", ")) over $(config.start_date) .. $(config.stop_date)"
    registrations = Any[]
    for (from, to) in calendar_years(config.start_date, config.stop_date)
        query = Dict(
            "BoundingBox" => box,
            "FromDateSamplingTime" => string(from),
            "ToDateSamplingTime" => string(to),
            "WaterLocationCategoryFilter" => config.water_categories,
            "ParameterIDFilter" => parameters,
        )
        found = vannmiljo_request("GetRegistrations", query, key)
        @info "  $from .. $to: $(length(found)) values"
        append!(registrations, found)
    end

    identifiers = unique(registration["WaterLocationID"] for registration in registrations)
    locations = isempty(identifiers) ? Any[] :
        vannmiljo_request("GetWaterLocations", Dict("WaterLocationIDFilter" => identifiers), key)

    return write_vannmiljo_observations(config, registrations, locations, longitude, latitude)
end

"""
    vannmiljo_api_key()

The key from the `$VANNMILJO_API_KEY_VARIABLE` environment variable. Raises naming it, rather than
letting the request come back as an HTTP 401.
"""
function vannmiljo_api_key()
    key = get(ENV, VANNMILJO_API_KEY_VARIABLE, "")
    isempty(key) && error(
        "No Vannmiljø API key. Set the $VANNMILJO_API_KEY_VARIABLE environment variable; keys are " *
        "free from the Norwegian Environment Agency (vannmiljoapi.miljodirektoratet.no).",
    )
    return key
end

function vannmiljo_parameter(variable)
    haskey(VANNMILJO_PARAMETERS, variable) || throw(ArgumentError(
        "Vannmiljø serves $(join(keys(VANNMILJO_PARAMETERS), " and ")), not $(repr(variable)).",
    ))
    return VANNMILJO_PARAMETERS[variable]
end

"""
    vannmiljo_request(endpoint, query, key)

The `Result` list of one POST to the API. Over the Oslofjord box (10-11.2°E, 59-59.9°N), two years
of T and S were 106 471 values in one 39 s request, so a request per calendar year stays well
inside the timeout.
"""
function vannmiljo_request(endpoint, query, key)
    output = IOBuffer()
    response = Downloads.request(
        VANNMILJO_URL * endpoint;
        method = "POST",
        headers = ["vannmiljoWebAPIKey" => key, "Content-Type" => "application/json"],
        input = IOBuffer(JSON.json(query)),
        output,
        timeout = 900,
    )
    body = String(take!(output))
    response.status == 200 || error("Vannmiljø's $endpoint answered HTTP $(response.status): $(first(body, 300))")
    return JSON.parse(body)["Result"]
end

"""
    calendar_years(start, stop)

`(from, to)` `Date`s covering `start` to `stop`, one pair per calendar year, both ends inclusive.
"""
calendar_years(start, stop) = [
    (max(Date(start), Date(year, 1, 1)), min(Date(stop), Date(year, 12, 31)))
    for year in Dates.year(start):Dates.year(stop)
]

"""
    utm33_box(longitude, latitude)

The bounding box, in UTM zone 33N metres (EPSG:25833, the frame the API takes), of the
`(minimum, maximum)` longitude/latitude box. From points along all four edges, since a parallel
is curved in UTM.
"""
function utm33_box(longitude, latitude)
    steps = range(0, 1; length = 11)
    along(range_) = [range_[1] + step * (range_[2] - range_[1]) for step in steps]
    x = [along(longitude); along(longitude); fill(longitude[1], 11); fill(longitude[2], 11)]
    y = [fill(latitude[1], 11); fill(latitude[2], 11); along(latitude); along(latitude)]

    ArchGDAL.importEPSG(4326; order = :trad) do source
        ArchGDAL.importEPSG(25833; order = :trad) do target
            ArchGDAL.createcoordtrans(source, target) do transform
                ArchGDAL.transform!(x, y, zeros(length(x)), transform)
            end
        end
    end

    return Dict(
        "xmin" => floor(minimum(x)),
        "ymin" => floor(minimum(y)),
        "xmax" => ceil(maximum(x)),
        "ymax" => ceil(maximum(y)),
    )
end

"""
    write_vannmiljo_observations(config, registrations, locations, longitude, latitude)

Write the `registrations` (the API's `GetRegistrations` records) of the locations (its
`GetWaterLocations` records) inside the `longitude`/`latitude` box, as `config` selects them: the
contractors, sampling methods and period it names, and its stations. Only a location with a kept
registration is a candidate for a station, so another programme's location nearby cannot take it.
Returns the directory.
"""
function write_vannmiljo_observations(config::VannmiljoObservations, registrations, locations, longitude, latitude)
    variables = Dict(VANNMILJO_PARAMETERS[v] => v for v in config.variables)
    samples = Dict{Tuple{Int,String},Vector{ObservationRow}}()
    for registration in registrations
        kept_registration(config, registration) || continue
        variable = get(variables, registration["ParameterID"], nothing)
        isnothing(variable) && continue
        time = DateTime(first(registration["SamplingTime"], 19))
        config.start_date <= time <= config.stop_date || continue
        isnothing(registration["UpperDepth"]) && continue
        upper = Float64(registration["UpperDepth"])
        lower = Float64(something(registration["LowerDepth"], upper))
        row = (time = time, depth = (upper + lower) / 2, value = Float64(registration["RegValue"]))
        push!(get!(samples, (registration["WaterLocationID"], variable), ObservationRow[]), row)
    end

    for identifier in unique(first(key) for key in keys(samples))
        keys_here = [key for key in keys(samples) if first(key) == identifier]
        starts = cast_starts([row.time for key in keys_here for row in samples[key]], VANNMILJO_CAST_SPAN)
        for key in keys_here
            samples[key] = [(time = starts[row.time], depth = row.depth, value = row.value) for row in samples[key]]
        end
    end

    sampled = Set(first(key) for key in keys(samples))
    candidates = [
        location for location in locations
        if location["WaterLocationID"] in sampled &&
           longitude[1] <= location["CoordX_dg"] <= longitude[2] &&
           latitude[1] <= location["CoordY_dg"] <= latitude[2]
    ]
    stations, sources = vannmiljo_stations(config, candidates)

    directory = observations_directory(config)
    mkpath(directory)
    written = Station[]
    labels = String[]
    for (station, location) in zip(stations, sources)
        found = false
        for variable in config.variables
            rows = get(samples, (location["WaterLocationID"], variable), ObservationRow[])
            isempty(rows) && continue
            rows = unique(row -> (row.time, row.depth), rows)
            write_observation_file(observation_file(directory, station, variable), rows)
            found = true
        end
        found || continue
        push!(written, station)
        push!(labels, location_label(location))
    end
    write_stations(joinpath(directory, "stations.csv"), written, labels)

    @info "Wrote $(length(written)) Vannmiljø stations to $directory"
    return directory
end

"""
    cast_starts(times, span)

Each of `times` mapped to the start of the cast it belongs to: the first time of a run of times
each less than `span` after that first one. Vannmiljø registers a CTD value at the hour it was
measured in, so one cast crossing an hour has two times — NIVA's at En-1 on 2024-09-26 started at
08:01 and was registered at 09:00 down to 4.5 m and at 10:00 below — and would otherwise be read as
two casts, each of half the column.
"""
function cast_starts(times, span)
    sorted = sort(unique(times))
    starts = Dict{DateTime,DateTime}()
    isempty(sorted) && return starts
    start = first(sorted)
    for time in sorted
        time - start > span && (start = time)
        starts[time] = start
    end
    return starts
end

# Whether a registration is one `config` keeps: its contractor and sampling method, and an exact
# number (not one registered as below a detection limit).
function kept_registration(config::VannmiljoObservations, registration)
    contractor = something(get(registration, "Contractor", nothing), "")
    method = something(get(registration, "SamplingMethodName", nothing), "")
    isempty(config.contractors) || contractor in config.contractors || return false
    isempty(config.sampling_methods) || method in config.sampling_methods || return false
    something(get(registration, "ValueOperator", nothing), "=") == "=" || return false
    return registration["RegValue"] isa Real
end

"""
    vannmiljo_stations(config, locations)

The stations to write and the location each one's data comes from. With `config.stations`, each
station paired with the location nearest it within `config.match_distance`; a station with none is
reported and left out. Without, one station per location, named by `location_name`, and a name two
locations share made unique with the location's id.
"""
function vannmiljo_stations(config::VannmiljoObservations, locations)
    if isempty(config.stations)
        names = [location_name(location["Name"], get(location, "Betegnelse", nothing)) for location in locations]
        names = [count(==(name), names) > 1 ? "$name-$(location["WaterLocationID"])" : name
                 for (name, location) in zip(names, locations)]
        stations = [Station(name, location["CoordX_dg"], location["CoordY_dg"])
                    for (name, location) in zip(names, locations)]
        order = sortperm(names)
        return stations[order], locations[order]
    end

    stations = Station[]
    sources = Any[]
    for station in config.stations
        distances = [haversine_distance(station.longitude, station.latitude,
                                        location["CoordX_dg"], location["CoordY_dg"])
                     for location in locations]
        nearest = isempty(distances) ? nothing : argmin(distances)
        if isnothing(nearest) || distances[nearest] > config.match_distance
            @info "  no Vannmiljø location within $(config.match_distance) m of $(station.name)"
            continue
        end
        push!(stations, station)
        push!(sources, locations[nearest])
    end
    return stations, sources
end

"""
    location_name(name, code)

What to call a location found without a station list: its code, or its name when it has none, or
the two joined when the code is only a number (`"Tisler-1"`).
"""
function location_name(name, code)
    code = strip(something(code, ""))
    isempty(code) && return String(name)
    return all(isdigit, code) ? "$name-$code" : String(code)
end

location_label(location) = string(location["WaterLocationCode"], " ", location["Name"])
