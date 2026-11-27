# The Oslofjord setup's own config types and hooks, included by `oslofjorden.jl`: the boundary
# tide and salinity corrections and the cast-corrected initial state. Why each exists: README.md.

using FjordSim: BoundaryData, Grids, InitialState, Plotting
using FjordSim.Configs: AbstractBoundaryDataConfig, AbstractInitialStateConfig,
    AbstractObservationConfig
using FjordSim.Observations: fetch_observations, observation_series, observation_stations
using FjordSim.Validation: TIDAL_CONSTITUENTS
using Dates: Date, DateTime, Dates, Day, Hour
using LinearAlgebra: qr
using NCDatasets: NCDataset
using Printf: @sprintf
using Statistics: mean, median

# --- The boundary tide ------------------------------------------------------------------------
#
# Norkyst-v3's boundary tide is too strong and nearly uniform along the edge, so one gain per
# constituent corrects the whole edge (README.md).

"""
    TidalGain(constituent, amplitude_factor, phase_shift)

The correction to one constituent: its amplitude is multiplied by `amplitude_factor` and its phase
delayed by `phase_shift` degrees (negative is earlier). `constituent` is a name from
`TIDAL_CONSTITUENTS`.
"""
struct TidalGain
    constituent::String
    amplitude_factor::Float64
    phase_shift::Float64
end

"""
    TideCorrectedBoundaries(source, gains)

`source`'s open-boundary data with its tide corrected per constituent by `gains`.

`prepare_boundary_data` writes the file exactly as `source` would, then corrects it in place
(`correct_boundary_tide!`). Every other hook is `source`'s: where the file is, which edges are open,
the download, and how the run reads the file.
"""
struct TideCorrectedBoundaries{C <: AbstractBoundaryDataConfig} <: AbstractBoundaryDataConfig
    source::C
    gains::Vector{TidalGain}
end

BoundaryData.boundary_data_path(config::TideCorrectedBoundaries) =
    BoundaryData.boundary_data_path(config.source)
BoundaryData.boundary_data_directory(config::TideCorrectedBoundaries) =
    BoundaryData.boundary_data_directory(config.source)
Plotting.plot_path(config::TideCorrectedBoundaries) = Plotting.plot_path(config.source)
Grids.open_edges(config::TideCorrectedBoundaries) = Grids.open_edges(config.source)
BoundaryData.boundary_variable_names(config::TideCorrectedBoundaries) =
    BoundaryData.boundary_variable_names(config.source)
BoundaryData.fetch_boundary_data(target_grid, config::TideCorrectedBoundaries) =
    BoundaryData.fetch_boundary_data(target_grid, config.source)

function BoundaryData.regrid_boundary_data(target_grid, config::TideCorrectedBoundaries)
    isempty(config.gains) && error(
        "No boundary tide gains are stated. After `download_boundary_data`, run " *
        "`julia --project scripts/measure_boundary_tide.jl examples/oslofjorden/oslofjorden.jl " *
        "Viker` and paste its output into `VIKER_TIDE_GAINS`.",
    )
    result = BoundaryData.regrid_boundary_data(target_grid, config.source)
    correct_boundary_tide!(result.output_file, Grids.open_edges(config), config.gains)
    return result
end

"""
    tidal_design(hours)

The least-squares design matrix for a mean plus every constituent in `TIDAL_CONSTITUENTS`, as
`harmonic_fit` builds it: column 1 is the mean, columns `2k` and `2k + 1` the cosine and sine of
constituent `k`.
"""
function tidal_design(hours)
    design = ones(length(hours), 1 + 2 * length(TIDAL_CONSTITUENTS))
    for (k, constituent) in enumerate(TIDAL_CONSTITUENTS)
        ω = 2π / constituent.period_hours
        design[:, 2k] = cos.(ω .* hours)
        design[:, 2k + 1] = sin.(ω .* hours)
    end
    return design
end

"""
    tidal_increment(factorization, design, series, gains)

What to add to each column of `series` (time × nodes) to apply `gains`: for each corrected
constituent, the corrected harmonic minus the fitted one.

The fit includes every constituent in `TIDAL_CONSTITUENTS`, not only the corrected ones, so the
others are not aliased into them. A constituent `a cos ωt + b sin ωt` has amplitude `hypot(a, b)`
and phase `atan(b, a)`. Scaling the amplitude by `g` and delaying the phase by `Δ` turns `(a, b)`
into `g` times `(a, b)` rotated by `Δ`.
"""
function tidal_increment(factorization, design, series, gains)
    coefficients = factorization \ series
    change = zero(coefficients)
    for gain in gains
        k = findfirst(c -> c.name == gain.constituent, TIDAL_CONSTITUENTS)
        isnothing(k) && throw(ArgumentError("$(gain.constituent) is not in TIDAL_CONSTITUENTS"))
        a = coefficients[2k, :]
        b = coefficients[2k + 1, :]
        Δ = deg2rad(gain.phase_shift)
        change[2k, :] = gain.amplitude_factor .* (a .* cos(Δ) .- b .* sin(Δ)) .- a
        change[2k + 1, :] = gain.amplitude_factor .* (a .* sin(Δ) .+ b .* cos(Δ)) .- b
    end
    return design * change
end

"""
    correct_boundary_tide!(path, edges, gains)

Correct the tide in the prepared boundary file at `path` on every edge in `edges`, in place.

`eta`, `ubar` and `vbar` are fitted and corrected at every wet boundary node over the file's whole
time axis. What the fit does not take out is left as it was: the mean, the surges, the subtidal
exchange. The velocities get the same gains as the sea level because the Flather condition combines
the exterior transport with the exterior elevation, and correcting only one would change the tide
the boundary lets in by the wrong amount. The 3D `u` and `v` get their depth mean's correction at
every wet level, so the barotropic tide in them stays the one in `ubar` and `vbar`. Dry nodes, all
missing in the file, are skipped. A file this has already corrected is refused rather than corrected
twice.
"""
function correct_boundary_tide!(path, edges, gains)
    NCDataset(path, "a") do ds
        haskey(ds.attrib, "tide_correction") && error(
            "$path already has its boundary tide corrected ($(ds.attrib["tide_correction"])). " *
            "Re-run `prepare_boundary_data` to start from the uncorrected data.",
        )
        times = ds["time"][:]
        design = tidal_design([Dates.value(t - first(times)) / 3_600_000 for t in times])
        factorization = qr(design)

        for edge in edges, (barotropic, full) in (("eta", nothing), ("ubar", "u"), ("vbar", "v"))
            name = "$(edge)_$(barotropic)"
            haskey(ds, name) || continue
            values = ds[name][:, :]
            wet = [i for i in axes(values, 1) if !any(ismissing, view(values, i, :))]
            series = Float64.(permutedims(values[wet, :]))
            increment = tidal_increment(factorization, design, series, gains)
            for (j, i) in enumerate(wet)
                values[i, :] .+= Float32.(view(increment, :, j))
            end
            ds[name][:, :] = values
            @info "Corrected the tide in $name at $(length(wet)) wet nodes"

            isnothing(full) && continue
            field = "$(edge)_$(full)"
            haskey(ds, field) || continue
            values = ds[field][:, :, :]
            for (j, i) in enumerate(wet), k in axes(values, 2)
                level = view(values, i, k, :)
                any(ismissing, level) && continue
                level .+= Float32.(view(increment, :, j))
            end
            ds[field][:, :, :] = values
            @info "Applied the $barotropic correction to every wet level of $field"
        end

        ds.attrib["tide_correction"] = join(
            ("$(g.constituent) x$(g.amplitude_factor) $(g.phase_shift)°" for g in gains), ", ")
    end
    return path
end

# --- The boundary salinity --------------------------------------------------------------------
#
# Norkyst-v3 800 m lacks the brackish surface layer at the mouth, so the edge's salinity is scaled
# to NIVA's FerryBox where the ships cross it (README.md).

"""
    FerryBoxCorrectedBoundaries(source; directory, band, smoothing_days, taper_depth)

`source`'s open-boundary data with the southern edge's salinity corrected to the FerryBox files in
`directory`.

`prepare_boundary_data` writes the file exactly as `source` would, then corrects it in place
(`correct_boundary_salinity!`). Every other hook is `source`'s: where the file is, which edges are
open, the download, and how the run reads the file.
"""
struct FerryBoxCorrectedBoundaries{C <: AbstractBoundaryDataConfig} <: AbstractBoundaryDataConfig
    source::C
    directory::String
    band::Float64
    smoothing_days::Int
    taper_depth::Float64
end

FerryBoxCorrectedBoundaries(source; directory, band, smoothing_days, taper_depth) =
    FerryBoxCorrectedBoundaries(source, String(directory), Float64(band), Int(smoothing_days),
                                Float64(taper_depth))

BoundaryData.boundary_data_path(config::FerryBoxCorrectedBoundaries) =
    BoundaryData.boundary_data_path(config.source)
BoundaryData.boundary_data_directory(config::FerryBoxCorrectedBoundaries) =
    BoundaryData.boundary_data_directory(config.source)
Plotting.plot_path(config::FerryBoxCorrectedBoundaries) = Plotting.plot_path(config.source)
Grids.open_edges(config::FerryBoxCorrectedBoundaries) = Grids.open_edges(config.source)
BoundaryData.boundary_variable_names(config::FerryBoxCorrectedBoundaries) =
    BoundaryData.boundary_variable_names(config.source)
BoundaryData.fetch_boundary_data(target_grid, config::FerryBoxCorrectedBoundaries) =
    BoundaryData.fetch_boundary_data(target_grid, config.source)

function BoundaryData.regrid_boundary_data(target_grid, config::FerryBoxCorrectedBoundaries)
    records = ferrybox_minutes(config.directory)
    result = BoundaryData.regrid_boundary_data(target_grid, config.source)
    correct_boundary_salinity!(result.output_file, records, config)
    return result
end

"""
    ferrybox_minutes(directory)

Every minute with a salinity in the FerryBox files `<ship>_<yyyymm>.csv` in `directory`, as
`(ship, time, longitude, latitude, salinity)`. Each file's columns are `time` (UTC), `lon`, `lat`,
`temperature` and `salinity`.
"""
function ferrybox_minutes(directory)
    pattern = r"^[A-Za-z0-9]+_\d{6}\.csv$"
    files = isdir(directory) ? filter(f -> occursin(pattern, f), readdir(directory)) : String[]
    isempty(files) && error(
        "No FerryBox files `<ship>_<yyyymm>.csv` in $directory. Put NIVA's FerryBox minutes there " *
        "(see examples/oslofjorden/README.md).",
    )
    records = @NamedTuple{
        ship::String, time::DateTime, longitude::Float64, latitude::Float64, salinity::Float64,
    }[]
    for file in files
        ship = first(split(file, "_"))
        for line in Iterators.drop(eachline(joinpath(directory, file)), 1)
            timestamp, longitude, latitude, _, salinity = split(line, ",")
            isempty(salinity) && continue
            push!(records, (ship = ship, time = DateTime(timestamp), longitude = parse(Float64, longitude),
                            latitude = parse(Float64, latitude), salinity = parse(Float64, salinity)))
        end
    end
    return records
end

"""
    piecewise_linear(x, xs, ys)

The piecewise-linear interpolant through `(xs, ys)` at `x`, held at the end values outside `xs`.
"""
function piecewise_linear(x, xs, ys)
    x <= first(xs) && return first(ys)
    x >= last(xs) && return last(ys)
    j = searchsortedlast(xs, x)
    w = (x - xs[j]) / (xs[j + 1] - xs[j])
    return (1 - w) * ys[j] + w * ys[j + 1]
end

"""
    daily_ratio(salinity, longitudes, edge_latitude, times, records, config)

Each ship's daily ratio of the FerryBox's salinity to the boundary's `salinity` (longitude × time,
at the ships' depth), over its minutes within `config.band` degrees north of the edge, and where it
crosses the edge (the median longitude of those minutes). The daily ratios are running means over
`config.smoothing_days`, with what gaps remain filled by linear interpolation. Outside a ship's
record its series is that of the ship with the longest record plus their mean difference.
"""
function daily_ratio(salinity, longitudes, edge_latitude, times, records, config)
    days = Date(first(times)):Day(1):Date(last(times))
    series = Dict{String, Vector{Float64}}()
    crossing = Dict{String, Float64}()
    for ship in unique(m.ship for m in records)
        boundary_sums, ferry_sums = zeros(length(days)), zeros(length(days))
        counts, crossed = zeros(Int, length(days)), Float64[]
        for m in records
            (m.ship == ship && edge_latitude <= m.latitude < edge_latitude + config.band &&
             first(longitudes) <= m.longitude <= last(longitudes) &&
             first(times) <= m.time <= last(times)) || continue
            n = searchsortedfirst(times, m.time)
            n = n > 1 && times[n] - m.time > m.time - times[n - 1] ? n - 1 : n
            boundary = salinity[argmin(abs.(longitudes .- m.longitude)), n]
            isnan(boundary) && continue
            d = Dates.value(Date(m.time) - first(days)) + 1
            boundary_sums[d] += boundary
            ferry_sums[d] += m.salinity
            counts[d] += 1
            push!(crossed, m.longitude)
        end
        isempty(crossed) && continue
        half = config.smoothing_days ÷ 2
        smoothed = map(eachindex(days)) do d
            window = max(1, d - half):min(length(days), d + half)
            known = filter(e -> counts[e] > 0, window)
            isempty(known) ? NaN : mean(ferry_sums[e] / boundary_sums[e] for e in known)
        end
        known = findall(isfinite, smoothed)
        inside = first(known):last(known)
        series[ship] = [k in inside ? piecewise_linear(k, known, smoothed[known]) : NaN for k in eachindex(days)]
        crossing[ship] = median(crossed)
    end
    isempty(series) && error("No FerryBox minute lies within $(config.band)° of the southern edge.")

    reference = argmax(ship -> count(isfinite, series[ship]), collect(keys(series)))
    series[reference] = [piecewise_linear(k, findall(isfinite, series[reference]),
                                     filter(isfinite, series[reference])) for k in eachindex(days)]
    for (ship, values) in series
        ship == reference && continue
        known = isfinite.(values)
        offset = mean(values[known] .- series[reference][known])
        values[.!known] .= series[reference][.!known] .+ offset
    end
    return days, series, crossing
end

"""
    correct_boundary_salinity!(path, records, config)

Correct the southern edge's salinity in the prepared boundary file at `path` to the FerryBox
`records`, in place.

At each hour the ratio is each ship's `daily_ratio`, linear in time between the days' noons and
along the edge between the ships' crossings, held beyond the outermost ones. It scales the salinity
in full at the surface, tapering linearly to no change at `config.taper_depth`: a ratio, as fresh
water the source lacks, since subtracting the difference took the source's fresh pulses below zero.
Dry nodes, all missing in the file, stay missing. A file this has already corrected is refused
rather than corrected twice.
"""
function correct_boundary_salinity!(path, records, config)
    NCDataset(path, "a") do ds
        haskey(ds.attrib, "salinity_correction") && error(
            "$path already has its boundary salinity corrected ($(ds.attrib["salinity_correction"])). " *
            "Re-run `prepare_boundary_data` to start from the uncorrected data.",
        )
        longitudes = ds["Nx"][:]
        depths = -ds["Nz"][:]
        times = ds["time"][:]
        # The ships' intake is at about 4 m.
        ferry_level = argmin(abs.(depths .- 4))
        salinity = Float64.(coalesce.(ds["south_S"][:, ferry_level, :], NaN))
        days, series, crossing = daily_ratio(salinity, longitudes, ds["Ny"][1], times, records, config)

        ships = sort(collect(keys(series)), by = ship -> crossing[ship])
        crossings = [crossing[ship] for ship in ships]
        first_noon = DateTime(first(days)) + Hour(12)
        weights = max.(1 .- depths ./ config.taper_depth, 0)
        corrected_levels = findall(>(0), weights)
        lowest = Inf
        for chunk in Iterators.partition(eachindex(times), 1000)
            values = ds["south_S"][:, :, chunk]
            for (c, n) in enumerate(chunk)
                day = Dates.value(times[n] - first_noon) / 86_400_000 + 1
                at_ships = [piecewise_linear(day, eachindex(days), series[ship]) for ship in ships]
                for i in axes(values, 1)
                    ratio = piecewise_linear(longitudes[i], crossings, at_ships)
                    for k in corrected_levels
                        ismissing(values[i, k, c]) && continue
                        values[i, k, c] *= 1 - weights[k] * (1 - ratio)
                        lowest = min(lowest, values[i, k, c])
                    end
                end
            end
            ds["south_S"][:, :, chunk] = values
        end

        for ship in ships
            @info @sprintf("FerryBox %s crosses the edge at %.2f°E: its salinity over the boundary's %.3f",
                           ship, crossing[ship], mean(series[ship]))
        end
        @info "Corrected south_S down to $(config.taper_depth) m; its lowest corrected value is $lowest g/kg"
        ds.attrib["salinity_correction"] =
            "FerryBox " * join(ships, ", ") * ", tapering to none at $(config.taper_depth) m"
    end
    return path
end

# --- The basin water --------------------------------------------------------------------------
#
# Norkyst-v3's water in the sill basins is wrong, so `CastCorrectedInitialState` replaces it with the
# NIVA casts nearest the start (README.md).

"""
    CastCorrectedInitialState(source, casts, maximum_cast_offset, deep_only)

`source`'s initial state with T and S in the columns over the sill basins replaced by observed casts.

`prepare_initial_state` writes the file exactly as `source` would, then corrects it in place
(`correct_basin_water!`). `download_initial_state` fetches `source`'s download and the casts. Every
other hook is `source`'s: where the file is, its date. `casts` is the observation source the casts
come from, naming their stations; a station is used when it has a cast within
`maximum_cast_offset` of the start. The stations named in `deep_only` replace only the water below
their sill, not whole columns.
"""
struct CastCorrectedInitialState{C <: AbstractInitialStateConfig, O <: AbstractObservationConfig} <:
       AbstractInitialStateConfig
    source::C
    casts::O
    maximum_cast_offset::Day
    deep_only::Vector{String}
end

InitialState.initial_state_path(config::CastCorrectedInitialState) =
    InitialState.initial_state_path(config.source)
InitialState.initial_state_directory(config::CastCorrectedInitialState) =
    InitialState.initial_state_directory(config.source)
InitialState.initial_state_date(config::CastCorrectedInitialState) =
    InitialState.initial_state_date(config.source)
function InitialState.fetch_initial_state(target_grid, config::CastCorrectedInitialState, date)
    fetch_observations(target_grid, config.casts)
    return InitialState.fetch_initial_state(target_grid, config.source, date)
end

function InitialState.regrid_initial_state(
    target_grid,
    config::CastCorrectedInitialState,
    date;
    edges = nothing,
)
    result = InitialState.regrid_initial_state(target_grid, config.source, date; edges)
    correct_basin_water!(result.output_file, edges, config, date)
    return result
end

"""
    trapped_water(wet, edges)

Label the water below a sill: at each level, the wet cells that cannot reach an open edge through
cells wet at that level. `0` where the water is open, otherwise a number per horizontally connected
group of such cells at a level.
"""
function trapped_water(wet, edges)
    Nx, Ny, Nz = size(wet)
    labels = zeros(Int, Nx, Ny, Nz)
    neighbours = ((1, 0), (-1, 0), (0, 1), (0, -1))
    on_edge(i, j) = (:south in edges && j == 1) || (:north in edges && j == Ny) ||
                    (:west in edges && i == 1) || (:east in edges && i == Nx)
    # Every cell reachable from `queue` through wet cells of level `k` not yet `marked`, marked.
    function spread!(marked, queue, k)
        while !isempty(queue)
            i, j = pop!(queue)
            for (di, dj) in neighbours
                a, b = i + di, j + dj
                (1 <= a <= Nx && 1 <= b <= Ny && wet[a, b, k] && marked[a, b] == 0) || continue
                marked[a, b] = marked[i, j]
                push!(queue, (a, b))
            end
        end
        return marked
    end
    group = 0
    for k in 1:Nz
        level = zeros(Int, Nx, Ny)
        open = [(i, j) for i in 1:Nx, j in 1:Ny if wet[i, j, k] && on_edge(i, j)]
        foreach(c -> level[c...] = -1, open)
        spread!(level, open, k)
        for i in 1:Nx, j in 1:Ny
            (wet[i, j, k] && level[i, j] == 0) || continue
            group += 1
            level[i, j] = group
            spread!(level, [(i, j)], k)
        end
        labels[:, :, k] = max.(level, 0)
    end
    return labels
end

"""
    nearest_cast(casts, station, variable, date, offset)

The cast of `variable` at `station` nearest `date` within `offset`, as `(time, depths, values)`
without its missing samples, or `nothing`.
"""
function nearest_cast(casts, station, variable, date, offset)
    series = observation_series(casts, station, variable; window = (date - offset, date + offset))
    isnothing(series) && return nothing
    n = argmin(abs.(series.times .- date))
    keep = findall(isfinite, series.values[n, :])
    isempty(keep) && return nothing
    return (time = series.times[n], depths = series.depths[keep], values = series.values[n, keep])
end

"The cast's value at `depth`, linear between samples and the nearest sample beyond them."
function cast_value(cast, depth)
    d, v = cast.depths, cast.values
    depth <= first(d) && return first(v)
    depth >= last(d) && return last(v)
    n = searchsortedlast(d, depth)
    w = (depth - d[n]) / (d[n+1] - d[n])
    return (1 - w) * v[n] + w * v[n+1]
end

"""
    correct_basin_water!(path, edges, config, date)

Replace T and S in the prepared initial state at `path`, whole columns, with the cast of the station
nearest each column: in every column that reaches water below a sill (`trapped_water`) whose trapped
group, at its level or a shallower one, holds one of `config.casts`' stations; and, up to
`fringe_cells` cells from those, in every column with no open water at or below the corrected
neighbour's topmost trapped cell: the coastal columns too shallow to reach it, and pockets trapped
apart from every station, such as the Oslo harbour basins. A station's column is the deepest wet one
within a cell of its position. A station in `config.deep_only` replaces only the trapped cells of its
columns and has no fringe. Other columns keep `config.source`'s values. A file this has already
corrected is refused.

Whole columns, because Norkyst's water above the sills is 2-3 g/kg saltier than the casts and lay up
to 6.5 kg/m³ denser than the cast water below it. The fringe, because the coastal columns otherwise
kept Norkyst's surface water, 6.9 g/kg (median) saltier than the corrected columns beside them.
"""
function correct_basin_water!(path, edges, config::CastCorrectedInitialState, date; fringe_cells = 8)
    isempty(something(edges, Symbol[])) && error("Telling trapped water from open needs an open edge.")
    NCDataset(path, "a") do ds
        haskey(ds.attrib, "basin_correction") && error(
            "$path already has its basin water corrected ($(ds.attrib["basin_correction"])). " *
            "Re-run `prepare_initial_state` to start from the uncorrected state.",
        )
        longitude, latitude, depth = ds["Nx"][:], ds["Ny"][:], -ds["Nz"][:]
        T = Float32.(coalesce.(ds["T"][:, :, :, 1], NaN32))
        S = Float32.(coalesce.(ds["S"][:, :, :, 1], NaN32))
        wet = isfinite.(T)
        labels = trapped_water(wet, edges)

        stations = []
        for station in observation_stations(config.casts)
            casts = [nearest_cast(config.casts, station, v, date, config.maximum_cast_offset)
                     for v in ("T", "S")]
            if any(isnothing, casts)
                @warn "No T and S cast at $(station.name) within $(config.maximum_cast_offset) of $date; " *
                      "not used"
                continue
            end
            i0 = argmin(abs.(longitude .- station.longitude))
            j0 = argmin(abs.(latitude .- station.latitude))
            column = argmax(c -> count(@view wet[c..., :]),
                            [(i, j) for i in max(i0 - 1, 1):min(i0 + 1, length(longitude)),
                                        j in max(j0 - 1, 1):min(j0 + 1, length(latitude))])
            push!(stations, (; station, column, T = casts[1], S = casts[2],
                               deep_only = station.name in config.deep_only))
        end

        distance(s, i, j) =
            hypot((longitude[i] - s.station.longitude) * cosd(latitude[j]), latitude[j] - s.station.latitude)
        Nx, Ny, Nz = size(wet)
        # Per column: the station whose cast it takes (0: none) and the depth of its topmost trapped cell.
        assigned = zeros(Int, Nx, Ny)
        top_depth = fill(Inf, Nx, Ny)
        for j in 1:Ny, i in 1:Nx
            top = findlast(>(0), @view labels[i, j, :])
            isnothing(top) && continue
            # The stations in the top trapped cell's group, at its level or the nearest shallower one
            # that has any.
            candidates = Int[]
            for level in top:Nz
                labels[i, j, level] == 0 && break
                candidates =
                    [n for (n, s) in enumerate(stations) if labels[s.column..., level] == labels[i, j, level]]
                isempty(candidates) || break
            end
            isempty(candidates) && continue
            assigned[i, j] = argmin(m -> distance(stations[m], i, j), candidates)
            top_depth[i, j] = depth[top]
        end
        closed_below(i, j, d) = all(k -> !wet[i, j, k] || depth[k] < d || labels[i, j, k] > 0, 1:Nz)
        front = [(i, j) for j in 1:Ny, i in 1:Nx if assigned[i, j] > 0 && !stations[assigned[i, j]].deep_only]
        for _ in 1:fringe_cells
            next = Tuple{Int,Int}[]
            for (i, j) in front, (a, b) in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1))
                (1 <= a <= Nx && 1 <= b <= Ny && assigned[a, b] == 0 && wet[a, b, Nz]) || continue
                closed_below(a, b, top_depth[i, j]) || continue
                assigned[a, b] = assigned[i, j]
                top_depth[a, b] = top_depth[i, j]
                push!(next, (a, b))
            end
            front = next
        end

        corrected = zeros(Int, length(stations))
        for j in 1:Ny, i in 1:Nx
            n = assigned[i, j]
            n == 0 && continue
            for k in 1:Nz
                wet[i, j, k] || continue
                stations[n].deep_only && labels[i, j, k] == 0 && continue
                T[i, j, k] = cast_value(stations[n].T, depth[k])
                S[i, j, k] = cast_value(stations[n].S, depth[k])
                corrected[n] += 1
            end
        end

        ds["T"][:, :, :, 1] = T
        ds["S"][:, :, :, 1] = S
        used = ["$(s.station.name) $(Dates.Date(s.T.time))" for s in stations]
        ds.attrib["basin_correction"] =
            "T and S in the columns over the sill basins from the casts of " * join(used, ", ")
        for (s, n) in zip(stations, corrected)
            @info "Basin water from the $(Dates.Date(s.T.time)) cast at $(s.station.name): $n cells"
        end
    end
    return path
end
