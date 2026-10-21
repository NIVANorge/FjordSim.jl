module Rivers

# Rivers: the generic half of the river pipeline. `add_rivers` writes every river as a volume inflow
# into a set of surface columns near its mouth, in a river file of its own on its own daily axis.
# Only the four hooks `river_locations`, `river_series`, `river_discharge` and `download_rivers` are
# dataset-specific — `nve_rivers.jl` is the template to copy for a new river dataset.
#
# The run reads that file back through `river_forcing`, at the bottom of this file: the inflow
# becomes the free-surface source `forcing.η = Q / A`, and each river variable a source in the top
# cell of the same columns. `RiverTracerSource` says why that pair is what makes the budgets close.

export add_rivers,
    river_forcing_path,
    download_rivers,
    river_locations,
    river_series,
    river_discharge,
    river_search_radius,
    river_maximum_inflow_velocity,
    river_forcing,
    river_date_range,
    read_river_columns,
    combine_forcing,
    RiverLocation,
    NVERiversConfig,
    NVERiver

using Oceananigans
using Oceananigans.Units: Time
using Oceananigans.OutputReaders: Cyclical, AbstractInMemoryBackend, FlavorOfFTS, time_indices
using Oceananigans.Operators: Az, volume
using Oceananigans: fill_halo_regions!, interior
using Oceananigans.Grids: x_domain, y_domain
using Dates: Date, DateTime, Day
using Adapt
using NCDatasets

import Oceananigans: on_architecture
import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend

using ..Configs:
    AbstractForcingConfig,
    AbstractRiverConfig,
    FjordConfig
using ..Utils: coverage_window
using ..Grids: domain_grid, simulation_grid
using ..Plotting: Plotting, plot_rivers
using ..Bathymetry: bathymetry_path
using ..Forcing: forcing_path
using ..Regridding:
    native_times_to_seconds,
    FORCING_DEFLATE_LEVEL,
    column_wet_levels,
    define_forcing_dimensions!

"""
    river_forcing_path(config)

Resolve `config.output_file` against `config.data_root`: the river file written by `add_rivers`.
Defined for every `AbstractRiverConfig`, so a new river dataset inherits path resolution; one without
those fields overloads it. An absolute `output_file` is returned unchanged.
"""
river_forcing_path(config::AbstractRiverConfig) = joinpath(config.data_root, config.output_file)
Plotting.plot_path(config::AbstractRiverConfig) = joinpath(config.data_root, config.plot_file)

const RIVER_NEIGHBOUR_OFFSETS = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))

# The four neighbours a river's water spreads to. Edge neighbours only, because a diagonal step
# between two wet columns can cross a land corner the model has no face through.
const RIVER_SPREAD_OFFSETS = ((-1, 0), (1, 0), (0, -1), (0, 1))

# Spacing of a river file's time axis. River records are daily and the river hooks match them by
# calendar date, so a finer axis would only repeat records.
const RIVER_TIME_SPACING = Day(1)

# The radius in cells of the disk `narrow_water` opens the water with: a channel up to twice this
# wide counts as a river bed. Glomma's two branches through Fredrikstad are one to four cells wide on
# the 200 m grids, and the open fjord beside them is far wider.
const RIVER_BED_RADIUS = 2

# The variables of a river file the run reads by name, besides one per river variable.
const RIVER_INFLOW_NAME = "inflow"
const RIVER_OWNER_NAME = "river"

"""
    RiverLocation(id, name, longitude, latitude)

One river outlet as the source dataset describes it, before it is snapped to a grid cell.
"""
struct RiverLocation
    id::Int
    name::String
    longitude::Float64
    latitude::Float64
end

"""
    RiverMouth(location, i, j, distance, columns, area)

A `RiverLocation` snapped to the coastal water cell `(i, j)`, and the surface columns its water
enters through. `distance` is how many cells the outlet had to move to reach a coastal water cell;
`columns` are the `(i, j)` pairs that carry its inflow, the mouth first unless another river had
already claimed it; `area` is their summed horizontal area in m², which the discharge is divided by.
See `grow_river_columns`.
"""
struct RiverMouth
    location::RiverLocation
    i::Int
    j::Int
    distance::Float64
    columns::Vector{Tuple{Int,Int}}
    area::Float64
end

"""
    river_locations(config)

The river outlets of `config`, as a `Vector{RiverLocation}`. Required hook for every
`AbstractRiverConfig`.
"""
function river_locations end

"""
    river_series(config, times)

The properties of the water each river carries, as a `Dict` mapping a FjordSim tracer name (`"T"`,
`"S"`, ...) to a `(river, time)` matrix with one row per `river_locations` entry and one column per
entry of `times`. Required hook for every `AbstractRiverConfig`.

`NaN` is a value, not a gap: a river whose temperature is `NaN` carries the temperature of the water
it enters, so it still adds its volume without setting what that volume is.
"""
function river_series end

"""
    river_discharge(config, times)

The discharge of each river, in m³/s, as a `(river, time)` matrix with one row per
`river_locations` entry and one column per entry of `times`. Required hook for every
`AbstractRiverConfig`: it is the volume a river adds, and so what moves the free surface and drives
the estuarine circulation. Must be finite and non-negative throughout.
"""
function river_discharge end

"""
    download_rivers(target_grid, config)

Fetch the source data of a river dataset. Only needed by a dataset that downloads.

`target_grid` is the bare domain grid, and is there for the same reason `download_forcing` and
`download_boundaries` take one: a dataset that has to *find* its rivers rather than being told
where they are needs the domain bounds, which `x_domain`/`y_domain` give for any grid. A dataset
whose outlets are stated in the config ignores it.
"""
function download_rivers end

"""
    river_search_radius(config)

How far, in grid cells, `add_rivers` looks for a coastal water cell when an outlet does not
land on one. Defaults to `config.search_radius`.
"""
river_search_radius(config::AbstractRiverConfig) = config.search_radius

"""
    river_maximum_inflow_velocity(config)

The fastest a river's water may cross the surface face of the columns it enters through, in m/s.
Defaults to `config.maximum_inflow_velocity`.

The volume a river adds enters through the top face of its columns, at `w = Q / A` for a discharge
`Q` spread over a summed area `A`, so this is a vertical CFL bound in disguise: `w Δt / Δz` of the
top cell. The top-cell tracer source runs at the same rate, `Q / V = w / Δz`. One 193 m column
taking Glomma's 470 m³/s is `w ≈ 1.2 cm/s`, which through a 1 m top cell at a 180 s step is a
Courant number of 2. `add_rivers` therefore spreads each river over as many columns as its *peak*
discharge needs to stay below this speed — see `grow_river_columns`.
"""
river_maximum_inflow_velocity(config::AbstractRiverConfig) = config.maximum_inflow_velocity

"""
    column_areas(target_grid)

The horizontal area of every surface cell of `target_grid` in m², as an `(i, j)` matrix — what a
river's discharge is divided by to give the speed it enters at.
"""
function column_areas(target_grid)
    Nx, Ny, Nz = size(target_grid)
    return [Az(i, j, Nz, target_grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
end

"""
    is_coastal_cell(mask, i, j)

Whether `(i, j)` is a water cell with at least one land cell among its eight neighbours. A
river has to enter the domain at the coastline: an outlet snapped to open water would inject
freshwater in the middle of the fjord.
"""
function is_coastal_cell(mask, i, j)
    checkbounds(Bool, mask, i, j) || return false
    mask[i, j] || return false

    for (di, dj) in RIVER_NEIGHBOUR_OFFSETS
        neighbour_i = i + di
        neighbour_j = j + dj
        checkbounds(Bool, mask, neighbour_i, neighbour_j) || continue
        mask[neighbour_i, neighbour_j] || return true
    end

    return false
end

"""
    nearest_coastal_cell(mask, i, j, radius)

The coastal water cell closest to `(i, j)`, as `(i, j, distance)`, searching rings of growing
radius out to `radius` cells. Returns `nothing` when no ring holds one.

Ties are broken by the iteration order — latitude offset ascending outermost, longitude offset
ascending inside it — which is the order the reference implementation uses, so both pick the
same cell.
"""
function nearest_coastal_cell(mask, i, j, radius)
    is_coastal_cell(mask, i, j) && return (i, j, 0.0)

    for ring = 1:radius
        best = nothing

        for dj = -ring:ring, di = -ring:ring
            distance = sqrt(di^2 + dj^2)
            ring - 1//2 <= distance <= ring + 1//2 || continue
            is_coastal_cell(mask, i + di, j + dj) || continue
            if isnothing(best) || distance < best[3]
                best = (i + di, j + dj, distance)
            end
        end

        isnothing(best) || return best
    end

    return nothing
end

"""
    narrow_water(water, radius)

The water cells of `water` that lie in a channel too narrow for a disk of `radius` cells: everything
a morphological opening with that disk removes. A disk fits wherever the open fjord is, and the
opening — erode to the cells the whole disk fits around, then dilate back by the same disk — gives
the open fjord back, so what is left over is the channels no wider than `2 radius` cells. Outside
the grid counts as water, so an open boundary does not make the water beside it narrow.
"""
function narrow_water(water, radius)
    disk = [
        (di, dj) for di = -radius:radius, dj = -radius:radius if di^2 + dj^2 <= radius * (radius + 1)
    ]
    wet(i, j) = !checkbounds(Bool, water, i, j) || water[i, j]
    core = [
        water[i, j] && all(wet(i + di, j + dj) for (di, dj) in disk)
        for i in axes(water, 1), j in axes(water, 2)
    ]

    open = falses(size(water))
    for index in findall(core), (di, dj) in disk
        cell = (index[1] + di, index[2] + dj)
        checkbounds(Bool, open, cell...) && (open[cell...] = true)
    end

    return water .& .!open
end

"""
    river_bed_head(water, narrow, i, j)

Where a river resolved by the grid starts: the cell of the `narrow_water` channel at `(i, j)`, or
beside it, farthest along the channel from open water — or `(i, j)` itself when there is no such
channel, or when the channel never reaches open water (the whole water body is narrow, so it has no
"inland").

The grid's water often runs up a river beyond where a river dataset puts its outlet: the outlet is
where the river meets the sea in the dataset's own map, while the bathymetry's land mask may carry
the lower river as water. An outlet there would enter at the seaward end of its own bed, leave
the bed a stagnant dead end, and put its fresh water straight into the fjord. Moving the mouth to
the head makes the bed carry the river, which is what it is.

Distance is counted in steps through the channel, from the channel cells that touch open water;
equally distant heads go to the one nearest `(i, j)`.
"""
function river_bed_head(water, narrow, i, j)
    seeds = [(i + di, j + dj) for (di, dj) in ((0, 0), RIVER_NEIGHBOUR_OFFSETS...)
             if checkbounds(Bool, narrow, i + di, j + dj) && narrow[i + di, j + dj]]
    isempty(seeds) && return (i, j)

    distance = Dict{Tuple{Int,Int},Int}()
    queue = Tuple{Int,Int}[]
    channel = Set([first(seeds)])
    pending = [first(seeds)]
    while !isempty(pending)
        cell = pop!(pending)
        for (di, dj) in RIVER_SPREAD_OFFSETS
            neighbour = (cell[1] + di, cell[2] + dj)
            checkbounds(Bool, narrow, neighbour...) || continue
            if narrow[neighbour...]
                neighbour in channel || (push!(channel, neighbour); push!(pending, neighbour))
            elseif water[neighbour...] && !haskey(distance, cell)
                distance[cell] = 0
                push!(queue, cell)
            end
        end
    end
    isempty(queue) && return (i, j)

    while !isempty(queue)
        cell = popfirst!(queue)
        for (di, dj) in RIVER_SPREAD_OFFSETS
            neighbour = (cell[1] + di, cell[2] + dj)
            neighbour in channel && !haskey(distance, neighbour) || continue
            distance[neighbour] = distance[cell] + 1
            push!(queue, neighbour)
        end
    end

    farthest = maximum(values(distance))
    heads = [cell for (cell, steps) in distance if steps == farthest]
    return heads[argmin([(hi - i)^2 + (hj - j)^2 for (hi, hj) in heads])]
end

"""
    grow_river_columns(water, areas, i, j, required_area, claimed)

The surface columns one river's water enters through, as `(columns, area)`: grown from the mouth
`(i, j)` through edge-connected water, always taking next the column nearest the mouth, until their
summed `areas` reach `required_area` m². The first column is taken whatever its area, so a river
always has one.

Connected rather than merely near, because water only reaches a column the model can carry it to:
the nearest cells by distance can be across a peninsula, in the next arm of the fjord. Nearest
first, because the water spreads from the mouth — the columns are where the plume would be.

`claimed` marks columns an earlier river took. The growth passes *through* them, so a river is not
cut off behind its neighbour, but never takes one, so no column carries two rivers and the file
can say which river owns each. A water body too small for `required_area` is used whole, and the
caller warns.
"""
function grow_river_columns(water, areas, i, j, required_area, claimed)
    columns = Tuple{Int,Int}[]
    area = 0.0
    visited = falses(size(water))
    visited[i, j] = true
    frontier = [(i, j)]

    while !isempty(frontier) && (isempty(columns) || area < required_area)
        _, index = findmin(((ci, cj),) -> (ci - i)^2 + (cj - j)^2, frontier)
        column = popat!(frontier, index)

        if !claimed[column...]
            push!(columns, column)
            area += areas[column...]
        end

        for (di, dj) in RIVER_SPREAD_OFFSETS
            neighbour = (column[1] + di, column[2] + dj)
            checkbounds(Bool, water, neighbour...) || continue
            water[neighbour...] || continue
            visited[neighbour...] && continue
            visited[neighbour...] = true
            push!(frontier, neighbour)
        end
    end

    return columns, area
end

"""
    river_mouths(target_grid, locations, radius, peak_discharges, maximum_velocity)

Place each `RiverLocation` on the grid: snap its outlet to a coastal water cell, then give it the
columns its water enters through. `peak_discharges` is one discharge in m³/s per entry of
`locations`, the largest the run will see, and each river gets columns enough that it enters no
faster than `maximum_velocity` m/s even then.

An outlet is located by independent nearest-node lookups in longitude and latitude, then moved to
the nearest coastal water cell within `radius`. Outlets outside the grid, and outlets with no such
cell in reach, are dropped with a warning — a river placed on land would be silently lost. An outlet
that lands at the seaward end of a river bed the grid resolves is then moved up to where the bed
starts — see `river_bed_head`.

Rivers are given their columns largest first, so where two mouths are close the larger river takes
the water nearest it and the smaller one spreads around it — the one whose inflow speed matters
most is the one that gets the choice. The result keeps the order of `locations`.

Deliberately config-free: it takes the discharges rather than the river config, so the placement
rules can be driven directly from a bare `Vector{RiverLocation}`.
"""
function river_mouths(target_grid, locations, radius, peak_discharges, maximum_velocity)
    water, _ = column_wet_levels(target_grid)
    narrow = narrow_water(water, RIVER_BED_RADIUS)
    areas = column_areas(target_grid)
    longitudes = Array(λnodes(target_grid, Center()))
    latitudes = Array(φnodes(target_grid, Center()))

    snapped = Tuple{Int,Int,Int,Float64}[]
    for (n, location) in enumerate(locations)
        label = "river $(location.id) ($(location.name))"

        if !(first(longitudes) < location.longitude < last(longitudes)) ||
           !(first(latitudes) < location.latitude < last(latitudes))
            @info "Skipping $label: outlet is outside the grid"
            continue
        end

        i = argmin(abs.(longitudes .- location.longitude))
        j = argmin(abs.(latitudes .- location.latitude))

        nearest = nearest_coastal_cell(water, i, j, radius)
        if isnothing(nearest)
            @warn "Skipping $label: no coastal water cell within $radius cells of ($i, $j)"
            continue
        end

        coastal_i, coastal_j, distance = nearest
        head_i, head_j = river_bed_head(water, narrow, coastal_i, coastal_j)
        (head_i, head_j) == (coastal_i, coastal_j) ||
            @info "Moving $label up the river bed the grid resolves, from ($coastal_i, $coastal_j) " *
                  "to its head at ($head_i, $head_j)"

        push!(snapped, (n, head_i, head_j, distance))
    end

    claimed = falses(size(water))
    placed = Dict{Int,RiverMouth}()

    for (n, i, j, distance) in sort(snapped; by = entry -> -peak_discharges[first(entry)])
        location = locations[n]
        required_area = peak_discharges[n] / maximum_velocity
        columns, area = grow_river_columns(water, areas, i, j, required_area, claimed)

        if isempty(columns)
            @warn "Skipping river $(location.id) ($(location.name)): every water column it can " *
                  "reach from ($i, $j) is already carrying another river"
            continue
        end

        area < required_area &&
            @warn "River $(location.id) ($(location.name)) reaches only " *
                  "$(round(area, digits = 0)) of the $(round(required_area, digits = 0)) m² its " *
                  "peak discharge needs, so it enters at up to " *
                  "$(round(peak_discharges[n] / area, sigdigits = 3)) m/s"

        for column in columns
            claimed[column...] = true
        end
        placed[n] = RiverMouth(location, i, j, distance, columns, area)
    end

    return [placed[n] for (n, _, _, _) in snapped if haskey(placed, n)]
end

"""
    add_rivers(config::FjordConfig)
    add_rivers(target_grid, rivers::AbstractRiverConfig; coverage = nothing)

Write the river inflow file at `river_forcing_path(rivers)`. Returns `nothing` when the setup names
no `river_config`.

Each river is a volume of water — `river_discharge` m³/s carrying the `river_series` values — that
enters through the surface of the columns `river_mouths` gives it. The file stores that already on
the grid: the inflow speed of every column, and the values of the water entering it, on a daily axis
spanning the run window. So the run only reads gridded series, and everything about *where* a river
enters is decided here, once, against the same land mask the run uses.

The file stands alone. It is not a copy of the prepared forcing, so this step needs only the
bathymetry and works with no forcing config at all, a change to the forcing needs no re-run of it,
and `FromForcing` initial conditions read a forcing file no river has touched. Its time axis comes from the run window, so it needs
`coverage` — a change of `start_date` or `stop_time` means re-running this step.

The `FjordConfig` method is the setup-level driver: it reads the grid from the processed bathymetry,
so the rivers are placed against the same land mask the model sees, downloads the river data first,
and takes `coverage` from the simulation config.
"""
function add_rivers(config::FjordConfig)
    rivers = config.river_config
    isnothing(rivers) && return nothing
    validate_river_output(rivers, config.forcing_config)

    bathymetry_file = bathymetry_path(config.bathymetry_config)
    isfile(bathymetry_file) || error(
        "Processed bathymetry $bathymetry_file does not exist. " *
        "Run `julia --project -m FjordSim prepare_bathymetry` for this setup first.",
    )

    # The bare domain grid, which is all a download needs: a dataset that discovers its outlets is
    # asking where the domain is, not where the coastline is — the coastline decides where an outlet
    # is *placed*, and that happens below against the immersed grid. Same grid `download_forcing`
    # and `download_boundaries` hand their own datasets, for the same reason.
    download_rivers(domain_grid(config.grid_config, CPU()), rivers)

    # The grid stays on the CPU: building the land masks walks `peripheral_node` cell by cell.
    grid = simulation_grid(config.grid_config, bathymetry_file, CPU())
    result = add_rivers(grid, rivers; coverage = coverage_window(config.simulation_config))

    plot_file = plot_rivers(grid, rivers, river_forcing_path(rivers), result.mouths)

    @info "Placed $(length(result.mouths)) of $(length(river_locations(rivers))) rivers"
    @info "River variables: $(join(result.variables, ", "))"
    @info "Time range: $(first(result.times)) to $(last(result.times)) ($(length(result.times)) steps)"
    @info "River inflow saved to $(result.output_file)"
    @info "Plot saved to $plot_file"

    return result
end

"""
    validate_river_output(rivers, forcing_config)

Refuse a river file that would overwrite the prepared forcing: the two configs name their files
independently, often under one `data_root`.
"""
validate_river_output(rivers::AbstractRiverConfig, ::Nothing) = nothing

function validate_river_output(rivers::AbstractRiverConfig, forcing_config::AbstractForcingConfig)
    abspath(river_forcing_path(rivers)) == abspath(forcing_path(forcing_config)) && error(
        "The river output would overwrite the prepared forcing at $(forcing_path(forcing_config)). " *
        "Give the river config an `output_file` that differs from the forcing config's.",
    )
    return nothing
end

function add_rivers(target_grid, rivers::AbstractRiverConfig; coverage = nothing)
    times = river_forcing_times(coverage)
    locations = river_locations(rivers)
    discharge = river_discharge(rivers, times)
    series = river_series(rivers, times)

    shape = (length(locations), length(times))
    size(discharge) == shape || error(
        "`river_discharge` returned a $(size(discharge)) matrix; it has to be (river, time) = $shape.",
    )
    for (name, values) in series
        size(values) == shape || error(
            "`river_series` returned a $(size(values)) matrix for $name; it has to be " *
            "(river, time) = $shape.",
        )
    end
    all(value -> isfinite(value) && value >= 0, discharge) ||
        error("`river_discharge` has to be finite and non-negative at every river and time.")

    maximum_velocity = river_maximum_inflow_velocity(rivers)
    maximum_velocity > 0 || throw(
        ArgumentError(
            "The river config's `maximum_inflow_velocity` has to be positive, got $maximum_velocity.",
        ),
    )

    peak_discharges = vec(maximum(discharge; dims = 2))
    mouths = river_mouths(
        target_grid, locations, river_search_radius(rivers), peak_discharges, maximum_velocity,
    )
    isempty(mouths) && error("None of the $(length(locations)) river outlets landed on the grid.")
    report_river_mouths(mouths, locations, discharge)

    rows = Dict(location.id => row for (row, location) in enumerate(locations))
    for mouth in mouths, (name, values) in series
        all(isfinite, view(values, rows[mouth.location.id], :)) || error(
            "`river_series` has a non-finite $name for river $(mouth.location.id) " *
            "($(mouth.location.name)). A river needs a value for everything it names at every time: " *
            "one carrying the water it enters is unstable — see `RiverTracerSource`.",
        )
    end

    output_file = river_forcing_path(rivers)
    isdir(dirname(output_file)) || mkpath(dirname(output_file))
    write_river_file(output_file, target_grid, mouths, locations, discharge, series, times)

    @info "Finished adding rivers"

    return (; output_file, mouths, variables = sort(collect(keys(series))), times)
end

"""
    river_forcing_times(coverage)

The time axis a river file is written on: daily, anchored at the run window's first instant and
extended until it passes the last, so `river_date_range` reports a span `validate_time_coverage`
accepts — the reader indexes time `Cyclical()`, which wraps rather than failing outside its data,
and the coverage check is what keeps that wrap unreachable.

Daily because a river dataset's records are daily and the river hooks match them by *calendar
date*, so a finer cadence would ask for the same record several times over. A sub-daily dataset would
want a hook here.
"""
function river_forcing_times(coverage)
    isnothing(coverage) && error(
        "A river file takes its time axis from the run window, but this setup names no " *
        "`simulation_config` to read `start_date` and `stop_time` from. Name one.",
    )

    first_date, last_date = coverage
    dates = collect(first_date:RIVER_TIME_SPACING:last_date)
    # Extend until the axis reaches past the window's end, and to at least two records: the reader
    # keeps two time indices in memory, so a single-record axis reads past the end of the file.
    while last(dates) < last_date || length(dates) < 2
        push!(dates, last(dates) + RIVER_TIME_SPACING)
    end

    return dates
end

"""
    report_river_mouths(mouths, locations, discharge)

Log where each river landed, how many columns carry it, and how fast its water enters at its mean
and its peak discharge. The column count and the peak speed are what `maximum_inflow_velocity`
trades against each other, and neither is derivable from the config alone — both depend on the
discharge series and on how much water lies around the mouth — so a run's log should state them.
"""
function report_river_mouths(mouths, locations, discharge)
    rows = Dict(location.id => row for (row, location) in enumerate(locations))

    for mouth in mouths
        series = view(discharge, rows[mouth.location.id], :)
        mean_discharge = sum(series) / length(series)
        peak_discharge = maximum(series)
        moved = iszero(mouth.distance) ? "" : " (moved $(round(mouth.distance, digits = 2)) cells)"
        @info "River $(mouth.location.id) ($(mouth.location.name)) at cell " *
              "($(mouth.i), $(mouth.j))$moved: $(length(mouth.columns)) column(s), " *
              "$(round(mouth.area / 1e6, sigdigits = 3)) km², " *
              "Q mean $(round(mean_discharge, sigdigits = 3)) / " *
              "peak $(round(peak_discharge, sigdigits = 3)) m³/s, " *
              "inflow mean $(round(mean_discharge / mouth.area, sigdigits = 3)) / " *
              "peak $(round(peak_discharge / mouth.area, sigdigits = 3)) m/s"
    end

    return mouths
end

"""
    write_river_file(output_file, target_grid, mouths, locations, discharge, series, dates)

Write the river inflow file. Its dimensions and coordinates are `define_forcing_dimensions!`'s, so
the reader checks it against the grid exactly as it checks a prepared forcing file, plus a `river`
dimension. It holds:

- `inflow` `(Nx, Ny, time)`, m/s: the speed each column's river water enters at, `Q(t) / area` on a
  river's columns and zero everywhere else. This is what the free surface is forced with.
- one `(Nx, Ny, time)` variable per `river_series` name: the value of the water entering, on a
  river's columns, and `NaN` everywhere else — which the run reads as "nothing enters here".
- `river` `(Nx, Ny)`: which river each column carries, as an index along `river`, `0` for none.
- `discharge` `(river, time)` m³/s, and `river_id`, `river_name`, `river_i`, `river_j` `(river)`:
  what each river is and where its mouth is, for a reader that wants the rivers rather than the
  grid.

The gridded variables are chunked one horizontal slab per record, which is what the run reads.
"""
function write_river_file(output_file, target_grid, mouths, locations, discharge, series, dates)
    rows = [findfirst(location -> location.id == mouth.location.id, locations) for mouth in mouths]
    names = sort(collect(keys(series)))

    @info "Writing river file $output_file with $(length(dates)) time steps: " *
          "$(first(dates)) to $(last(dates))"
    isfile(output_file) && rm(output_file; force = true)

    NCDataset(output_file, "c") do ds
        define_forcing_dimensions!(ds, target_grid, dates)
        defDim(ds, "river", length(mouths))
        Nx, Ny = ds.dim["Nx"], ds.dim["Ny"]

        owner = zeros(Int32, Nx, Ny)
        for (n, mouth) in enumerate(mouths), (i, j) in mouth.columns
            owner[i, j] = n
        end
        defVar(ds, RIVER_OWNER_NAME, owner, ("Nx", "Ny");
            attrib = ["long_name" => "index along `river` of the river a column carries, 0 for none"])
        defVar(ds, "river_id", Int32[mouth.location.id for mouth in mouths], ("river",))
        defVar(ds, "river_name", [mouth.location.name for mouth in mouths], ("river",))
        defVar(ds, "river_i", Int32[mouth.i for mouth in mouths], ("river",))
        defVar(ds, "river_j", Int32[mouth.j for mouth in mouths], ("river",))
        defVar(ds, "discharge", discharge[rows, :], ("river", "time"); attrib = ["units" => "m3 s-1"])

        inflow = defVar(ds, RIVER_INFLOW_NAME, Float32, ("Nx", "Ny", "time");
            chunksizes = [Nx, Ny, 1], deflatelevel = FORCING_DEFLATE_LEVEL,
            attrib = ["units" => "m s-1"])
        values = [
            defVar(ds, name, Float32, ("Nx", "Ny", "time");
                chunksizes = [Nx, Ny, 1], deflatelevel = FORCING_DEFLATE_LEVEL,
                attrib = ["_FillValue" => NaN32]) for name in names
        ]

        for index in eachindex(dates)
            inflow_slab = zeros(Float32, Nx, Ny)
            value_slabs = [fill(NaN32, Nx, Ny) for _ in names]

            for (n, mouth) in enumerate(mouths), (i, j) in mouth.columns
                inflow_slab[i, j] = discharge[rows[n], index] / mouth.area
                for (slab, name) in zip(value_slabs, names)
                    slab[i, j] = series[name][rows[n], index]
                end
            end

            inflow[:, :, index] = inflow_slab
            for (variable, slab) in zip(values, value_slabs)
                variable[:, :, index] = slab
            end
        end
    end

    return output_file
end

# --- Reading ---

"""
    RiverBackend

`FieldTimeSeries` backend for a river file: two time indices in memory, streamed from NetCDF,
exactly like `NetCDFBackend`. A separate type because a river variable is `(Nx, Ny, time)` where a
forcing variable is `(Nx, Ny, Nz, time)`, the same reason `BoundaryBackend` is one.
"""
struct RiverBackend <: AbstractInMemoryBackend{Int}
    start::Int
    length::Int
end

RiverBackend(length) = RiverBackend(1, length)
new_backend(::RiverBackend, start, length) = RiverBackend(start, length)

Base.length(backend::RiverBackend) = backend.length
Base.summary(backend::RiverBackend) = string("RiverBackend(", backend.start, ", ", backend.length, ")")

const RiverFTS = FlavorOfFTS{<:Any,<:Any,<:Any,<:Any,<:RiverBackend}

"""
    load_river_data(; path, variable_name, time_indices_in_memory, reference_date = nothing)

One page of a river variable as an `(Nx, Ny, 1, page)` array, plus the file's time axis in seconds
from `reference_date`. `NaN` is kept: it is what "nothing enters here" reads as.
"""
function load_river_data(; path, variable_name, time_indices_in_memory, reference_date = nothing)
    NCDataset(path) do ds
        variable = ds[variable_name]
        Nx, Ny = size(variable, 1), size(variable, 2)

        data = zeros(Float64, Nx, Ny, 1, length(time_indices_in_memory))
        for (page, index) in enumerate(time_indices_in_memory)
            @views data[:, :, 1, page] .= coalesce.(variable[:, :, index], NaN)
        end

        native_times = ds["time"][:]
        start_time = isnothing(reference_date) ? first(native_times) : reference_date
        times = convert.(Int64, native_times_to_seconds(native_times, start_time))

        return data, times
    end
end

"""
    set!(fts::RiverFTS, path, name)

Page a river series' in-memory window in from the file. Called by
`Oceananigans.OutputReaders.update_field_time_series!` as the run advances, which reaches these
series because `extract_field_time_series` walks the fields of the forcing structs that hold them.
"""
function set!(fts::RiverFTS, path::String = fts.path, name::String = fts.name)
    data, _ = load_river_data(; path, variable_name = name, time_indices_in_memory = time_indices(fts))

    copyto!(interior(fts, :, :, :, :), data)
    fill_halo_regions!(fts)

    return nothing
end

"""
    river_field_series(grid, path, name, reference_date)

The river variable `name` of the file at `path` as a `(Center, Center, Nothing)` `FieldTimeSeries`
on `grid`, two records in memory, indexed `Cyclical()` like every FjordSim reader.
"""
function river_field_series(grid, path, name, reference_date)
    backend = RiverBackend(2)
    time_indices_in_memory = (1, length(backend))
    data, times = load_river_data(; path, variable_name = name, time_indices_in_memory, reference_date)

    fts = FieldTimeSeries{Center,Center,Nothing}(
        grid,
        times;
        backend,
        time_indexing = Cyclical(),
        path,
        name,
    )
    copyto!(interior(fts, :, :, :, :), data)
    fill_halo_regions!(fts)

    return fts
end

"""
    RiverInflow(inflow)

The free-surface forcing of the rivers: `∂η/∂t` gains the inflow speed of every river column, in
m/s, which is `Q / A` — the volume the rivers add, spread over the area they enter through.

It reads nothing but its own series, which is not a choice: the split-explicit free surface calls
`forcing.η` inside every barotropic substep with only `(; η, U, V)` for fields, so a forcing that
wanted a tracer could not have one there.
"""
struct RiverInflow{F}
    inflow::F
end

Adapt.adapt_structure(to, forcing::RiverInflow) = RiverInflow(Adapt.adapt(to, forcing.inflow))
on_architecture(to, forcing::RiverInflow) = RiverInflow(on_architecture(to, forcing.inflow))

@inline (forcing::RiverInflow)(i, j, k, grid, clock, fields) =
    @inbounds forcing.inflow[i, j, 1, Time(clock.time)]

"""
    RiverTracerSource(inflow, value, field_name)

The tracer half of the rivers: in the top cell of a river column, `w (c_r − c) Az / V`, with `w` the
inflow speed, `c_r` the river's value and `c` the cell's own. Zero below the top cell and off the
river columns. `value` is `nothing` for a tracer the river file carries no variable for, which the
rivers then bring in at zero.

Why `c_r − c` and not `c_r`: the rivers' volume enters through the *surface face* — Oceananigans
diagnoses `w` from continuity with no way to put a source in a cell, so the free-surface forcing
shows up as a downward `w` through the top of the column. Tracer advection does not zero its flux
through that face, so the inflow arrives carrying the top cell's own value. This term replaces what
it carried with what the river carries, so the column gains exactly `Q c_r` — the same correction
NumericalEarth applies through its surface boundary condition on a z-star grid.

Every tracer gets one, because a river column left carrying its own value is unstable: the inflow
brings in `c` while upwind-biased WENO reconstructs the outflow faces of a local extremum at less
than `c`, so the cell gains on every step and runs away — 9.5 to 15 °C in five hours in a
stratified test column, and 59 °C beside Numedalslågen, which has no temperature gauge, before
the run crashed. A river with a value of its own pins the cell to it at the same rate, which keeps
it bounded.
"""
struct RiverTracerSource{F,V,N}
    inflow::F
    value::V
    field_name::N
end

Adapt.adapt_structure(to, source::RiverTracerSource) = RiverTracerSource(
    Adapt.adapt(to, source.inflow),
    Adapt.adapt(to, source.value),
    Adapt.adapt(to, source.field_name),
)

on_architecture(to, source::RiverTracerSource) = RiverTracerSource(
    on_architecture(to, source.inflow),
    on_architecture(to, source.value),
    on_architecture(to, source.field_name),
)

@inline river_value(value, i, j, time) = @inbounds value[i, j, 1, time]
@inline river_value(::Nothing, i, j, time) = 0

@inline function (source::RiverTracerSource)(i, j, k, grid, clock, fields)
    time = Time(clock.time)
    inflow = @inbounds source.inflow[i, j, 1, time]
    entering = river_value(source.value, i, j, time)
    value = @inbounds fields[i, j, k, source.field_name]
    rate = inflow * Az(i, j, k, grid, Center(), Center(), Center()) /
           volume(i, j, k, grid, Center(), Center(), Center())
    term = rate * (entering - value)
    return ifelse((k == grid.Nz) & isfinite(term), term, zero(grid))
end

"""
    river_forcing(rivers, grid, filepath, tracers, reference_date)

The forcing terms the river file at `filepath` adds to the model: `η = RiverInflow`, and a
`RiverTracerSource` for every tracer in `tracers`, all sharing one inflow series. A tracer the file
carries a variable for enters at the river's value; one it does not enters at zero, which is logged.
`reference_date` zeroes the time axis, as for every other reader.

A river variable naming no model tracer is ignored — a river dataset may offer more than a given
model carries. `nothing` rivers give an empty `NamedTuple`.
"""
function river_forcing(::AbstractRiverConfig, grid, filepath, tracers, reference_date)
    names = NCDataset(filepath) do ds
        haskey(ds, RIVER_INFLOW_NAME) || error(
            "River file $filepath has no `$RIVER_INFLOW_NAME` variable, so it predates rivers as a " *
            "volume source. Re-run `julia --project -m FjordSim add_rivers` for this setup.",
        )
        size(grid, 1) == ds.dim["Nx"] && size(grid, 2) == ds.dim["Ny"] ||
            throw(DimensionMismatch("river file dimensions not equal to grid dimensions"))
        [String(tracer) for tracer in tracers if haskey(ds, String(tracer))]
    end
    absent = [Symbol(tracer) for tracer in tracers if String(tracer) ∉ names]
    isempty(absent) || @info "The river file names no $(join(absent, ", ")), so the rivers bring " *
                             "them in at zero"

    inflow = river_field_series(grid, filepath, RIVER_INFLOW_NAME, reference_date)
    value(tracer) = String(tracer) in names ?
                    river_field_series(grid, filepath, String(tracer), reference_date) : nothing
    sources = (
        Symbol(tracer) => RiverTracerSource(inflow, value(tracer), Val(Symbol(tracer)))
        for tracer in tracers
    )

    return (; η = RiverInflow(inflow), sources...)
end

river_forcing(::Nothing, grid, filepath, tracers, reference_date) = (;)

"""
    river_date_range(rivers)

First and last date the river file covers, as a `(first, last)` tuple, or `nothing` for a setup with
no rivers. The mirror of `forcing_date_range` and `boundary_date_range`.
"""
river_date_range(rivers::AbstractRiverConfig) = NCDataset(river_forcing_path(rivers)) do ds
    dates = ds["time"][:]
    (first(dates), last(dates))
end

river_date_range(::Nothing) = nothing

"""
    read_river_columns(rivers)

Every column the river file carries a river in, as `(; target_i, target_j)` — the shape
NumericalEarth's `river_mouth_vertical_diffusivity` reads a river routing as — or `nothing` for a
setup with no rivers.
"""
function read_river_columns(rivers::AbstractRiverConfig)
    owner = NCDataset(ds -> ds[RIVER_OWNER_NAME][:, :], river_forcing_path(rivers))
    indices = findall(>(0), owner)
    return (; target_i = [index[1] for index in indices], target_j = [index[2] for index in indices])
end

read_river_columns(::Nothing) = nothing

"""
    combine_forcing(a, b)

Two forcing `NamedTuple`s as one: a field forced by only one keeps its term, and a field forced by
both gets the pair as a tuple, which Oceananigans adds up as `MultipleForcings`.
"""
function combine_forcing(a::NamedTuple, b::NamedTuple)
    shared = (name => (a[name], b[name]) for name in keys(a) if haskey(b, name))
    return merge(a, b, (; shared...))
end

include("nve_rivers.jl")

end # module
