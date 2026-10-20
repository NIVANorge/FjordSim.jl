module Regridding

# The regridding core shared by the forcing and open-boundary pipelines: source records and the time
# axes built from them, the water masks, the projected source grid, and the kernel that interpolates
# a source slab onto the target nodes. Nothing here names a dataset.

export native_times_to_seconds,
    FORCING_DEFLATE_LEVEL,
    interpolation_architecture,
    ProjectedSourceGrid,
    SourceFill,
    SourceRecord,
    ForcingTimeStep,
    PreparedVariable,
    uniform_time_steps,
    pad_time_steps,
    water_mask,
    column_wet_levels,
    projected_target_nodes,
    source_validity,
    read_last_dimension,
    as_source_slab,
    source_fill,
    fill_source!,
    nearest_valid_map,
    source_level_coordinate,
    source_field_grid,
    set_source_field!,
    interpolate_to_target!,
    SourceReader,
    source_slab,
    blended_slab,
    define_forcing_dimensions!

using Oceananigans
using Oceananigans: fill_halo_regions!
using Oceananigans.Grids: active_cell, peripheral_node
using Oceananigans.Fields: interpolate
using Oceananigans.Utils: launch!
using KernelAbstractions: @kernel, @index, @Const
using Dates: DateTime, Millisecond, Second
import Dates
using ArchGDAL
using CUDA  # needed for GPU() to be usable, and for CUDA.functional()
using NCDatasets

using ..Configs: lateral_edges

function native_times_to_seconds(native_times, start_time = native_times[1])
    return [Second(t - start_time).value for t in native_times]
end

const FORCING_DEFLATE_LEVEL = 5

"""
    interpolation_architecture(config)

Where `prepare_forcing` and `prepare_boundaries` run their interpolation kernel, from
`config.architecture`:

- `:auto`: the GPU when one is usable, else the CPU. The default, so one setup runs unchanged on
  a GPU machine and on a laptop.
- `:gpu`: the GPU, erroring when none is usable rather than silently running ~12x slower.
- `:cpu`: the CPU.

Only the kernel moves; `target_grid` and the masks stay on the CPU because building them walks
`peripheral_node` cell by cell.

`Forcing` and `Boundaries` each add the method reading their config's `architecture` field.
`Simulations.simulation_architecture` also reuses the `Val` methods below, which is why the
resolution and the field read are separate methods rather than one function.
"""
function interpolation_architecture end

interpolation_architecture(::Val{:cpu}) = CPU()
interpolation_architecture(::Val{:auto}) = CUDA.functional() ? GPU() : CPU()

function interpolation_architecture(::Val{:gpu})
    CUDA.functional() || error(
        "architecture = :gpu was requested but no usable GPU was found. " *
        "Use :auto to fall back to the CPU, or :cpu to ask for it explicitly.",
    )
    return GPU()
end

interpolation_architecture(::Val{selector}) where {selector} =
    throw(ArgumentError("architecture must be one of (:auto, :cpu, :gpu), got :$selector"))

"""
    ProjectedSourceGrid

Geometry of a source subset living on a regular grid in projected meters: the coordinates,
its depth levels and the projection they are defined in.

Nothing here is dataset-specific, so any source on a regular projected grid can return one of
these from `forcing_source_grid` and inherit `source_field_grid` and
`projected_target_nodes` unchanged.

# Fields
- `x`, `y`: Projected coordinate centers in meters, regularly spaced.
- `depths`: Depth levels in meters, positive down.
- `proj4`: PROJ.4 definition of the projection.
"""
struct ProjectedSourceGrid
    x::Vector{Float64}
    y::Vector{Float64}
    depths::Vector{Float64}
    proj4::String
end

"""
    SourceFill

Where every masked source cell takes its value from, so that a filled source slab can be handed
to `Oceananigans.Fields.interpolate` without a `NaN` anywhere in it.

Built once per variable from the first time step, because the mask is a property of the dataset
rather than of the step. `nearest` is a linear index into the `(NX, NY)` plane of the same level;
levels with no valid cell at all (source levels deeper than the fjord) are marked in
`level_valid` and filled from the deepest level that does have data.

This is not `NumericalEarth`'s `inpaint_mask!`: that re-`NaN`s every masked cell before
propagating horizontally, discarding its own `continue_downwards!` pass, so a fully masked level
has nothing to spread from. With the default `NearestNeighborInpainting(Inf)` it then never
terminates, and with a finite cap `_fill_nans!` turns the level into zeros — a plausible-looking
fabricated temperature or salinity.

# Fields
- `valid`: Whether each source cell holds data.
- `nearest`: `(NX, NY, levels)` linear index of the nearest valid cell in the same level.
- `level_valid`: Whether a level has any valid cell inside the subset.
"""
struct SourceFill
    valid::BitArray{3}
    nearest::Array{Int32,3}
    level_valid::Vector{Bool}
end

"""
    SourceRecord

One time record present in a downloaded source file.

# Fields
- `date`: Time of the record.
- `filepath`: Source file holding it.
- `index`: Time index within that file.
"""
struct SourceRecord
    date::DateTime
    filepath::String
    index::Int
end

"""
    ForcingTimeStep

One output time step and the source records it is built from.

A step present in the download reads a single record (`lower === upper`, `weight == 0`). A day
missing from the download is the linear interpolation of the records bracketing it, following
`NumericalEarth.DataWrangling.fill_gaps!`.

# Fields
- `date`: Time of the output step.
- `lower`, `upper`: Bracketing source records; equal when the step is present in the download.
- `weight`: Weight of `upper`; zero when the step is present in the download.
"""
struct ForcingTimeStep
    date::DateTime
    lower::SourceRecord
    upper::SourceRecord
    weight::Float32
end

"""
    PreparedVariable

Everything needed to write one forcing variable: where it comes from, the land mask and target
nodes of its Oceananigans location, how to fill its source mask, and its relaxation lambdas.

`x` and `y` are the target nodes projected into the source coordinate system, held as 2D arrays
so the interpolation kernel only indexes them.

Parametric in the number of NetCDF dimension names, which is the one thing an interior forcing
variable and an open-boundary one disagree about: four for the forcing file's
`(along, across, Nz, time)`, three or two for the boundary file, whose across-edge dimension is one
cell and whose surface variables have no depth axis. Everything else — the three-dimensional mask
with its one-cell slab, the projected nodes, the source fill — is shared, which is why
`prepare_boundaries` reuses this type rather than defining its own.
"""
struct PreparedVariable{N}
    source_name::String
    name::String
    dimensions::NTuple{N,String}
    mask::Array{Bool,3}
    x::Matrix{Float64}
    y::Matrix{Float64}
    z::Vector{Float64}
    mask_fill::SourceFill
    lambda::Array{Float32,3}
end

# A hand-written constructor because `Base.@kwdef`-style automatic conversion does not happen for a
# *parametric* struct: the generated constructor keeps the declared field types in its signature, so
# the `BitArray` `water_mask` returns would not match `Array{Bool,3}`.
PreparedVariable(source_name, name, dimensions::NTuple{N,String}, mask, x, y, z, mask_fill, lambda) where {N} =
    PreparedVariable{N}(
        String(source_name),
        String(name),
        dimensions,
        Array{Bool,3}(mask),
        Matrix{Float64}(x),
        Matrix{Float64}(y),
        Vector{Float64}(z),
        mask_fill,
        Array{Float32,3}(lambda),
    )

"""
    uniform_time_steps(records, spacing; max_gap, unit)

Complete `records` to a gap-free axis at `spacing`, linearly interpolating any missing step between
the records bracketing it. Downloaded months are occasionally short a record, which would otherwise
leave a hole in the cyclical forcing period.

A gap longer than `max_gap` steps is still interpolated, but warned about separately: MET's
NorKyst archive is missing whole weeks in 2017 and 2018, and papering over those silently would
misrepresent the forcing. Returned as single-record steps if the times are not on a whole-`spacing`
cadence.

`unit` only names the spacing in the warnings. Written generically because the interior forcing is
daily and the open-boundary data is hourly, and the gap-filling argument is identical for both.
"""
function uniform_time_steps(records, spacing::Dates.Period; max_gap, unit)
    single(record) = ForcingTimeStep(record.date, record, record, 0)

    differences = diff([record.date for record in records])
    if !all(difference -> difference % spacing == Millisecond(0), differences)
        @warn "Source times are not on a whole-$unit cadence; using them as downloaded."
        return map(single, records)
    end

    steps = ForcingTimeStep[]
    interpolated = DateTime[]
    step_milliseconds = Dates.value(Millisecond(spacing))

    for n = 1:length(records)-1
        lower = records[n]
        upper = records[n+1]
        push!(steps, single(lower))

        span = Dates.value(upper.date - lower.date) ÷ step_milliseconds
        span > max_gap + 1 && @warn "Source gap of $(span - 1) $unit(s) after $(lower.date) " *
                                    "exceeds max_gap = $max_gap; interpolating anyway"

        for missing_step = 1:span-1
            date = lower.date + missing_step * spacing
            push!(interpolated, date)
            push!(steps, ForcingTimeStep(date, lower, upper, missing_step / span))
        end
    end
    push!(steps, single(last(records)))

    isempty(interpolated) ||
        @warn "Interpolated $(length(interpolated)) missing $unit(s): $(join(interpolated, ", "))"

    return steps
end

"""
    pad_time_steps(steps, coverage)

Extend `steps` to span `coverage`, replicating the nearest step at each end that falls short. A
`nothing` coverage leaves `steps` alone.

Written after `daily_time_steps` so its whole-day cadence check sees only the downloaded times: a
pad step is deliberately off that cadence whenever the run's window is.

The pad carries the *same* source records and blend weight as the step it copies, so the written
field is identical to its neighbour rather than a re-derived approximation of it. That is what makes
the padding a statement about the time axis alone.

The reason this exists is that a run's window rarely lines up with a dataset's records — NorKyst's
first daily record here is at 12:00, so a run starting at 00:00 on the same day has no forcing for
its first half-day. `forcing_from_file` reads with `Cyclical()` time indexing, which does not fail
outside its data but wraps to the far end, so the shortfall would be filled with the following
December rather than reported. Padding puts a real record at both ends instead, which is what lets
`Simulations.validate_time_coverage` stay a hard check.

That last point is also why a pad may reach **at most one record spacing** past the downloaded axis.
Unbounded, this would defeat the very check it exists to serve: a window running a year past the
data would be "covered" by one replicated December, `forcing_date_range` would report the invented
span, `validate_time_coverage` would pass, and the run would interpolate twelve months between two
identical records — quietly worse than the wrap. Overshooting by more than a spacing is a missing
download, and is reported as one.
"""
pad_time_steps(steps, ::Nothing) = steps

function pad_time_steps(steps, coverage)
    first_needed, last_needed = coverage
    padded = collect(steps)

    if first(padded).date > first_needed
        head = first(padded)
        spacing = forcing_record_spacing(padded)
        head.date - first_needed <= spacing || error(
            "The run starts at $first_needed but the downloaded forcing only begins at " *
            "$(head.date), more than one record spacing " *
            "($(Dates.canonicalize(spacing))) later. Download the preceding period, or move the " *
            "simulation config's `start_date` to " *
            "$(head.date - spacing) or later.",
        )
        pushfirst!(padded, ForcingTimeStep(first_needed, head.lower, head.upper, head.weight))
        @info "Padded forcing start: replicated $(head.date) at $first_needed"
    end

    if last(padded).date < last_needed
        tail = last(padded)
        spacing = forcing_record_spacing(padded)
        last_needed - tail.date <= spacing || error(
            "The run ends at $last_needed but the downloaded forcing already stops at " *
            "$(tail.date), more than one record spacing " *
            "($(Dates.canonicalize(spacing))) earlier. Download the following period, or shorten " *
            "`stop_time` so the run ends by " *
            "$(tail.date + spacing).",
        )
        push!(padded, ForcingTimeStep(last_needed, tail.lower, tail.upper, tail.weight))
        @info "Padded forcing end: replicated $(tail.date) at $last_needed"
    end

    return padded
end

"""
    forcing_record_spacing(steps)

The axis's record spacing, taken from its first pair — `daily_time_steps` has already made it
uniform, so one pair speaks for the whole axis. A single-record axis has no spacing to bound a pad
by, and is an error rather than a guess.
"""
function forcing_record_spacing(steps)
    length(steps) >= 2 || error(
        "Cannot pad a forcing axis of $(length(steps)) record(s): its record spacing, which " *
        "bounds how far a pad may reach, is unknown.",
    )
    return steps[2].date - steps[1].date
end

"""
    water_mask(target_grid, LX, LY, edges)

Water flag per cell of `target_grid` at the location `(LX, LY, Center)`, taken from
`Oceananigans.Grids.peripheral_node` so that it is by construction the same predicate the model
uses. That covers the staggering: a velocity face is peripheral when either tracer cell it
separates is inactive, and `PartialCellBottom` decides "inactive" with its own
`minimum_fractional_cell_height` rather than a hand-rolled cell-center test.

`edges` names the open lateral boundaries, whose outermost face rows `open_boundary_water!` restores —
one `Symbol`, a collection of them, or `nothing`. A setup naming no open-boundary data restores none:
with no exterior state there is nothing to read at those rows, so every lateral boundary is a closed
wall and the peripheral mask stands as it is. A domain in the open ocean names all four and gets every
face row back.

Only the component staggered across a given edge is affected, so the four restores are independent and
their order does not matter — `u` is touched by `:west` and `:east`, `v` by `:south` and `:north`, and
a tracer by none of them.
"""
function water_mask(target_grid, LX, LY, edges)
    mask = peripheral_water_mask(target_grid, LX, LY)
    for edge in lateral_edges(edges)
        open_boundary_water!(mask, target_grid, LX, LY, Val(edge))
    end
    return mask
end

function peripheral_water_mask(target_grid, LX, LY)
    Nx, Ny, Nz = size(target_grid)
    nx = LX === Face ? Nx + 1 : Nx
    ny = LY === Face ? Ny + 1 : Ny

    mask = falses(nx, ny, Nz)
    for k = 1:Nz, j = 1:ny, i = 1:nx
        mask[i, j, k] = !peripheral_node(i, j, k, target_grid, LX(), LY(), Center())
    end

    return mask
end

"""
    open_boundary_water!(mask, target_grid, LX, LY, ::Val{edge})

Mark the velocity faces lying on the open lateral boundary `edge` as water wherever the tracer
cell just inside them is active.

`peripheral_node` treats those faces as peripheral because the tracer cell outside the domain is
a halo cell, and halo cells are inactive in `Bounded` directions. For a closed wall that is
right, but `edge` is where the domain is open and the normal velocity carries a
`NormalFlowBoundaryCondition`, so leaving that row masked would drop the data exactly where it is
needed — both for the prepared forcing's initial conditions and for the boundary file. Only the
component staggered across `edge` is affected, which is why each method tests one location.

One method per edge rather than a four-branch `if`: the edge is a `Val` everywhere it decides
anything in this package.
"""
open_boundary_water!(mask, target_grid, LX, LY, ::Val{edge}) where {edge} = mask

function open_boundary_water!(mask, target_grid, LX, LY, ::Val{:south})
    LY === Face || return mask
    for k = 1:size(target_grid, 3), i in axes(mask, 1)
        mask[i, 1, k] = active_cell(i, 1, k, target_grid)
    end
    return mask
end

function open_boundary_water!(mask, target_grid, LX, LY, ::Val{:north})
    LY === Face || return mask
    Ny = size(target_grid, 2)
    for k = 1:size(target_grid, 3), i in axes(mask, 1)
        mask[i, Ny+1, k] = active_cell(i, Ny, k, target_grid)
    end
    return mask
end

function open_boundary_water!(mask, target_grid, LX, LY, ::Val{:west})
    LX === Face || return mask
    for k = 1:size(target_grid, 3), j in axes(mask, 2)
        mask[1, j, k] = active_cell(1, j, k, target_grid)
    end
    return mask
end

function open_boundary_water!(mask, target_grid, LX, LY, ::Val{:east})
    LX === Face || return mask
    Nx = size(target_grid, 1)
    for k = 1:size(target_grid, 3), j in axes(mask, 2)
        mask[Nx+1, j, k] = active_cell(Nx, j, k, target_grid)
    end
    return mask
end

"""
    column_wet_levels(target_grid)

`(surface, levels)`: the surface-level water mask of `target_grid` as an `(i, j)` matrix, and the
number of wet cells in each column. Both come from the same `water_mask` that `prepare_forcing`
uses, so "wet" means exactly what it means when the forcing file is written, including
`PartialCellBottom` cells. A column's wet cells are contiguous and end at the surface, so
`levels[i, j]` also names where they start: `Nz - levels + 1`.
"""
function column_wet_levels(target_grid)
    mask = water_mask(target_grid, Center, Center, nothing)
    surface = mask[:, :, size(mask, 3)]
    levels = dropdims(sum(mask; dims = 3); dims = 3)
    return surface, levels
end

"""
    projected_target_nodes(longitude, latitude, source::ProjectedSourceGrid)

Project the target node grid into the source projection, flattened with longitude varying
fastest. The whole point cloud is transformed in a single GDAL call, as in
`FjordSim.Bathymetry.sample_bathymetry_points!`.
"""
function projected_target_nodes(longitude, latitude, source::ProjectedSourceGrid)
    point_count = length(longitude) * length(latitude)
    x = Vector{Float64}(undef, point_count)
    y = Vector{Float64}(undef, point_count)

    point = 1
    for φ in latitude, λ in longitude
        x[point] = λ
        y[point] = φ
        point += 1
    end

    ArchGDAL.importEPSG(4326; order = :trad) do source_srs
        ArchGDAL.importPROJ4(source.proj4) do target_srs
            ArchGDAL.createcoordtrans(source_srs, target_srs) do transform
                ArchGDAL.transform!(x, y, zeros(Float64, point_count), transform)
            end
        end
    end

    return x, y
end

"""
    source_validity(filepath, source_name)

Whether each source cell holds data, taken from the first time step. Source land cells and
the points the download masked out as being outside the requested box are `NaN`, and that
pattern is a property of the dataset rather than of the time step, so one read serves every
step.

Always `(NX, NY, levels)`: a surface variable with no depth axis (`zeta`, `ubar`, `vbar`, which the
boundary pipeline prepares) reads as a plane and is reshaped to one level, so `SourceFill`,
`fill_source!` and `set_source_field!` stay three-dimensional for every source variable.
"""
function source_validity(filepath, source_name)
    return NCDataset(filepath) do ds
        return isfinite.(as_source_slab(coalesce.(read_last_dimension(ds[source_name], 1), NaN32)))
    end
end

"""
    read_last_dimension(variable, index)

`variable` at `index` along its last (time) dimension, with `:` for every dimension before it —
`(X, Y, depth)` for a 3D source variable, `(X, Y)` for a surface one. The trailing index cannot be
written literally because the two ranks differ.
"""
read_last_dimension(variable, index) =
    variable[ntuple(_ -> Colon(), ndims(variable) - 1)..., index]

"""
    as_source_slab(data)

`data` as a `(NX, NY, levels)` array: a 3D slab unchanged, a 2D surface plane reshaped to a single
level. Everything downstream of the source read — `SourceFill`, `fill_source!`,
`set_source_field!`, the interpolation kernel — is written for three dimensions, and one reshape
here is cheaper than a rank parameter threaded through all of them.
"""
as_source_slab(data::AbstractArray{<:Any,3}) = data
as_source_slab(data::AbstractArray{<:Any,2}) = reshape(data, size(data, 1), size(data, 2), 1)

"""
    source_fill(valid)

Build the `SourceFill` for a source mask: the nearest valid cell per masked cell per level.
"""
function source_fill(valid)
    NX, NY, ND = size(valid)
    level_valid = [any(@view valid[:, :, level]) for level = 1:ND]
    any(level_valid) || error("Source data has no valid cell inside the requested region.")

    nearest = zeros(Int32, NX, NY, ND)
    for level = 1:ND
        level_valid[level] || continue
        nearest[:, :, level] = nearest_valid_map(valid, level)
    end

    return SourceFill(BitArray(valid), nearest, level_valid)
end

"""
    fill_source!(slab, mask_fill::SourceFill)

Replace every masked value in `slab` with the nearest valid value in the same level, and every
fully masked level with the deepest level that has data. Mutates and returns `slab`, which is then
free of `NaN` so `interpolate` can never propagate one.

Filling a level from the deepest level with data reproduces what the previous hand-written
vertical interpolation achieved by clamping: source levels below the fjord's own bathymetry carry
no information, so the deepest real value is the best available.
"""
function fill_source!(slab, mask_fill::SourceFill)
    NX, NY, ND = size(slab)
    size(mask_fill.valid) == (NX, NY, ND) || throw(DimensionMismatch("source fill does not match the slab"))
    flat = reshape(slab, NX * NY, ND)

    for level = 1:ND
        mask_fill.level_valid[level] || continue
        for j = 1:NY, i = 1:NX
            @inbounds mask_fill.valid[i, j, level] || (flat[i+(j-1)*NX, level] = flat[mask_fill.nearest[i, j, level], level])
        end
    end

    deepest = findlast(mask_fill.level_valid)
    for level = 1:ND
        mask_fill.level_valid[level] && continue
        @inbounds flat[:, level] .= @view flat[:, deepest]
    end

    return slab
end

"""
    nearest_valid_map(valid, level)

For every cell of source `level`, the linear index of the nearest valid cell in that level.

Computed by relaxing squared Euclidean distance out from all valid cells at once, propagating
the originating cell over 8 neighbors, which costs one near-linear pass per level. This plays the
role of the `scipy.ndimage.distance_transform_edt` nearest-fill used to prepare the Oslofjord
forcing, and of `NumericalEarth`'s `inpaint_mask!`, which cannot be used here — see `SourceFill`.
"""
function nearest_valid_map(valid, level)
    NX, NY, _ = size(valid)
    nearest = zeros(Int32, NX, NY)
    origin_i = zeros(Int32, NX, NY)
    origin_j = zeros(Int32, NX, NY)
    distance = fill(typemax(Int32), NX, NY)
    queue = Tuple{Int32,Int32}[]

    for j = 1:NY, i = 1:NX
        valid[i, j, level] || continue
        nearest[i, j] = i + (j - 1) * NX
        origin_i[i, j] = i
        origin_j[i, j] = j
        distance[i, j] = 0
        push!(queue, (i, j))
    end

    head = 1
    while head <= length(queue)
        i, j = queue[head]
        head += 1
        source_i = origin_i[i, j]
        source_j = origin_j[i, j]

        for dj = -1:1, di = -1:1
            di == 0 && dj == 0 && continue
            neighbor_i = i + di
            neighbor_j = j + dj
            1 <= neighbor_i <= NX && 1 <= neighbor_j <= NY || continue

            candidate = (neighbor_i - source_i)^2 + (neighbor_j - source_j)^2
            candidate < distance[neighbor_i, neighbor_j] || continue

            distance[neighbor_i, neighbor_j] = candidate
            origin_i[neighbor_i, neighbor_j] = source_i
            origin_j[neighbor_i, neighbor_j] = source_j
            nearest[neighbor_i, neighbor_j] = source_i + (source_j - 1) * NX
            push!(queue, (neighbor_i, neighbor_j))
        end
    end

    return nearest
end

"""
    source_level_coordinate(source::ProjectedSourceGrid, z)

The heights `z` (metres, negative down) as fractional positions on the source's level axis — the
vertical coordinate of `source_field_grid`, whose cell `k` is the `k`-th source level counted from
the deepest. A height between two levels maps linearly between their indices, so the trilinear
interpolation `interpolate_to_target!` does in the projected horizontal and this coordinate is
linear in depth between the bracketing source levels. A height beyond the shallowest or deepest
level maps beyond that index, where the halo repeats the end level.

The levels cannot be the grid's own vertical coordinate. An Oceananigans grid is specified by its
faces, with centres at face midpoints, and a source's `depth` values *are* centres of a strongly
non-uniform axis: solving for faces that put centres on them leaves one free parameter, which each
level bounds from one side, alternating with parity. The interval is about a metre wide for the
operational NorKyst's 16 levels and empty for the Norkyst-v3 hindcast's 25.
"""
function source_level_coordinate(source::ProjectedSourceGrid, z)
    heights = -reverse(source.depths)
    return [level_position(heights, height) for height in z]
end

function level_position(heights, height)
    k = clamp(searchsortedlast(heights, height), 1, length(heights) - 1)
    return k + (height - heights[k]) / (heights[k+1] - heights[k])
end

"""
    source_field_grid(source::ProjectedSourceGrid, architecture = CPU())

The source subset as a `RectilinearGrid` in projected meters, on `architecture`.

The projected coordinates are regular, so this is a legal `RectilinearGrid` and
`Oceananigans.Fields.interpolate` applies. It is *not* expressible as a `LatitudeLongitudeGrid`:
NorKyst's grid is rotated about 59 degrees from east in the Oslofjord region, which also rules
out the `NumericalEarth` dataset path, whose `native_grid` is always a `LatitudeLongitudeGrid`.

The vertical coordinate is the level index, with cell `k` centred on `k`; target heights reach it
through `source_level_coordinate`.

A source on a differently-shaped grid overloads this, and `projected_target_nodes`, on its own
source-grid type.
"""
function source_field_grid(source::ProjectedSourceGrid, architecture = CPU())
    Δx = source.x[2] - source.x[1]
    Δy = source.y[2] - source.y[1]

    return RectilinearGrid(
        architecture;
        size = (length(source.x), length(source.y), length(source.depths)),
        x = collect(range(source.x[1] - Δx / 2, step = Δx, length = length(source.x) + 1)),
        y = collect(range(source.y[1] - Δy / 2, step = Δy, length = length(source.y) + 1)),
        z = collect(range(1 // 2, step = 1, length = length(source.depths) + 1)),
        topology = (Bounded, Bounded, Bounded),
    )
end

"""
    set_source_field!(source_field, slab, mask_fill::SourceFill)

Load one source slab into `source_field`, filling its mask via `fill_source!` and flipping the
depth axis so it increases upwards like the grid. `slab` is mutated.
"""
function set_source_field!(source_field, slab, mask_fill::SourceFill)
    fill_source!(slab, mask_fill)
    set!(source_field, reverse(slab, dims = 3))
    fill_halo_regions!(source_field)

    return source_field
end

@kernel function _interpolate_to_target!(output, @Const(x), @Const(y), @Const(z), @Const(mask),
                                         source_field, source_grid, source_location)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        node = (x[i, j], y[i, j], z[k])
        output[i, j, k] = ifelse(
            mask[i, j, k],
            convert(Float32, interpolate(node, source_field, source_location, source_grid)),
            NaN32,
        )
    end
end

"""
    interpolate_to_target!(output, source_field, x, y, z, mask, architecture)

Interpolate `source_field` onto the target nodes `x`, `y`, `z`, writing `NaN` where `mask` says the
target cell is land. Every array must live on `architecture`.

The target nodes were projected into the source coordinate system by `projected_target_nodes`,
which is what lets `Oceananigans.Fields.interpolate` do the work: the source is a
`RectilinearGrid` in the same projected meters, so one trilinear call covers both the horizontal
and the vertical.

`launch!` takes its work size from `size(output)` and its device from `architecture`; the grid
argument is unused for a size-tuple work spec, so the source grid is passed for consistency.
"""
function interpolate_to_target!(output, source_field, x, y, z, mask, architecture)
    launch!(
        architecture,
        source_field.grid,
        size(output),
        _interpolate_to_target!,
        output,
        x,
        y,
        z,
        mask,
        source_field,
        source_field.grid,
        (Center(), Center(), Center()),
    )

    return output
end

"""
    SourceReader(filepath)

Holds the one source file currently open, reopening on demand. Output steps are written
in date order so consecutive reads almost always hit the same file; only an interpolated day
straddling a month boundary reopens, which happens at most once per gap.
"""
mutable struct SourceReader{D}
    filepath::String
    dataset::D
end

SourceReader(filepath::String) = SourceReader(filepath, open_source_file(filepath))

# The chunk cache each variable of a source file gets, and how many chunks it may hold. Sized for a
# downloaded month written with netCDF's default chunking: NorKyst-v3 hourly boundary months come as
# `(15, 22, 7, 180)` chunks, so one hour of a full-depth variable touches 64 compressed chunks — about
# 106 MB once inflated — that each hold 180 hours. Below that, every hour decompresses the same 64
# chunks again: 88 ms per slab measured, against 1.0 ms with the block cached.
const SOURCE_CHUNK_CACHE_BYTES = 256 * 2^20
const SOURCE_CHUNK_CACHE_SLOTS = 1009
const SOURCE_CHUNK_CACHE_PREEMPTION = 0.75

"""
    open_source_file(filepath)

Open a downloaded source file for reading one time record at a time, with a per-variable chunk cache
large enough to hold every chunk one record touches — see `SOURCE_CHUNK_CACHE_BYTES`.

netCDF takes a variable's cache from the process-wide default at the moment the file is opened, so
the default is raised for the open and put back straight after: every other file this process opens
keeps the library's own setting.
"""
function open_source_file(filepath)
    size, nelems, preemption = NCDatasets.get_chunk_cache()
    NCDatasets.set_chunk_cache(
        size = SOURCE_CHUNK_CACHE_BYTES,
        nelems = SOURCE_CHUNK_CACHE_SLOTS,
        preemption = SOURCE_CHUNK_CACHE_PREEMPTION,
    )
    try
        return NCDataset(filepath)
    finally
        NCDatasets.set_chunk_cache(; size, nelems, preemption)
    end
end

Base.close(reader::SourceReader) = close(reader.dataset)

"""
    source_slab(reader::SourceReader, record::SourceRecord, source_name)

One `(X, Y, depth)` slab of `source_name`, with masked points as `NaN`. A surface variable with no
depth axis comes back as a single level, per `as_source_slab`.
"""
function source_slab(reader::SourceReader, record::SourceRecord, source_name)
    if reader.filepath != record.filepath
        close(reader.dataset)
        reader.dataset = open_source_file(record.filepath)
        reader.filepath = record.filepath
    end

    plane = read_last_dimension(reader.dataset[source_name], record.index)
    return as_source_slab(Float32.(coalesce.(plane, NaN32)))
end

"""
    blended_slab(reader::SourceReader, step::ForcingTimeStep, source_name)

The source slab for `step`: a single record when the step is present in the download, otherwise
the linear blend of the records bracketing an interpolated day. Land stays `NaN` because it is
`NaN` in both records.
"""
function blended_slab(reader::SourceReader, step::ForcingTimeStep, source_name)
    lower = source_slab(reader, step.lower, source_name)
    iszero(step.weight) && return lower

    upper = source_slab(reader, step.upper, source_name)
    @. lower = (1 - step.weight) * lower + step.weight * upper

    return lower
end

"""
    define_forcing_dimensions!(ds, target_grid, dates)

Define the dimensions and coordinate variables `forcing_from_file` checks against the grid, on the
time axis `dates`.

Takes the dates rather than the `ForcingTimeStep`s they come from, because `write_river_file` writes
this same layout for a file that has no source records behind it.
"""
function define_forcing_dimensions!(ds, target_grid, dates)
    coordinates = (
        "Nx" => Array(λnodes(target_grid, Center())),
        "Ny" => Array(φnodes(target_grid, Center())),
        "Nz" => Array(znodes(target_grid, Center())),
        "Nx_faces" => Array(λnodes(target_grid, Face())),
        "Ny_faces" => Array(φnodes(target_grid, Face())),
        "Nz_faces" => Array(znodes(target_grid, Face())),
    )

    for (name, values) in coordinates
        defDim(ds, name, length(values))
    end
    defDim(ds, "time", length(dates))

    for (name, values) in coordinates
        defVar(ds, name, Float64, (name,))[:] = values
    end
    defVar(ds, "time", collect(dates), ("time",))

    return ds
end

end # module
