# The run's 3D snapshots, read back for what the station output cannot give: a column at any station,
# whether or not a station writer recorded it, and the state of the whole domain at every record.
#
# Both come from one pass over the files. NetCDF stores a snapshot file one record per chunk, so a
# column's time series read straight from the file reads every chunk once per station.

"""
    snapshot_files(directory, stem)

Every file a `SnapshotWriter` with output stem `stem` wrote into `directory`, in the order their
records were written: by run tag, then loop, then part. A file that was never split counts as part 1.
"""
function snapshot_files(directory, stem::AbstractString)
    isdir(directory) || return String[]

    pattern = Regex("^" * stem * "_(\\d{8}T\\d{6})(?:_loop(\\d+))?(?:_part(\\d+))?\\.nc\$")
    found = Tuple{String,Int,Int,String}[]
    for file in readdir(directory)
        match_result = match(pattern, file)
        isnothing(match_result) && continue
        tag, loop, part = match_result.captures
        push!(
            found,
            (
                tag,
                isnothing(loop) ? 0 : parse(Int, loop),
                isnothing(part) ? 1 : parse(Int, part),
                joinpath(directory, file),
            ),
        )
    end

    return last.(sort!(found))
end

"""
    snapshot_records(files)

`(file, record, time)` for every record of `files`, in time order, with each time read once.

A resumed run writes under a new run tag from its checkpoint on, and each loop of a looped run starts
its model time again, so where two files hold the same time it is the later file that holds the
run's record of it: each file replaces every earlier record from its own first time on.
"""
function snapshot_records(files)
    records = Tuple{Int,Int,DateTime}[]
    for (file, path) in enumerate(files)
        times = NCDataset(path) do ds
            reference = station_reference_date(Dict{String,Any}(ds.attrib), path)
            [reference + Dates.Millisecond(round(Int, 1000 * s)) for s in ds["time"][:]]
        end
        isempty(times) && continue
        filter!(record -> record[3] < first(times), records)
        append!(records, ((file, index, time) for (index, time) in enumerate(times)))
    end
    return records
end

"""
    snapshot_grid(path)

The geometry a snapshot file carries, as a `NamedTuple`: cell-centre `longitudes` and `latitudes`,
the domain's face bounds, cell-centre `depths` (positive down, in the file's bottom-to-top order),
each column's `water_depth`, the `wet` cells and each cell's `volume`.
"""
function snapshot_grid(path)
    return NCDataset(path) do ds
        λ_faces = Float64.(ds["λ_faa"][:])
        φ_faces = Float64.(ds["φ_afa"][:])
        wet = Array(ds["inactive_nodes_ccc"][:, :, :]) .== 0
        area = Float64.(ds["Δx_cca"][:, :]) .* Float64.(ds["Δy_cca"][:, :])
        thickness = Float64.(ds["Δz_aac"][:])
        (
            longitudes = Float64.(ds["λ_caa"][:]),
            latitudes = Float64.(ds["φ_aca"][:]),
            longitude_bounds = (first(λ_faces), last(λ_faces)),
            latitude_bounds = (first(φ_faces), last(φ_faces)),
            depths = -Float64.(ds["z_aac"][:]),
            water_depth = -Float64.(ds["bottom_height"][:, :]),
            wet,
            volume = [area[i, j] * thickness[k] for i in axes(wet, 1), j in axes(wet, 2), k in axes(wet, 3)],
            dimensions = NCDatasets.dimnames(ds["inactive_nodes_ccc"]),
        )
    end
end

"""
    snapshot_column(grid, station, deepest, search_radius)

The wet column a station's observations are scored at, as `(i, j, offset)` with the offset in grid
cells, or `nothing` when the station is outside the grid or has no wet column within
`search_radius` cells.

The nearest wet column is often the wrong one in a fjord: a cast to 300 m in a narrow basin sits beside
cells a fraction of that deep. So each wet column within `search_radius` (a disc, as for
`StationWriter`) is scored by the share of the deepest observation, `deepest` metres, it is too
shallow to reach, plus its distance as a share of `search_radius`, and the lowest score wins.
Missing the whole cast costs as much as lying at the edge of the search. A surface record,
`deepest == 0`, takes the nearest wet column.
"""
function snapshot_column(grid, station::Station, deepest, search_radius)
    (; longitudes, latitudes, water_depth, wet) = grid
    inside =
        grid.longitude_bounds[1] <= station.longitude <= grid.longitude_bounds[2] &&
        grid.latitude_bounds[1] <= station.latitude <= grid.latitude_bounds[2]
    inside || return nothing

    i₀ = argmin(abs.(longitudes .- station.longitude))
    j₀ = argmin(abs.(latitudes .- station.latitude))
    top = size(wet, 3)

    best = nothing
    for dj = -search_radius:search_radius, di = -search_radius:search_radius
        i, j = i₀ + di, j₀ + dj
        checkbounds(Bool, water_depth, i, j) && wet[i, j, top] || continue
        offset = sqrt(di^2 + dj^2)
        offset <= search_radius + 1//2 || continue

        shortfall = deepest > 0 ? max(deepest - water_depth[i, j], 0) / deepest : 0.0
        score = shortfall + offset / max(search_radius, 1)
        if isnothing(best) || score < best.score
            best = (; i, j, offset, score)
        end
    end

    return isnothing(best) ? nothing : (best.i, best.j, best.offset)
end

"""
    read_snapshots(files, records, grid, columns, variables, ranges)

One pass over `records` of `files`. Returns the `(time, depth)` series of each of `variables` at each
of `columns`, as `series[variable][:, :, column]` shallowest level first and `NaN` in dry cells, and
one row of domain statistics per record and variable (`domain_extremes.csv`), counting the wet cells
outside that variable's `(low, high)` in `ranges`.
"""
function read_snapshots(files, records, grid, columns, variables, ranges)
    series = Dict(v => fill(NaN, length(records), length(grid.depths), length(columns)) for v in variables)
    statistics = DomainStatistics[]
    total_volume = sum(grid.volume[grid.wet])

    for (file, path) in enumerate(files)
        NCDataset(path) do ds
            for (row, (record_file, record, time)) in enumerate(records)
                record_file == file || continue
                for variable in variables
                    # The stored values, fill values and all: only the wet cells are read.
                    field = ds[variable].var[:, :, :, record]
                    copy_columns!(series[variable], field, grid.wet, columns, row)
                    low, high = get(ranges, variable, (-Inf, Inf))
                    push!(statistics, domain_statistics(time, variable, field, grid, total_volume, low, high))
                end
                row % 30 == 0 && @info "  snapshots read to $time"
            end
        end
    end

    return series, statistics
end

# One record of each column into `series[row, :, column]`, shallowest level first.
function copy_columns!(series, field, wet, columns, row)
    Nz = size(field, 3)
    for (c, (i, j)) in enumerate(columns), k = 1:Nz
        level = Nz + 1 - k
        wet[i, j, level] && (series[row, k, c] = field[i, j, level])
    end
    return series
end

const DomainStatistics = @NamedTuple{
    time::DateTime,
    variable::String,
    mean::Float64,
    minimum::Float64,
    minimum_at::NTuple{3,Float64},
    maximum::Float64,
    maximum_at::NTuple{3,Float64},
    below::Int,
    above::Int,
}

"""
    domain_statistics(time, variable, field, grid, total_volume, low, high)

The volume mean of `field` over the wet cells, its extremes with the `(longitude, latitude, depth)`
of each, and how many wet cells lie below `low` and above `high`. A non-finite value in a wet cell
makes the mean `NaN`, which is what a blow-up should look like.
"""
function domain_statistics(time, variable, field, grid, total_volume, low, high)
    total = 0.0
    smallest, largest = Inf, -Inf
    smallest_at = largest_at = first(CartesianIndices(field))
    below = above = 0

    for index in CartesianIndices(field)
        grid.wet[index] || continue
        x = Float64(field[index])
        total += x * grid.volume[index]
        x < smallest && (smallest = x; smallest_at = index)
        x > largest && (largest = x; largest_at = index)
        below += x < low
        above += x > high
    end

    position(index) = (grid.longitudes[index[1]], grid.latitudes[index[2]], grid.depths[index[3]])
    return DomainStatistics((
        time, variable, total / total_volume,
        smallest, position(smallest_at), largest, position(largest_at), below, above,
    ))
end

"""
    write_domain_statistics(path, statistics)

`domain_extremes.csv`: one row per snapshot record and variable.
"""
function write_domain_statistics(path, statistics)
    open(path, "w") do io
        println(
            io,
            "time,variable,mean,minimum,minimum_longitude,minimum_latitude,minimum_depth," *
            "maximum,maximum_longitude,maximum_latitude,maximum_depth,cells_below,cells_above",
        )
        for s in statistics
            println(
                io,
                join(
                    (
                        Dates.format(s.time, "yyyy-mm-ddTHH:MM:SS"), s.variable,
                        @sprintf("%.4f", s.mean), @sprintf("%.4f", s.minimum),
                        @sprintf("%.5f", s.minimum_at[1]), @sprintf("%.5f", s.minimum_at[2]),
                        @sprintf("%.2f", s.minimum_at[3]), @sprintf("%.4f", s.maximum),
                        @sprintf("%.5f", s.maximum_at[1]), @sprintf("%.5f", s.maximum_at[2]),
                        @sprintf("%.2f", s.maximum_at[3]), s.below, s.above,
                    ),
                    ",",
                ),
            )
        end
    end
    return path
end

"""
    snapshot_series(run, validation, unrecorded, directory)

The observations no station writer recorded, each beside its column of the run's `SnapshotWriter`
output, as `(observed, modelled, tolerance)`; and, as a side effect, the domain statistics of every
snapshot record in `domain_extremes.csv` and `domain_extremes.png`, and where each station was placed
in `snapshot_columns.csv`.

`unrecorded` holds `(station, observed)` pairs. Only the variables the snapshots hold at cell
centres can be scored (`T`, `S`): a velocity is on a cell face. `tolerance` is half the writer's
interval plus an hour, since a record lands on the first time step past its scheduled time.

With `validation` `nothing` there are no snapshots to read, and every unrecorded observation is
reported as having no model output.
"""
function snapshot_series(run, ::Nothing, unrecorded, directory)
    for (station, observed) in unrecorded
        @info "  no model output for $(station.name) $(observed.variable)"
    end
    return Tuple{ObservationSeries,ModelSeries,Dates.Millisecond}[]
end

function snapshot_series(run, validation::ValidationConfig, unrecorded, directory)
    found = Tuple{ObservationSeries,ModelSeries,Dates.Millisecond}[]
    writer, files = snapshot_output(run)
    if isempty(files)
        @warn "No SnapshotWriter output in $(results_root(run)), so no snapshot to score or summarise."
        return snapshot_series(run, nothing, unrecorded, directory)
    end

    records = snapshot_records(files)
    if isempty(records)
        @warn "The SnapshotWriter files in $(results_root(run)) hold no record yet."
        return snapshot_series(run, nothing, unrecorded, directory)
    end
    grid = snapshot_grid(last(files))
    variables = NCDataset(last(files)) do ds
        [v for v in VALIDATION_VARIABLES if haskey(ds, v) && NCDatasets.dimnames(ds[v])[1:3] == grid.dimensions]
    end
    @info "Reading $(length(records)) snapshot records of $(join(variables, ", ")) from $(length(files)) files"

    # Every station's column, chosen once from the deepest of its observations of any variable.
    deepest = Dict{Station,Float64}()
    for (station, observed) in unrecorded
        reached = [observed.depths[k] for k in axes(observed.values, 2) if any(isfinite, view(observed.values, :, k))]
        deepest[station] = max(get(deepest, station, 0.0), isempty(reached) ? 0.0 : maximum(reached))
    end
    placed = Dict(station => snapshot_column(grid, station, depth, validation.search_radius) for (station, depth) in deepest)
    write_snapshot_columns(joinpath(directory, "snapshot_columns.csv"), placed, deepest, grid)

    columns = unique(cell[1:2] for cell in values(placed) if !isnothing(cell))
    series, statistics = read_snapshots(files, records, grid, columns, variables, validation.plausible_ranges)
    write_domain_statistics(joinpath(directory, "domain_extremes.csv"), statistics)
    plot_domain_extremes(joinpath(directory, "domain_extremes.png"), statistics, validation.plausible_ranges)

    times = last.(records)
    depths = reverse(grid.depths)
    tolerance = Dates.Millisecond(round(Int, 500 * writer.interval)) + Dates.Hour(1)
    for (station, observed) in unrecorded
        cell = placed[station]
        if isnothing(cell) || !(observed.variable in variables)
            @info "  no model output for $(station.name) $(observed.variable)"
            continue
        end
        column = findfirst(==(cell[1:2]), columns)
        metadata = Dict{String,Any}(
            "i" => cell[1], "j" => cell[2], "offset_cells" => cell[3],
            "water_depth" => grid.water_depth[cell[1], cell[2]],
        )
        modelled = ModelSeries(
            station.name, observed.variable, times, depths, series[observed.variable][:, :, column], metadata,
        )
        push!(found, (observed, modelled, tolerance))
    end

    return found
end

"""
    snapshot_output(run)

The first `SnapshotWriter` the run names and its files, or `(nothing, String[])` when it names none.
"""
function snapshot_output(run)
    for writer in writers(run)
        stem = snapshot_stem(writer)
        isnothing(stem) || return writer, snapshot_files(results_root(run), stem)
    end
    return nothing, String[]
end

# The stem of a `SnapshotWriter`'s files; any other writer writes no 3D snapshots.
snapshot_stem(writer::SnapshotWriter) = first(splitext(writer.output_file))
snapshot_stem(writer) = nothing

"""
    write_snapshot_columns(path, placed, deepest, grid)

`snapshot_columns.csv`: where each station's snapshot column is, and how deep, beside how deep its
observations go.
"""
function write_snapshot_columns(path, placed, deepest, grid)
    open(path, "w") do io
        println(io, "station,longitude,latitude,deepest_observation,i,j,column_longitude,column_latitude,water_depth,offset_cells")
        for station in sort(collect(keys(placed)); by = s -> s.name)
            cell = placed[station]
            head = @sprintf("%s,%.5f,%.5f,%.1f", station.name, station.longitude, station.latitude, deepest[station])
            if isnothing(cell)
                println(io, head, ",,,,,,")
            else
                i, j, offset = cell
                println(
                    io, head,
                    @sprintf(",%d,%d,%.5f,%.5f,%.1f,%.2f", i, j, grid.longitudes[i], grid.latitudes[j], grid.water_depth[i, j], offset),
                )
            end
        end
    end
    return path
end
