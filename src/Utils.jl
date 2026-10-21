module Utils

export progress,
    fjord_data_root,
    fjord_results_root

using Oceananigans
using Oceananigans.Fields: interior, location
using Oceananigans.OutputReaders: FieldTimeSeries, OnDisk
using Oceananigans.Utils: prettytime
using Oceananigans.Advection: cell_advection_timescale
using JLD2: @save
using Printf: @sprintf
using Dates: Second

using ..Configs: AbstractSimulationConfig

cell_advection_timescale_coupled_model(coupled_model) = cell_advection_timescale(coupled_model.ocean.model)

# Generic Base.getindex for any variable using Val types
# fields[:VARIABLE] doesn't work in CUDA kernels
@inline Base.getindex(fields, i, j, k, ::Val{name}) where {name} = @inbounds getfield(fields, name)[i, j, k]

function compute_faces(centers)
    spacing = diff(centers)[1]  # Assuming uniform spacing
    faces = vcat([centers[1] - spacing / 2], (centers[1:end-1] .+ centers[2:end]) / 2, [centers[end] + spacing / 2])
    return faces
end

const WALL_TIME = Ref(time_ns())

"""
    progress(sim)

The default `ProgressCallback` report: iteration, model time, time step, the per-component velocity
maxima and the temperature range.

Every reduction is taken on the `Field`, not on `interior(field)`. That is the whole difference
between a useful report and a misleading one: a reduction over a `Field` on an `ImmersedBoundaryGrid`
excludes the immersed periphery for free (`Oceananigans.ImmersedBoundaries.NotImmersed`), while
`interior` hands back a bare array and throws that away. Oceananigans writes `zero(eltype)` into
every immersed peripheral tracer cell at the top of each `update_state!`, and on a fjord grid most of
the block is land — so the array form reported the land mask as the minimum temperature, for the
whole run. It also keeps the reduction on the GPU rather than indexing a `CuArray` cell by cell.

Reaches `sim.model.ocean.model.tracers.T`, so a model whose tracers omit `:T` needs a `report` of its
own on its `ProgressCallback` rather than this one.
"""
function progress(sim)
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T

    Tmin = minimum(T)
    Tmax = maximum(T)

    umax = (maximum(abs, u), maximum(abs, v), maximum(abs, w))

    step_time = 1e-9 * (time_ns() - WALL_TIME[])

    msg = @sprintf("Iter: %d, time: %s, Δt: %s", iteration(sim), prettytime(sim), prettytime(sim.Δt))
    msg *= @sprintf(
        ", max|u|: (%.2e, %.2e, %.2e) m s⁻¹, extrema(T): (%.2f, %.2f) ᵒC, wall time: %s",
        umax...,
        Tmin,
        Tmax,
        prettytime(step_time)
    )

    @info msg

    WALL_TIME[] = time_ns()
end

function safe_execute(callable)
    return function (args...)
        if callable === nothing || args === nothing
            return nothing
        elseif isa(callable, Function)
            return callable(args...)
        else
            return nothing
        end
    end
end

function extract_z_faces(grid)
    bar = grid["zᵃᵃᶜ"]
    zero_index = findfirst(x -> x > 0.0, bar)
    n = grid["Nz"] + 1
    if zero_index > 1
        start_index = max(1, zero_index - n)
        z = bar[start_index:zero_index-1]
    else
        z = Float64[]
    end
    return z
end

function netcdf_to_jld2(netcdf_file::String, jld2_file::String)
    ds = NCDataset(netcdf_file, "r")
    data_dict = Dict()
    for variable_name in keys(ds)
        array = convert(Array, ds[variable_name])
        data_dict[variable_name] = array
        print(size(array))
    end

    @save jld2_file data_dict
    close(ds)
    println("Conversion completed: NetCDF to JLD2")
end

function save_fts(; jld2_filepath, fts_name, fts, grid, times, boundary_conditions)
    isfile(jld2_filepath) && rm(jld2_filepath)
    LX, LY, LZ = location(fts)
    on_disk_fts = FieldTimeSeries{LX,LY,LZ}(
        grid,
        times;
        boundary_conditions,
        backend = OnDisk(),
        path = jld2_filepath,
        name = fts_name,
    )
    for i = 1:size(fts)[end]
        set!(on_disk_fts, fts[i], i, times[i])
    end
end

function recursive_merge(a::NamedTuple, b::NamedTuple)
    # Get all unique keys from both NamedTuples
    all_keys = union(keys(a), keys(b))

    # Initialize an empty NamedTuple for the result
    result_pairs = Pair{Symbol,Any}[]

    for key in all_keys
        a_value = get(a, key, nothing)
        b_value = get(b, key, nothing)

        if a_value isa NamedTuple && b_value isa NamedTuple
            # If both values are NamedTuples, recursively merge them
            push!(result_pairs, key => recursive_merge(a_value, b_value))
        elseif b_value !== nothing
            # If only b_value exists or is not a NamedTuple, use b_value
            push!(result_pairs, key => b_value)
        else
            # Otherwise, use a_value (if it exists)
            push!(result_pairs, key => a_value)
        end
    end
    return (; result_pairs...)
end

"""
    recursive_merge(a, b, rest...)

Fold `recursive_merge` left to right over any number of named tuples, so later arguments win per
leaf exactly as they do for two. A boundary config contributing several independent groups — an open
edge contributes one per part of the state — merges them in one call rather than nesting pairs.
"""
recursive_merge(a::NamedTuple, b::NamedTuple, rest::NamedTuple...) =
    foldl(recursive_merge, rest; init = recursive_merge(a, b))


"""
    coverage_window(config)

The calendar interval a run needs its prepared inputs to span, as `(first, last)` — or `nothing`
for a setup naming no simulation config.

Passed to `prepare_forcing` and `prepare_atmosphere` as their `coverage`, so each pads its own time
axis to reach both ends. A plain tuple of dates rather than the config itself, so neither pipeline
ever sees a simulation config.

A looped run reuses one window rather than extending it, so this is `stop_time` from `start_date`
and not `loops * stop_time`.
"""
coverage_window(::Nothing) = nothing

coverage_window(config::AbstractSimulationConfig) =
    (config.start_date, config.start_date + Second(round(Int, config.stop_time)))

"""
    fjord_data_root(fjord)
    fjord_results_root(fjord)

The directory a setup keeps `fjord`'s input data in, and the one it writes that fjord's results to.

The environment names the *parent* — `FJORDSIM_DATA_ROOT` and `FJORDSIM_RESULTS_ROOT` — and the
argument names the fjord under it, so one variable relocates every setup at once and a setup still
states which fjord it is. Unset, they reproduce `~/FjordSim_data/<fjord>/` and
`~/FjordSim_results/<fjord>/` exactly, which is what every setup spelled with `homedir()` before.

This is what lets the same setup run against a laptop's home directory and a cloud VM's staging
disk without editing the setup: the roots are the only thing that differs between the two, since
every other path resolves against them.

A fjord that reads another's downloads passes that other name — `examples/drammensfjorden.jl` shares
Oslofjord's bathymetry source as `fjord_data_root("oslofjorden")` — so a relocation moves both.
"""
fjord_data_root(fjord) =
    joinpath(get(ENV, "FJORDSIM_DATA_ROOT", joinpath(homedir(), "FjordSim_data")), fjord)

fjord_results_root(fjord) =
    joinpath(get(ENV, "FJORDSIM_RESULTS_ROOT", joinpath(homedir(), "FjordSim_results")), fjord)

end # module
