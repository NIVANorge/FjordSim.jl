# Operational Norkyst-v3 adapter for the initial state: one hourly instant of a 160 m nest, merged
# into an 800 m member that fills the rest of the domain. The download does the merge, so what it
# writes is one file on one regular 160 m projected grid and the generic regridding reads it as it
# reads any NorKyst subset.

using Dates: Date, DateTime, @dateformat_str
import Dates
using Oceananigans.Grids: x_domain, y_domain

using ..Regridding: nested_slab
using ..NorKyst:
    NORKYST_V3_OPENDAP_URL,
    NORKYST_VARIABLE_NAMES,
    NORKYST_PROJECTION_VARIABLE,
    norkyst_source_grid,
    norkyst_v3_url,
    record_index,
    open_v3_file,
    read_remote,
    subset_ranges,
    concrete_float_data,
    decoded_attributes,
    define_subset_variable!

# The two resolutions a nested member is published at.
const NORKYST_V3_NEST_RESOLUTION = "160"
const NORKYST_V3_PARENT_RESOLUTION = "800"

"""
    NorKystV3InitialStateConfig

The initial state from the operational Norkyst-v3 (from 2024-01-01) at `start_date`, which must be
one of its hourly instants: a 160 m nest wherever the nest reaches, and an 800 m member everywhere
else in the domain box.

The download merges them on the nest's own 160 m grid, extended over the 800 m member's subset:
inside the nest every node is a nest node and its coastline decides what is water; outside, the
800 m member is interpolated bilinearly. The nest's own 800 m parent, the same member at 800 m and
two-way coupled to it, makes the seam a change of resolution only; the Oslofjord nest's (`m71`) is
published from 2024-11-23, and before that the reference run `m00` is the fill.

# Fields
- `data_root`: Directory holding this setup's initial-state files. Required.
- `output_directory`: Name of the directory the merged download is written to. Required.
- `output_file`: Name of the prepared NetCDF written by `prepare_initial_state`.
- `opendap_url`: The collection's public endpoint. Defaulted.
- `nest_member`: The member whose 160 m nest is read: `"m70"` (Sulafjorden) or `"m71"`
  (Oslofjord). Required.
- `parent_member`: The member whose 800 m file fills outside the nest. Required.
- `parameters`: Variable names to extract (e.g. `["temperature", "salinity"]`). Required.
- `start_date`: The instant the state is for, independent of the run's; a run starting from it
  with `FromInitialState` must start at the same instant. Required.
"""
Base.@kwdef mutable struct NorKystV3InitialStateConfig <: AbstractInitialStateConfig
    data_root::String
    output_directory::String
    output_file::String = "initial_state.nc"
    opendap_url::String = NORKYST_V3_OPENDAP_URL
    nest_member::String
    parent_member::String
    parameters::Vector{String}
    start_date::DateTime
end

"""
    initial_state_filename(config::NorKystV3InitialStateConfig, date)

Name of the merged NetCDF downloaded for the instant `date`.
"""
initial_state_filename(config::NorKystV3InitialStateConfig, date) =
    "Norkyst-v3-operational_$(config.nest_member)_in_$(config.parent_member)_" *
    "$(Dates.format(date, dateformat"yyyymmddTHH")).nc"

initial_state_variable_names(config::NorKystV3InitialStateConfig) = NORKYST_VARIABLE_NAMES

initial_state_source_grid(config::NorKystV3InitialStateConfig, filepath) = norkyst_source_grid(filepath)

function initial_state_record(config::NorKystV3InitialStateConfig, date)
    filepath = joinpath(initial_state_directory(config), initial_state_filename(config, date))
    isfile(filepath) || error(
        "No Norkyst-v3 initial state for $date at $filepath. " *
        "Run `julia --project -m FjordSim download_initial_state` for this setup first.",
    )

    return NCDataset(filepath) do ds
        SourceRecord(DateTime(first(ds["time"][:])), filepath, 1)
    end
end

"""
    fetch_initial_state(target_grid, config::NorKystV3InitialStateConfig, date)

Download the record at `date` of the 160 m nest and of the 800 m member filling outside it, the
latter subset to the lon/lat box of `target_grid`, and write them merged by `nested_slab` into
`initial_state_directory(config)`. An instant already downloaded is skipped.
"""
function fetch_initial_state(target_grid, config::NorKystV3InitialStateConfig, date)
    output_path = joinpath(initial_state_directory(config), initial_state_filename(config, date))
    if isfile(output_path)
        @info "Skipping $date (already exists)"
        return output_path
    end
    mkpath(initial_state_directory(config))

    parent_url = norkyst_v3_url(config.opendap_url, NORKYST_V3_PARENT_RESOLUTION, config.parent_member, date)
    nest_url = norkyst_v3_url(config.opendap_url, NORKYST_V3_NEST_RESOLUTION, config.nest_member, date)
    @info "Downloading the Norkyst-v3 state at $date from $nest_url and $parent_url"

    open_v3_file(parent_url, date) do parent
        open_v3_file(nest_url, date) do nest
            projection(ds) = NCDatasets.variable(ds, NORKYST_PROJECTION_VARIABLE).attrib["proj4"]
            projection(parent) == projection(nest) ||
                error("The nest and its parent are on different projections.")
            parent["depth"][:] == nest["depth"][:] ||
                error("The nest and its parent have different depth levels.")

            subset = subset_ranges(parent, x_domain(target_grid), y_domain(target_grid), config.parameters)
            coarse_x = Float64.(parent["X"][subset.ranges["X"]])
            coarse_y = Float64.(parent["Y"][subset.ranges["Y"]])
            fine_x = Float64.(nest["X"][:])
            fine_y = Float64.(nest["Y"][:])
            x = collect(range(first(coarse_x), last(coarse_x); step = fine_x[2] - fine_x[1]))
            y = collect(range(first(coarse_y), last(coarse_y); step = fine_y[2] - fine_y[1]))

            parent_index = record_index(parent, date, parent_url)
            nest_index = record_index(nest, date, nest_url)
            @info "  Merged grid $(length(x)) x $(length(y)) at $(fine_x[2] - fine_x[1]) m"

            # Written under a temporary name and renamed once complete, so an interrupted download
            # is redone on the next run rather than skipped as already present.
            temporary_path = output_path * ".tmp"
            NCDataset(temporary_path, "c") do output
                defDim(output, "X", length(x))
                defDim(output, "Y", length(y))
                defDim(output, "depth", NCDatasets.dim(parent, "depth"))
                defDim(output, "time", 1)
                defVar(output, "X", x, ("X",); attrib = decoded_attributes(parent["X"]))
                defVar(output, "Y", y, ("Y",); attrib = decoded_attributes(parent["Y"]))
                defVar(output, "time", [DateTime(date)], ("time",))
                define_subset_variable!(output, parent, "depth", subset; deflatelevel = 0)
                define_subset_variable!(output, parent, NORKYST_PROJECTION_VARIABLE, subset; deflatelevel = 0)

                # Every variable is (X, Y, depth, time) in both files.
                for name in config.parameters
                    coarse = concrete_float_data(read_remote(
                        parent[name], (subset.ranges["X"], subset.ranges["Y"], :, parent_index), parent,
                    ))
                    fine = concrete_float_data(read_remote(nest[name], (:, :, :, nest_index), nest))
                    merged = nested_slab(coarse, coarse_x, coarse_y, fine, fine_x, fine_y, x, y)
                    variable = defVar(output, name, Float32, ("X", "Y", "depth", "time");
                        deflatelevel = 5, attrib = decoded_attributes(parent[name]))
                    variable[:, :, :, 1] = merged
                end
            end
            mv(temporary_path, output_path; force = true)
        end
    end

    @info "Initial-state download saved to $output_path"
    return output_path
end
