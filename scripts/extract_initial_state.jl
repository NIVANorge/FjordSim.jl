# Copy one record of a run's ocean snapshot file into a small file `FromResults` can start a new
# run from:
#
#     julia --project scripts/extract_initial_state.jl SNAPSHOTS DATE OUTPUT
#
# e.g. the 2015-04-01 state of the full Oslofjord validation run, as the start of a 2014 run:
#
#     julia --project scripts/extract_initial_state.jl \
#         ~/FjordSim_results/oslofjorden_validation/snapshots_ocean_20261014T142318.nc \
#         2015-04-01 ~/FjordSim_data/oslofjorden_validation/initial_state_20150401.nc
#
# and in the config, `initial_conditions = FromResults(joinpath(fjord_data_root(
# "oslofjorden_validation"), "initial_state_20150401.nc"))`, which reads the file's only record.
#
# `FromResults(SNAPSHOTS, DATE)` reads the same record straight from the snapshot file, which is
# enough for a local run. This copy is for a run on GCP: a fresh `gcp/fjordsim-gcp run` stages the
# top-level *.nc of the data directory, not an earlier run's results, so a snapshot file of tens of
# GB never reaches the VM. This file is ~100 MB and `push-data` uploads it like any prepared input.
# The record is chosen the way `FromResults` chooses it: the nearest one within half the snapshot
# interval, since an adaptive time step writes each a little late.
#
# The output keeps the snapshot's own dimension and variable names, which is what `results_state`
# reads. Every variable with a time axis is copied, so the turbulence field `e` rides along although
# `FromResults` reads only the tracers and `u`, `v`. Salinity is floored at zero: a run leaves a few
# cells at river mouths slightly negative, which no start state should carry.
using NCDatasets
using Dates: DateTime, Millisecond

const START_DATE_ATTRIBUTE = "start_date"

"The record of `ds` nearest `date`, checked to be within half the snapshot interval."
function nearest_record(ds, date)
    source_start = DateTime(ds.attrib[START_DATE_ATTRIBUTE])
    seconds = (date - source_start).value / 1000
    times = ds["time"][:]
    index = argmin(abs.(times .- seconds))
    offset = times[index] - seconds
    abs(offset) <= ds.attrib["interval"] / 2 || error(
        "No record within half an interval of $date; the file runs $source_start to " *
        "$(source_start + Millisecond(round(Int, 1000 * last(times)))).",
    )
    return index, source_start + Millisecond(round(Int, 1000 * times[index]))
end

function extract_initial_state(source, date, destination)
    isfile(destination) && error("$destination exists; delete it first to rewrite it.")

    NCDataset(source) do ds
        index, record_date = nearest_record(ds, date)
        names = [name for name in keys(ds) if "time" in dimnames(ds[name]) && name != "time"]
        spatial_dimensions = unique(dimension for name in names for dimension in dimnames(ds[name])
                                    if dimension != "time")

        NCDataset(destination, "c") do out
            out.attrib[START_DATE_ATTRIBUTE] = string(record_date)
            out.attrib["source"] = abspath(source)
            out.attrib["source_record"] = index
            out.attrib["source_start_date"] = ds.attrib[START_DATE_ATTRIBUTE]

            for dimension in spatial_dimensions
                defVar(out, dimension, ds[dimension][:], (dimension,); attrib = ds[dimension].attrib)
            end
            defVar(out, "time", [0.0], ("time",); attrib = ds["time"].attrib)

            for name in names
                data = ds[name][:, :, :, index]
                name == "S" && (data = max.(data, 0))
                defVar(out, name, reshape(data, size(data)..., 1), dimnames(ds[name]);
                       attrib = ds[name].attrib)
            end
        end

        @info "Wrote record $index ($record_date) of $source to $destination: $(join(names, ", "))"
    end
    return destination
end

length(ARGS) == 3 || error("usage: extract_initial_state.jl SNAPSHOTS DATE OUTPUT")
extract_initial_state(expanduser(ARGS[1]), DateTime(ARGS[2]), expanduser(ARGS[3]))
