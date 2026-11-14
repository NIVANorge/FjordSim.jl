# Measure the per-constituent gains that correct a setup's open-boundary tide to a Kartverket tide
# gauge, from the downloaded boundary data alone, before any run:
#
#     julia --project scripts/measure_boundary_tide.jl CONFIG GAUGE
#
# e.g. for the Oslofjord setup and Viker,
#
#     julia --project scripts/measure_boundary_tide.jl examples/oslofjorden/oslofjorden.jl Viker
#
# GAUGE is a Kartverket gauge of the setup's: one its tide-gauge writer records `η` at, which
# `validate_simulation` scores, or one of its `observation_configs`. Run it after
# `download_boundary_data`. It averages the downloaded `zeta` over every wet source point
# of the boundary band, fetches the gauge's 10-minute record over the same window (cached under the
# setup's data root), fits both with `harmonic_fit` on one reference time, and prints one
# `TidalGain(constituent, observed / boundary amplitude, observed - boundary phase)` per constituent,
# ready to paste into the setup. Averaging along the band is sound where the boundary tide is nearly
# uniform along the edge, and gauge-over-boundary stands for gauge-over-model where the model carries
# the boundary tide to the gauge unchanged; both hold for the Oslofjord (see its setup).
#
# A config whose boundary data config wraps another one in a `source` field, as the Oslofjord's tide
# correction does, is measured on that source.
using FjordSim
using FjordSim: BoundaryData, Configs
using FjordSim.Observations: KartverketSeaLevel, fetch_observations, observation_series, observation_stations
using FjordSim.Validation: default_observations, harmonic_fit
using NCDatasets
using Dates: DateTime

const CONSTITUENTS = ("M2", "S2", "N2", "K2", "O1", "K1")

"The band-mean `zeta` of every downloaded record, by time, duplicates dropped."
function boundary_sea_level(source)
    records = BoundaryData.boundary_time_steps(source)
    times = DateTime[]
    values = Float64[]
    for path in unique(record.filepath for record in records)
        NCDataset(path) do ds
            zeta = coalesce.(ds["zeta"][:, :, :], NaN)
            for (index, time) in enumerate(ds["time"][:])
                wet = filter(isfinite, view(zeta, :, :, index))
                isempty(wet) && continue
                push!(times, DateTime(time))
                push!(values, sum(wet) / length(wet))
            end
        end
    end
    keep = unique(i -> times[i], sortperm(times))
    return times[keep], values[keep]
end

gauges(source::KartverketSeaLevel) = observation_stations(source)
gauges(source) = Station[]

function main(arguments)
    length(arguments) == 2 || error("usage: measure_boundary_tide.jl CONFIG GAUGE")
    config_path, name = arguments
    config = fjord_config(config_path)
    boundary = Configs.boundary_data_config(config)
    source = hasproperty(boundary, :source) ? boundary.source : boundary

    sources = (
        default_observations(Configs.run_config(config), Configs.bathymetry_config(config))...,
        Configs.observation_configs(config)...,
    )
    found = [station for candidate in sources for station in gauges(candidate) if station.name == name]
    isempty(found) && error("$name is none of the setup's Kartverket gauges.")
    station = first(found)

    times, model = boundary_sea_level(source)
    window = (first(times), last(times))
    @info "Boundary sea level: $(length(times)) records, $(window[1]) to $(window[2])"

    gauge = KartverketSeaLevel(
        data_root = source.data_root,
        stations = [station],
        start_date = window[1],
        stop_date = window[2],
    )
    fetch_observations(nothing, gauge)
    observed = observation_series(gauge, station, "eta"; window)
    @info "$name: $(length(observed.times)) samples"

    reference = first(times)
    boundary_fit = harmonic_fit(times, model; reference)
    gauge_fit = harmonic_fit(observed.times, vec(observed.values); reference)

    println("[")
    for constituent in CONSTITUENTS
        k = findfirst(c -> c.name == constituent, boundary_fit.constituents)
        gain = gauge_fit.amplitude[k] / boundary_fit.amplitude[k]
        shift = mod(gauge_fit.phase[k] - boundary_fit.phase[k] + 180, 360) - 180
        println(
            "    TidalGain(\"$constituent\", $(round(gain; digits = 3)), $(round(shift; digits = 1))),",
            "   # $(round(100 * boundary_fit.amplitude[k]; digits = 2)) cm at the boundary, ",
            "$(round(100 * gauge_fit.amplitude[k]; digits = 2)) observed, ",
            "± $(round(100 * gauge_fit.amplitude_error[k]; digits = 2))",
        )
    end
    println("]")
    return nothing
end

main(ARGS)
