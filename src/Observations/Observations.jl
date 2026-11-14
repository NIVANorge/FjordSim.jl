"""
    Observations

Observed data about a fjord: the sources, the `download_observations` step that fetches them over
the setup's domain, and the one record type every source returns. `validate_simulation` scores a run
against them; a setup can also start its initial state from them or measure its boundary tide
against them.

The same adapter shape every other dataset in FjordSim uses: an abstract config supertype, one
subtype per source, and hooks dispatching on it, so adding a programme is a new subtype and its
methods, never an edit to the code that reads them.

Four files below it:

- `stations.jl` — `Station`, and the files every source writes and reads: `stations.csv` and one
  `<station tag>_<variable>.csv` per station and variable.
- `kartverket.jl` — Kartverket's tide gauges.
- `csv.jl` — `CsvObservations`, files someone else downloaded.
- `vannmiljo.jl` — the Norwegian Environment Agency's Vannmiljø database.
"""
module Observations

export Station, download_observations

using Dates: DateTime, Date
import Dates
using Downloads
using JSON
using ArchGDAL
using Printf: @sprintf
using Oceananigans.Architectures: CPU

using ..Configs: AbstractFjordConfig, AbstractObservationConfig, grid_config, observation_configs
using ..Utils: months_in
using ..Grids: domain_grid, download_box

observations_directory(config::AbstractObservationConfig) =
    joinpath(config.data_root, config.output_directory)

"""
    fetch_observations(target_grid, config)

Download everything `config` can supply inside `target_grid`'s longitude/latitude box over the
source's own period.

Defaults to a no-op: a source read from files downloaded some other way has nothing to fetch, and
should not have to say so.
"""
fetch_observations(target_grid, ::AbstractObservationConfig) = nothing

"""
    observation_stations(config)

The `Station`s this source can supply, as positions rather than names, so the same lon/lat that
placed the model's output places the observation beside it.
"""
function observation_stations end

"""
    observation_series(config, station, variable; window)

One station's record of one variable over `window`, as an `ObservationSeries`, or `nothing` when the
source has nothing for that pair.

Returning `nothing` rather than raising is deliberate. A validation sweep asks every source for every
station and variable it might hold; most combinations are empty, and that is not an error.
"""
function observation_series end

"""
    ObservationSeries

One station's observations of one variable: a time axis, a depth axis, and a value per pair.

`values` is `(length(times), length(depths))`. A surface or point record has `depths == [0.0]` and
one column. A profile record — a CTD programme, an ADCP — has one row per cast or sample time.

**Profiles are interpolated onto one depth axis at read time.** Instruments do not sample the same
depths twice: a CTD cast returns whatever the winch gave, and an ADCP's bins depend on its range.
Carrying that through would make every comparison a special case, and a Hovmøller diagram needs a
rectangular array anyway. The interpolation is linear and does not extrapolate — a depth outside a
given cast's range is `NaN` for that cast, not the nearest value — so a shallow cast cannot invent
deep water.

# Fields
- `station`: the station name, as the source calls it.
- `variable`: the FjordSim name, matching the model field it is scored against (`"T"`, `"S"`, `"u"`,
  `"v"`, `"eta"`).
- `times`: one `DateTime` per row.
- `depths`: metres below the surface, positive down, **ascending**, one per column.
  `depth_column` and `match_series` both bisect it, so an unsorted axis would silently mis-sample.
- `values`: the observations, `NaN` where missing.
- `units`, `source`: provenance, carried so a figure can be labelled without the config.
"""
struct ObservationSeries
    station::String
    variable::String
    times::Vector{DateTime}
    depths::Vector{Float64}
    values::Matrix{Float64}
    units::String
    source::String
end

"""
    interpolate_to_depths(sample_depths, sample_values, depths)

Linear interpolation of one cast onto a common depth axis, with `NaN` outside the cast's own range.
Used by every profile source so they all present the same shape.
"""
function interpolate_to_depths(sample_depths, sample_values, depths)
    order = sortperm(sample_depths)
    d = Float64[sample_depths[i] for i in order if isfinite(sample_values[i])]
    v = Float64[sample_values[i] for i in order if isfinite(sample_values[i])]

    output = fill(NaN, length(depths))
    length(d) >= 2 || return output

    for (index, depth) in enumerate(depths)
        (depth < first(d) || depth > last(d)) && continue
        upper = searchsortedfirst(d, depth)
        if upper == 1 || d[upper] == depth
            output[index] = v[upper]
        else
            lower = upper - 1
            weight = (depth - d[lower]) / (d[upper] - d[lower])
            output[index] = (1 - weight) * v[lower] + weight * v[upper]
        end
    end

    return output
end

include("stations.jl")
include("kartverket.jl")
include("csv.jl")
include("vannmiljo.jl")

"""
    download_observations(config)

Fetch every source in the setup's `observation_configs` over its domain: the longitude/latitude box
of `domain_grid`, each source over its own period. Only picks out the grid and the sources and hands
them to `fetch_observation_sources`.
"""
download_observations(config::AbstractFjordConfig) =
    fetch_observation_sources(domain_grid(grid_config(config), CPU()), observation_configs(config))

"""
    fetch_observation_sources(target_grid, sources)

`fetch_observations` for each of `sources` in turn. A source that fails — no network, no key, a
service down — is warned about and the next one is tried: each source's files stand on their own.
"""
function fetch_observation_sources(target_grid, sources)
    for source in sources
        try
            fetch_observations(target_grid, source)
        catch exception
            exception isa InterruptException && rethrow()
            @warn "Could not fetch $(nameof(typeof(source)))" exception
        end
    end
    return nothing
end

end # module
