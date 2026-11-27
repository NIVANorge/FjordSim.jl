# The Oslofjord setup, 160 m, October 2024 to October 2026. What it covers, how to prepare and run
# it, and why its values are what they are: README.md beside this file.

using FjordSim
using Oceananigans:
    CATKEVerticalDiffusivity,
    HydrostaticSphericalCoriolis,
    SeawaterBuoyancy,
    WENO,
    WENOVectorInvariant
using Oceananigans.TurbulenceClosures: HorizontalScalarBiharmonicDiffusivity
using Oceananigans.Units: day, days, hour, hours, minutes, second
using Dates: Date, DateTime, Dates, Day, Hour, Month, Second
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using FjordSim: BoundaryData, Grids, InitialState, Plotting
using FjordSim.Configs: AbstractBoundaryDataConfig, AbstractInitialStateConfig,
    AbstractObservationConfig
using FjordSim.Observations: VannmiljoObservations, fetch_observations, observation_series,
    observation_stations
using FjordSim.Validation: TIDAL_CONSTITUENTS, ValidationConfig
using LinearAlgebra: qr
using NCDatasets: NCDataset
using NumericalEarth: FreezingLimitedOceanTemperature
using Printf: @sprintf
using Statistics: mean, median

# --- Stations ---------------------------------------------------------------------------------
#
# Each group states where its positions come from: a published table is exact, prose or a map is not.

# Kartverket's permanent tide gauges, from its station register. Exact.
const OSLOFJORD_TIDE_GAUGES = [
    Station(name = "Viker",      longitude = 10.949769, latitude = 59.036046),
    Station(name = "Oscarsborg", longitude = 10.604861, latitude = 59.678073),
    Station(name = "Oslo",       longitude = 10.734510, latitude = 59.908559),
]

# The ADCP moorings of September-November 2014, METreport 11/2017 Table 1, with its depths. Exact.
const OSLOFJORD_MOORINGS = [
    Station(name = "Km1", longitude = 10.627372, latitude = 59.582064),  # Filtvedt,       153 m
    Station(name = "Kn2", longitude = 10.646087, latitude = 59.581803),  # Brenntangen,     54 m
    Station(name = "Ri1", longitude = 10.497661, latitude = 59.350124),  # Småskjær,        20 m
    Station(name = "Rl1", longitude = 10.581023, latitude = 59.343452),  # Laksetrappa,     75 m
    Station(name = "Rm1", longitude = 10.626822, latitude = 59.352375),  # Botnegrunnen,    96 m
    Station(name = "Rn1", longitude = 10.653576, latitude = 59.363182),  # Evje,            64 m
    # ExxonMobil's ADCP at the Slagen refinery, from the report's prose. APPROXIMATE.
    Station(name = "Slagen", longitude = 10.517, latitude = 59.383),
]

# NIVA's Ytre Oslofjord stations, METreport 11/2017 Table 2. Exact. AquaMonitor puts OF-1 5 km west,
# in 76 m of water, where its casts reach 440 m; the table's position is 462 m deep.
const OSLOFJORD_CTD_STATIONS = [
    Station(name = "D-2",  longitude = 10.4210, latitude = 59.6280),  # Inner Drammensfjord
    Station(name = "D-3",  longitude = 10.3140, latitude = 59.7060),  # Solumstrand
    Station(name = "LA-1", longitude = 10.0520, latitude = 59.0190),  # Larviksfjord
    Station(name = "MO-2", longitude = 10.6780, latitude = 59.4840),  # Kippenes
    Station(name = "OF-1", longitude = 10.7540, latitude = 59.0410),  # Torbjørnskjær
    Station(name = "OF-5", longitude = 10.4580, latitude = 59.4870),  # Breiangen
    Station(name = "S-9",  longitude = 11.1620, latitude = 59.1140),  # Haslau, Singlefjord
    Station(name = "TØ-1", longitude = 10.3550, latitude = 59.2030),  # Vestfjord
    Station(name = "Ø-1",  longitude = 10.8340, latitude = 59.1370),  # Leira, Vesterelva
]

# The programme's other stations with casts in 2024-2026, at AquaMonitor's positions. Exact.
const OSLOFJORD_PROGRAMME_STATIONS = [
    Station(name = "BO-1",  longitude = 10.536233, latitude = 59.189933),  # Bolærne,          cast to 58 m
    Station(name = "I-1",   longitude = 11.001966, latitude = 59.109432),  # Ramsø, Østerelva, cast to 51 m
    Station(name = "I-4",   longitude = 10.951380, latitude = 59.184500),  # Kallera,          cast to 11 m
    Station(name = "KF-1",  longitude = 10.729045, latitude = 59.288682),  # Krokstadfjorden,  cast to 15 m
    Station(name = "SAN-5", longitude = 10.269300, latitude = 59.549400),  # Sandebukta,       cast to 51 m
    Station(name = "SF-3",  longitude = 10.247500, latitude = 59.068000),  # Kvernberget,      cast to 83 m
    Station(name = "VAL-1", longitude = 10.525400, latitude = 59.270100),  # Vallø,            cast to 129 m
]

# NIVA's Inner Oslofjord stations, at AquaMonitor's positions. Exact.
const INNER_OSLOFJORD_STATIONS = [
    Station(name = "Dk1", longitude = 10.569384, latitude = 59.814999),  # Steilene,         cast to 99 m
    Station(name = "Ep1", longitude = 10.723783, latitude = 59.786301),  # Bunnefjorden,     cast to 149 m
    Station(name = "Im2", longitude = 10.628217, latitude = 59.622032),  # Elle,             cast to 202 m
    Station(name = "Fl1", longitude = 10.574734, latitude = 59.754066),  # Vestfjorden,      cast to 166 m
    Station(name = "Gk1", longitude = 10.543083, latitude = 59.701218),  # Gråøyrenna,       cast to 112 m
    Station(name = "Bl4", longitude = 10.568850, latitude = 59.880833),  # Bærumsbassenget,  cast to 31 m
    Station(name = "Cj1", longitude = 10.506983, latitude = 59.840483),  # Holmenfjorden,    cast to 56 m
]

# The Drammensfjord's deepest path, from the Drammenselva mouth over the Svelvik sill to Breiangen
# (README.md): distance along it and depth on the 193 m grid it was found on.
const DRAMMENSFJORD_SECTION = [
    Station(name = "DS-00", longitude = 10.21026, latitude = 59.73954),  #  0.00 km,   5.0 m
    Station(name = "DS-01", longitude = 10.21709, latitude = 59.73597),  #  0.55 km,   4.3 m
    Station(name = "DS-02", longitude = 10.22393, latitude = 59.73418),  #  1.02 km,  13.0 m
    Station(name = "DS-03", longitude = 10.23077, latitude = 59.73239),  #  1.49 km,   6.2 m
    Station(name = "DS-04", longitude = 10.23761, latitude = 59.73239),  #  1.88 km,   7.0 m
    Station(name = "DS-05", longitude = 10.24444, latitude = 59.73418),  #  2.35 km,  35.7 m
    Station(name = "DS-06", longitude = 10.25128, latitude = 59.73418),  #  2.73 km,  44.1 m
    Station(name = "DS-07", longitude = 10.25812, latitude = 59.73239),  #  3.20 km,  52.8 m
    Station(name = "DS-08", longitude = 10.26496, latitude = 59.72881),  #  3.76 km,  64.0 m
    Station(name = "DS-09", longitude = 10.27179, latitude = 59.72524),  #  4.31 km,  70.1 m
    Station(name = "DS-10", longitude = 10.27863, latitude = 59.72166),  #  4.86 km,  76.4 m
    Station(name = "DS-11", longitude = 10.28547, latitude = 59.71808),  #  5.42 km,  79.6 m
    Station(name = "DS-12", longitude = 10.29231, latitude = 59.71451),  #  5.97 km,  83.2 m
    Station(name = "DS-13", longitude = 10.29915, latitude = 59.71093),  #  6.53 km,  86.5 m
    Station(name = "DS-14", longitude = 10.30598, latitude = 59.70735),  #  7.08 km,  94.6 m
    Station(name = "DS-15", longitude = 10.31282, latitude = 59.70735),  #  7.47 km,  96.6 m
    Station(name = "DS-16", longitude = 10.31966, latitude = 59.70557),  #  7.94 km,  98.4 m
    Station(name = "DS-17", longitude = 10.32650, latitude = 59.70557),  #  8.32 km,  99.4 m
    Station(name = "DS-18", longitude = 10.33333, latitude = 59.70378),  #  8.79 km, 102.8 m
    Station(name = "DS-19", longitude = 10.34017, latitude = 59.70378),  #  9.18 km,  95.6 m
    Station(name = "DS-20", longitude = 10.34701, latitude = 59.70378),  #  9.56 km,  64.4 m
    Station(name = "DS-21", longitude = 10.35385, latitude = 59.70378),  #  9.95 km,  90.8 m
    Station(name = "DS-22", longitude = 10.36068, latitude = 59.70020),  # 10.50 km, 108.8 m
    Station(name = "DS-23", longitude = 10.36752, latitude = 59.69662),  # 11.05 km, 110.9 m
    Station(name = "DS-24", longitude = 10.37436, latitude = 59.69305),  # 11.61 km, 113.4 m
    Station(name = "DS-25", longitude = 10.38120, latitude = 59.68947),  # 12.16 km, 114.3 m
    Station(name = "DS-26", longitude = 10.38803, latitude = 59.68589),  # 12.72 km, 114.1 m
    Station(name = "DS-27", longitude = 10.38803, latitude = 59.68232),  # 13.11 km, 114.0 m
    Station(name = "DS-28", longitude = 10.39145, latitude = 59.67874),  # 13.59 km, 114.0 m
    Station(name = "DS-29", longitude = 10.39829, latitude = 59.67516),  # 14.14 km, 114.7 m
    Station(name = "DS-30", longitude = 10.39829, latitude = 59.67159),  # 14.54 km, 114.0 m
    Station(name = "DS-31", longitude = 10.40171, latitude = 59.66801),  # 15.02 km, 114.0 m
    Station(name = "DS-32", longitude = 10.40855, latitude = 59.66443),  # 15.57 km, 118.3 m
    Station(name = "DS-33", longitude = 10.41197, latitude = 59.66086),  # 16.05 km, 119.9 m
    Station(name = "DS-34", longitude = 10.41197, latitude = 59.65728),  # 16.45 km, 120.8 m
    Station(name = "DS-35", longitude = 10.41197, latitude = 59.65370),  # 16.84 km, 121.1 m
    Station(name = "DS-36", longitude = 10.41197, latitude = 59.65013),  # 17.24 km, 121.6 m
    Station(name = "DS-37", longitude = 10.41197, latitude = 59.64655),  # 17.64 km, 121.9 m
    Station(name = "DS-38", longitude = 10.41197, latitude = 59.64297),  # 18.04 km, 121.9 m
    Station(name = "DS-39", longitude = 10.41538, latitude = 59.63940),  # 18.51 km, 121.9 m
    Station(name = "DS-40", longitude = 10.41880, latitude = 59.63582),  # 18.99 km, 121.9 m
    Station(name = "DS-41", longitude = 10.42222, latitude = 59.63224),  # 19.47 km, 116.8 m
    Station(name = "DS-42", longitude = 10.42222, latitude = 59.62867),  # 19.86 km, 116.8 m
    Station(name = "DS-43", longitude = 10.42222, latitude = 59.62509),  # 20.26 km, 114.8 m
    Station(name = "DS-44", longitude = 10.41538, latitude = 59.62151),  # 20.82 km,  85.0 m
    Station(name = "DS-45", longitude = 10.41197, latitude = 59.61973),  # 21.21 km,  33.9 m
    Station(name = "DS-46", longitude = 10.41197, latitude = 59.61615),  # 21.61 km,  13.5 m
    Station(name = "DS-47", longitude = 10.41197, latitude = 59.61257),  # 22.00 km,  11.6 m
    Station(name = "DS-48", longitude = 10.41197, latitude = 59.60900),  # 22.40 km,  11.6 m
    Station(name = "DS-49", longitude = 10.41538, latitude = 59.60542),  # 22.88 km,  11.6 m
    Station(name = "DS-50", longitude = 10.41880, latitude = 59.60184),  # 23.36 km,  12.4 m
    Station(name = "DS-51", longitude = 10.42222, latitude = 59.59827),  # 23.83 km,  12.2 m
    Station(name = "DS-52", longitude = 10.42564, latitude = 59.59469),  # 24.31 km,  12.5 m
    Station(name = "DS-53", longitude = 10.42564, latitude = 59.59290),  # 24.51 km,  10.8 m, Svelvik sill
    Station(name = "DS-54", longitude = 10.42564, latitude = 59.59111),  # 24.71 km,  12.3 m
    Station(name = "DS-55", longitude = 10.42222, latitude = 59.58754),  # 25.18 km,  14.4 m
    Station(name = "DS-56", longitude = 10.42222, latitude = 59.58396),  # 25.58 km,  17.0 m
    Station(name = "DS-57", longitude = 10.42222, latitude = 59.58038),  # 25.98 km,  18.1 m
    Station(name = "DS-58", longitude = 10.42222, latitude = 59.57681),  # 26.38 km,  19.0 m
    Station(name = "DS-59", longitude = 10.42222, latitude = 59.57323),  # 26.77 km,  22.5 m
    Station(name = "DS-60", longitude = 10.41880, latitude = 59.56965),  # 27.25 km,  27.7 m
    Station(name = "DS-61", longitude = 10.41197, latitude = 59.56608),  # 27.80 km,  35.7 m
    Station(name = "DS-62", longitude = 10.41197, latitude = 59.56250),  # 28.20 km,  25.0 m
    Station(name = "DS-63", longitude = 10.41197, latitude = 59.55892),  # 28.60 km,  22.7 m
    Station(name = "DS-64", longitude = 10.41197, latitude = 59.55535),  # 29.00 km,  22.7 m
    Station(name = "DS-65", longitude = 10.41197, latitude = 59.55177),  # 29.40 km,  23.9 m
    Station(name = "DS-66", longitude = 10.41197, latitude = 59.54819),  # 29.79 km,  27.6 m
    Station(name = "DS-67", longitude = 10.40855, latitude = 59.54462),  # 30.27 km,  32.0 m
    Station(name = "DS-68", longitude = 10.40855, latitude = 59.54104),  # 30.67 km,  37.4 m
    Station(name = "DS-69", longitude = 10.40513, latitude = 59.53746),  # 31.15 km,  56.0 m
    Station(name = "DS-70", longitude = 10.40513, latitude = 59.53389),  # 31.54 km,  73.6 m
    Station(name = "DS-71", longitude = 10.40171, latitude = 59.53031),  # 32.02 km,  88.9 m
    Station(name = "DS-72", longitude = 10.40171, latitude = 59.52673),  # 32.42 km, 111.0 m
    Station(name = "DS-73", longitude = 10.40171, latitude = 59.52316),  # 32.82 km, 125.6 m
    Station(name = "DS-74", longitude = 10.40513, latitude = 59.51958),  # 33.29 km, 157.2 m
    Station(name = "DS-75", longitude = 10.40513, latitude = 59.51600),  # 33.69 km, 156.0 m
    Station(name = "DS-76", longitude = 10.40513, latitude = 59.51422),  # 33.89 km, 153.4 m
]

# IMR's ØkoKyst and Tisler stations, at Vannmiljø's positions. Exact. Torbjørnskjær, Breiangen vest,
# Bastøy and Missingene are NIVA's OF-1, OF-5, OF-4 and OF-2.
const IMR_STATIONS = [
    Station(name = "Bastøy",         longitude = 10.590569, latitude = 59.358670),  # cast to 280 m
    Station(name = "Breiangen vest", longitude = 10.458333, latitude = 59.486668),  # cast to 190 m
    Station(name = "Færder",         longitude = 10.550012, latitude = 58.966672),  # cast to 150 m
    Station(name = "Missingene",     longitude = 10.653178, latitude = 59.191948),  # cast to 350 m
    Station(name = "Tisler-1",       longitude = 10.975667, latitude = 58.992333),  # cast to 130 m
    Station(name = "Tisler-5",       longitude = 10.958000, latitude = 58.997333),  # cast to 120 m
    Station(name = "Torbjørnskjær",  longitude = 10.760812, latitude = 59.040700),  # cast to 440 m
]

const NIVA_CTD_STATIONS =
    vcat(OSLOFJORD_CTD_STATIONS, OSLOFJORD_PROGRAMME_STATIONS, INNER_OSLOFJORD_STATIONS)
const OSLOFJORD_HYDROGRAPHY_STATIONS = vcat(NIVA_CTD_STATIONS, IMR_STATIONS)

# Fixed-point temperature sites, METreport 11/2017 §3.4, placed from its prose and maps. APPROXIMATE.
const OSLOFJORD_TEMPERATURE_STATIONS = [
    Station(name = "Åsgårdstrand", longitude = 10.470, latitude = 59.322),  # Scanmar, 1 m
    Station(name = "Sjøstrand",    longitude = 10.471, latitude = 59.811),
    Station(name = "Hvalstrand",   longitude = 10.483, latitude = 59.837),
    Station(name = "Storøyodden",  longitude = 10.630, latitude = 59.894),
]

# The NIVA stations whose casts start the water below the Drøbak and Svelvik sills
# (`CastCorrectedInitialState`), at AquaMonitor's positions. Exact.
const OSLOFJORD_BASIN_CAST_STATIONS = [
    Station(name = "Im2", longitude = 10.628217, latitude = 59.622032),  # Elle, Drøbak Sound
    Station(name = "D-2", longitude = 10.420983, latitude = 59.627766),  # Drammensfjord
    Station(name = "D-3", longitude = 10.313517, latitude = 59.705883),  # Drammensfjord
    Station(name = "Hm4", longitude = 10.612483, latitude = 59.675934),  # Oscarsborg
    Station(name = "En-1", longitude = 10.553669, latitude = 59.679659), # Engene
    Station(name = "Gk1", longitude = 10.543083, latitude = 59.701218),  # Gråøyrenna
    Station(name = "Gl2", longitude = 10.579634, latitude = 59.717384),  # Håøya
    Station(name = "Fl1", longitude = 10.574734, latitude = 59.754066),  # Vestfjorden
    Station(name = "Ej1", longitude = 10.516217, latitude = 59.793167),  # VEAS
    Station(name = "Dk1", longitude = 10.569384, latitude = 59.814999),  # Steilene
    Station(name = "Dm1", longitude = 10.617300, latitude = 59.827200),  # Steilene Nord
    Station(name = "Ep1", longitude = 10.723783, latitude = 59.786301),  # Bunnefjorden
    Station(name = "Cp2", longitude = 10.701384, latitude = 59.865002),  # Oksvalflua
    Station(name = "Cq1", longitude = 10.736667, latitude = 59.872849),  # Bekkelagsbassenget
    Station(name = "Bn1", longitude = 10.646733, latitude = 59.880634),  # Lysakerfjorden
]

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

# Viker's tide over the boundary's, from `scripts/measure_boundary_tide.jl` over 2024-10-01 to
# 2026-09-30.
const VIKER_TIDE_GAINS = [
    TidalGain("M2", 0.704, 3.3),     # 17.41 cm at the boundary, 12.26 observed, ± 0.08
    TidalGain("S2", 1.377, -27.6),   # 1.72 cm at the boundary, 2.37 observed, ± 0.08
    TidalGain("N2", 2.05, 1.7),      # 1.44 cm at the boundary, 2.95 observed, ± 0.08
    TidalGain("K2", 1.315, -10.0),   # 0.52 cm at the boundary, 0.68 observed, ± 0.08
    TidalGain("O1", 2.175, 21.0),    # 1.24 cm at the boundary, 2.7 observed, ± 0.08
    TidalGain("K1", 0.244, -95.1),   # 1.11 cm at the boundary, 0.27 observed, ± 0.08
]

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


data_root = fjord_data_root("oslofjorden")
start_date = DateTime(2024, 10, 1)
# Two years, two whole winters.
stop_date = DateTime(2026, 10, 1)
# The full-field writers start a new file each month, so a finished month can leave the VM.
monthly_splits = [(start_date + Month(month) - start_date) / Second(1) for month in 1:23]
FT = Float32
# A geometric stretch from a 1 m surface cell to the deepest sounding, with 41-132 m in 5.06 m levels
# against WENO's vertical diffusion in the sill basins (README.md).
z_faces = vcat(
    [-467.0, -433.5, -400.0, -366.5, -333.0, -299.5, -266.0, -232.5, -199.0, -165.5],
    range(-132.0, -41.0; length = 19),
    [-32.0, -25.0, -19.0, -14.5, -10.8, -7.9, -5.5, -3.7, -2.2, -1.0, 0.0],
)

FjordConfig(
    grid_config = EvenGrid(
        # 160 m both ways at 59.44°N. Every coefficient below tied to Δx is scaled to it.
        size      = (424, 681, length(z_faces) - 1),
        halo      = (7, 7, 7),
        longitude = (10.00, 11.20),
        latitude  = (58.95, 59.93),
        z_faces   = z_faces,
        float_type = FT,
    ),
    # Basin and sill depths each knob gives: README.md.
    bathymetry_config = DybdedataConfig(
        data_root             = data_root,
        output_file           = "bathymetry.nc",
        plot_file             = "bathymetry.png",
        # A 40 m native grid, which holds Svelvik's 11 m.
        raw_resolution_factor = 4,
        padding_cells         = 2,
        include_contours      = false,
        contour_stride        = 10,
        interpolation_passes  = 1,
        major_basins          = 1,
        minimum_depth         = 2.0,
        # No headland within five rows of the open edge.
        open_boundary_land_cells = 5,
        # Every island kept, Kaholmen in the Drøbak Sound among them.
        max_island_cells      = 0,
        close_narrow_passages = true,
        spike_ratio           = 0.5,
        minimum_cell_fraction = 0.2,
        # 0.25 took 34 m off Fl1 and 24 m off Hm4; 0.5 takes 8 and 4.
        max_slope_factor      = 0.5,
        # Regridding took Svelvik to 5.6 m against 11 m real.
        preserve_sills        = true,
        # The Drøbak jetty, whose ~6 m openings the grid cannot hold. Without it the sill is the
        # 23.5 m channel the jetty closes rather than the 21 m east channel.
        depth_caps            = [
            DepthCap(longitude = [10.6096, 10.6100, 10.6005], latitude = [59.6725, 59.6648, 59.6641], depth = 6.0),
        ],
        # Walled-off cells in the boundary row sat below −1.5 °C.
        raise_boundary_pockets = true,
        geonorge_cache        = true,
        regrid_cache          = false,
        # Removes the water the gridding put where the FileGDB has no data (README.md).
        max_sounding_distance = 2000.0,
    ),
    # The 160 m nest where it reaches, the 800 m run elsewhere, and the sill basins from the NIVA
    # casts within 14 days of the start, Im2's only below its sill (README.md).
    initial_state_config = CastCorrectedInitialState(
        NorKystV3InitialStateConfig(
            data_root        = data_root,
            output_directory = "norkyst_v3_initial",
            output_file      = "initial_state.nc",
            nest_member      = "m71",
            parent_member    = "m00",
            parameters       = ["temperature", "salinity", "u_eastward", "v_northward"],
            start_date       = start_date,
        ),
        VannmiljoObservations(
            data_root        = data_root,
            output_directory = joinpath("observations", "vannmiljo_niva_initial"),
            start_date       = start_date - Day(14),
            stop_date        = start_date + Day(14),
            contractors      = ["NIVA"],
            stations         = OSLOFJORD_BASIN_CAST_STATIONS,
        ),
        Day(14),
        ["Im2"],
    ),
    river_config = NVERiversConfig(
        data_root  = data_root,
        output_file = "rivers_nve.nc",
        plot_file   = "rivers_nve.png",
        start_date  = start_date,
        stop_date   = stop_date,
        # Every mouth this size or larger, discovered from NVE's river network (README.md).
        minimum_discharge = 0.5,
        # A vertical Courant number of 0.18 on the 1 m top cell at the 3-minute time step.
        maximum_inflow_velocity = 1.0e-3,
        # Gauges attached to discovered mouths, keyed by `vassdragsnr`. The rest run at their
        # catchment's normal.
        outlets = [
            # Glomma's two mouths share Solbergfoss's discharge, split by assumption (README.md),
            # with Sarpfossen's temperature.
            NVERiver(
                vassdragsnr = "002.A21", name = "Glomma (Osterelva)",
                discharge_station = "2.605.0", temperature_station = "2.1087.0",
                discharge_fraction = 2 // 3,
            ),
            NVERiver(
                vassdragsnr = "002.2A", name = "Glomma (Vesterelva)",
                discharge_station = "2.605.0", temperature_station = "2.1087.0",
                discharge_fraction = 1 // 3,
            ),
            NVERiver(
                vassdragsnr = "012.A2", name = "Drammenselva",
                discharge_station = "12.534.0", temperature_station = "12.534.0",
            ),
            # Mossevassdraget. No temperature series.
            NVERiver(vassdragsnr = "003.A4", discharge_station = "3.23.0"),
            # Akerselva's temperature gauge is a sun-exposed sensor, so none is attached.
            NVERiver(
                vassdragsnr = "006.A10", name = "Akerselva",
                discharge_station = "6.38.0",
            ),
            # Lysakerelva, as Akerselva.
            NVERiver(vassdragsnr = "007.A0", discharge_station = "7.29.0"),
            NVERiver(
                vassdragsnr = "008.A2",
                discharge_station = "8.2.0", temperature_station = "8.2.0",
            ),
            # Lierelva. No temperature series.
            NVERiver(vassdragsnr = "011.A0", discharge_station = "11.6.0"),
        ],
    ),
    # Hourly, so the Flather condition sees the tide. The reference run `m00`, the one member
    # published over the whole period; its levels stop at 300 m, the edge reaches 317 m.
    boundary_data_config = FerryBoxCorrectedBoundaries(
        TideCorrectedBoundaries(
            NorKystV3BoundaryDataConfig(
                data_root        = data_root,
                output_directory = "norkyst_v3_hourly",
                output_file      = "boundaries.nc",
                plot_file        = "boundaries.png",
                # The western and eastern walls stay closed (README.md).
                open_edges       = :south,
                margin           = 0.05,
                architecture     = :auto,
                member           = "m00",
                # The download writes `ubar_eastward` and `vbar_northward` as depth means.
                parameters       = [
                    "temperature", "salinity", "u_eastward", "v_northward", "zeta",
                    "ubar_eastward", "vbar_northward",
                ],
                start_date       = start_date,
                stop_date        = stop_date,
            ),
            VIKER_TIDE_GAINS,
        );
        # NIVA's FerryBox minutes (README.md).
        directory      = joinpath(data_root, "observations", "niva_ferrybox_2024-2026"),
        # Minutes within 0.05° (5.6 km) north of the edge's cell centres.
        band           = 0.05,
        # A week, so a day with one crossing does not set the edge alone.
        smoothing_days = 7,
        # Fitted to the Færder casts (README.md).
        taper_depth    = 25.0,
    ),
    atmosphere_config = NORA3Config(
        data_root        = data_root,
        output_directory = "nora3",
        output_file      = "atmosphere.nc",
        plot_file        = "atmosphere.png",
        resolution       = 0.02,
        padding          = 0.1,
        start_date       = start_date,
        stop_date        = stop_date,
        # Snow falls into the water, which supplies the heat that melts it.
        snowfall         = true,
    ),
    model_config = CoupledHydrostaticModel(
        buoyancy           = SeawaterBuoyancy(FT, equation_of_state = TEOS10EquationOfState(FT)),
        # A sponge damps what the open boundary radiates (README.md).
        closure            = BoundarySponge(
            # `minimum_tke` is the vertical mixing below the surface layer (README.md).
            base = (
                CATKEVerticalDiffusivity(FT; minimum_tke = 7e-7),
                # The 2Δx mode damped in 72 min; κ = 0, as diffusion overshot at fronts (README.md).
                HorizontalScalarBiharmonicDiffusivity(FT; ν = 9.4e3, κ = 0),
            ),
            width_cells = 20,
            # 30 and 15 on 193 m cells, scaled by (160/193)².
            viscosity   = 20.0,
            diffusivity = 10.0,
        ),
        # Upwind at the river mouths; T's WENO falls back to upwind near walls in the top four cells
        # only; S's WENO is bounded (README.md). `e` is named, or it gets `Centered()`.
        tracer_advection   = (
            T = RiverMouthAdvection(
                interior = SurfaceLayerAdvection(
                    surface = WENO(FT; minimum_buffer_upwind_order = 1),
                    interior = WENO(FT),
                    surface_levels = 4,
                ),
                width_cells = 3,
                narrow_cells = false,
            ),
            S = RiverMouthAdvection(
                interior = WENO(FT; bounds = (0, 40)),
                width_cells = 3,
                narrow_cells = false,
            ),
            e = WENO(FT; minimum_buffer_upwind_order = 1),
        ),
        momentum_advection = WENOVectorInvariant(FT),
        tracers            = (:T, :S),
        coriolis           = HydrostaticSphericalCoriolis(FT),
        # Clamped at freezing, with no ice (README.md).
        sea_ice            = FreezingLimitedOceanTemperature(FT),
        # The blue-green band's 3 m from the fjord's 4-6 m Secchi depths (README.md).
        radiative_forcing  = TwoColorRadiationConfig(
            first_color_fraction = 0.58,
            first_decay_scale    = 0.35,
            second_decay_scale   = 3.0,
        ),
        free_surface       = SplitExplicitFreeSurfaceConfig(cfl = 0.7),
        extra_kwargs       = (;),
    ),
    # Nudged hard on inflow, radiating on outflow (Marchesiello et al., 2001).
    boundary_condition_configs = (
        AirSeaFluxes(),
        QuadraticBottomDrag(coefficient = 0.003),
        OpenLateralBoundaryFromData(
            inflow_timescale  = 3hours,
            outflow_timescale = 360days,
        ),
    ),
    run_config = RunConfig(
        results_root       = fjord_results_root("oslofjorden"),
        architecture       = :auto,
        writers = (
            # Daily, a file a month; `e` for diagnosing the mixing.
            SnapshotWriter(
                name               = :ocean,
                output_file        = "snapshots_ocean.nc",
                variables          = (:T, :S, :u, :v, :e),
                interval           = 1day,
                split_times        = monthly_splits,
                overwrite_existing = true,
            ),
            # JLD2, as η cannot go in a NetCDF writer (`FieldSnapshotWriter`).
            SurfaceSnapshotWriter(
                name               = :surface,
                output_file        = "snapshots_surface.jld2",
                variables          = (:η, :T, :S, :u, :v),
                interval           = 1hour,
                split_times        = monthly_splits,
                overwrite_existing = true,
            ),
            # Ten minutes resolves the shortest constituents for the harmonic analysis.
            FieldStationWriter(
                name               = :tidegauge,
                output_file        = "stations_tidegauge.jld2",
                variables          = (:η,),
                stations           = OSLOFJORD_TIDE_GAUGES,
                interval           = 10minutes,
                overwrite_existing = true,
            ),
            StationWriter(
                name               = :moorings,
                output_file        = "stations_moorings.nc",
                variables          = (:u, :v),
                stations           = OSLOFJORD_MOORINGS,
                interval           = 1hour,
                overwrite_existing = true,
            ),
            # `e` for the vertical mixing.
            StationWriter(
                name               = :hydrography,
                output_file        = "stations_hydrography.nc",
                variables          = (:T, :S, :e),
                stations           = OSLOFJORD_HYDROGRAPHY_STATIONS,
                interval           = 1hour,
                overwrite_existing = true,
            ),
            SectionWriter(
                name               = :section,
                output_file        = "stations_section.nc",
                variables          = (:T, :S, :e, :u, :v),
                stations           = DRAMMENSFJORD_SECTION,
                interval           = 1hour,
                overwrite_existing = true,
            ),
            # The whole column, for the surface layer's diurnal cycle.
            StationWriter(
                name               = :surface_temperature,
                output_file        = "stations_temperature.nc",
                variables          = (:T,),
                stations           = OSLOFJORD_TEMPERATURE_STATIONS,
                interval           = 1hour,
                overwrite_existing = true,
            ),
            CheckpointWriter(interval = 12hours, cleanup = true),
        ),
        # The downslope flow takes dense water over the Drøbak sill down the slope to its own level
        # rather than letting it convect (README.md).
        callbacks = (
            ProgressCallback(name = :progress, interval = 1hour, report = progress),
            DownslopeFlow(name = :downslope, tracers = (:T, :S), friction_timescale = 1000),
        ),
        time_stepping = AdaptiveTimeStep(
            initial_time_step    = 1second,
            cfl                  = 0.3,
            max_time_step        = 3minutes,
            max_time_step_change = 1.01,
        ),
        initial_conditions = FromInitialState(),
        start_date         = start_date,
        stop_date          = stop_date,
        loops              = 1,
        pickup             = false,
    ),
    # One source per contractor, so each cast is scored once. Vannmiljø places Ø-1 and Cj1 about
    # 600 m from the positions above, hence 1 km.
    observation_configs = (
        VannmiljoObservations(
            data_root        = data_root,
            output_directory = joinpath("observations", "vannmiljo_niva"),
            start_date       = start_date,
            stop_date        = stop_date,
            contractors      = ["NIVA"],
            stations         = NIVA_CTD_STATIONS,
            match_distance   = 1000.0,
            source           = "NIVA",
        ),
        VannmiljoObservations(
            data_root        = data_root,
            output_directory = joinpath("observations", "vannmiljo_imr"),
            start_date       = start_date,
            stop_date        = stop_date,
            contractors      = ["Havforskningsinstituttet"],
            stations         = IMR_STATIONS,
            source           = "IMR",
        ),
    ),
    # Below -1 °C is colder than the freezing clamp allows; above 36 g/kg saltier than any inflow.
    validation_config = ValidationConfig(
        spin_up          = Day(31),
        depth_bands      = [(0, 5), (5, 20), (20, 50), (50, Inf)],
        station_groups   = [
            "Outer fjord" => [
                s.name for s in vcat(OSLOFJORD_CTD_STATIONS, OSLOFJORD_PROGRAMME_STATIONS, IMR_STATIONS)
                if !(s.name in ("D-2", "D-3"))
            ],
            "Drammensfjord" => ["D-2", "D-3"],
            "Inner Oslofjord" => [s.name for s in INNER_OSLOFJORD_STATIONS],
        ],
        search_radius    = 3,
        plausible_ranges = Dict("T" => (-1.0, 25.0), "S" => (0.0, 36.0)),
        stratification   = ("S", 1.0, 25.0),
    ),
)
