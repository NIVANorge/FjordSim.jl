# The Oslofjord setup: Geonorge Sjøkart Dybdedata bathymetry, an operational Norkyst-v3 initial
# state (its 160 m Oslofjord nest, filled out by the 800 m reference run), NVE rivers, hourly
# operational Norkyst-v3 data on the open southern edge with its tide corrected to Viker, NORA3
# atmosphere with snowfall and a slab sea ice, on a 424 x 681 x 39 grid of 160 m cells covering
# 10.00-11.20°E, 58.95-59.93°N. It runs from 1 October 2024 to 1 October 2026, two winters, and
# writes what `validate_simulation` scores it with: tide gauges, CTD stations, current moorings, a
# section along the Drammensfjord and fixed-point temperature sites.
#
# The box is the one MET Norway's ROMS model FjordOs CL covers (Røed et al., *A high-resolution,
# curvilinear ROMS model for the Oslofjord*, METreport 4/2016, evaluated in Hjelmervik et al.,
# *Evaluation of the FjordOs-model*, METreport 11/2017), and most station positions below come from
# those reports. Three extensions beyond the inner fjord each bring something in:
#
# - *East to 11.20°E*: the Hvaler archipelago, both arms of Glomma's delta and CTD station S-9.
#   Glomma is the largest river in the domain, and its outlet being in the *outer* fjord rather than
#   at the head is what METreport 4/2016 §1.1 calls out as making the Oslofjord's estuarine
#   circulation "deviate considerably from a classical textbook example".
# - *West to 10.00°E*: Numedalslågen, 6514 km², the third largest catchment in the FjordOs river
#   table, and CTD station LA-1.
# - *South to 58.95°N*: LA-1, at 59.019°N, about 7.7 km inside the southern edge rather than inside
#   its 3.2 km `BoundarySponge`.
#
# The sea ice is there because the surface freezes. On this box's 193 m run of 2014-15, a mild
# winter, with the ocean temperature only clamped at freezing, surface cells sat at freezing on 115
# of 641 days, up to 5.6 % of the wet surface (2014-12-27): Drammensfjorden above all, then the
# Glomma estuary and the inner fjord off Oslo, all brackish at 1-10 g/kg.
#
# The physics was measured against the NIVA casts of 2014-15 on this box at 193 m, in the runs
# archived as `oslofjorden_validation_2026-10-16` and `oslofjorden_validation_2026-10-21`: 5.06 m
# levels between 41 and 132 m, CATKE's `minimum_tke` at 7e-7, no biharmonic tracer diffusion and
# bounded WENO for salinity. Each is argued beside its value.
#
# The rivers are **discovered**, not stated: `minimum_discharge = 0.5` has `add_rivers` read NVE's
# ELVIS river network and REGINE catchments for this domain and place every mouth at or above that
# size, so no river coordinate appears in this file. On the smaller 10.2-11.02°E, 59.0-59.93°N box
# that was 21 mouths carrying 99.6 % of its 1179 m³/s; this one adds Numedalslågen in the west and
# Glomma's Singlefjord side in the east. Every discovered mouth enters whether or not it appears
# below: an outlet discovery finds already gets a catchment-normal discharge on its own, with or
# without an override. The `outlets` list is eight overrides, keyed by `vassdragsnr`, and their job
# is attaching a real HydAPI gauge — `discharge_station`, sometimes `temperature_station` — to a
# mouth discovery already found, so its inflow follows an observed series instead of its catchment
# normal; splitting a discharge with `discharge_fraction` is the exception Glomma's two mouths need,
# not what the list is for. Every other discovered mouth runs at its catchment's constant normal,
# carrying the temperature of the water it enters. Only the gauge half needs a credential, an NVE
# HydAPI key in `NVE_API_KEY` (free at https://hydapi.nve.no/Users); the map services are open, and
# nothing else here needs one either.
#
# The river config gets the same `data_root` as the rest of the setup, so the cached NVE responses
# land alongside it rather than being shared from elsewhere. It writes `rivers_nve.nc`, which
# `build_simulation` reads as a volume inflow.
#
# It also names a model, its boundary conditions and a `run_config`, so `run_simulation`
# works once the preparation steps have run. Nothing about the run has a default, so every knob is
# stated here and nowhere else — split by what dispatches on it: `CoupledHydrostaticModel` is what
# the model is, the `boundary_condition_configs` tuple is what bounds it, the `writers` tuple is what
# it writes, and `AdaptiveTimeStep` is how its clock advances. Results go to `~/FjordSim_results/oslofjorden/`,
# separate from the input data root, and each writer's file carries the launch tag so runs do not
# overwrite each other.
#
# Every dataset states its own `start_date` and `stop_date`, the period it is downloaded and prepared
# for, and the run states its own; here they are all `start_date` and `stop_date` below. A run inside
# the datasets' period needs nothing re-prepared. Changing a dataset's period means re-running its
# download and prepare steps (`add_rivers` for the rivers), and the initial state is the hourly
# instant of its config's `start_date`, which `FromInitialState` requires to be the run's.
#
# # What it can and cannot be expected to show
#
# At 160 m this grid is finer than NorKyst-800 everywhere and much finer in the vertical, but about
# twice as coarse as FjordOs CL in exactly the places FjordOs was built for. The Drøbak Sound is
# 1-2 km wide, six to twelve cells here against fifteen to twenty-five there; the Drøbak Jetty's two
# openings are ~6 m deep and tens of metres wide and are not representable at all; Svelvik is 180 m
# wide and 11 m deep, one cell. So the sill-controlled tidal jets METreport 4/2016 features, and the
# Drammensfjord exchange through Svelvik, are out of reach on this grid whatever is tuned.
#
# # Scoring it
#
# `validate_simulation` fetches what it can, scores the run's station files against it and writes
# tables and figures to `<results_root>/validation/`. Only **Kartverket water level** is public, and
# `default_observations` builds it automatically from the stations the `:tidegauge` writer records
# `η` at. The NIVA CTD series, the Fagrådet Inner Oslofjord programme, the current moorings and the
# temperature sites are held by their owners rather than published, so each wants a reader
# subtyping `AbstractObservationConfig` (or a `CsvObservations` over an export) passed to
# `validate_simulation` explicitly — see `docs/adding-a-source.md`. Their station names must be spelt
# as below, since stations are matched between a source and the run by `name`.
#
# # Preparing and running it
#
#   julia --project -m FjordSim prepare_bathymetry     --config examples/oslofjorden/oslofjorden.jl
#   julia --project -m FjordSim download_initial_state --config examples/oslofjorden/oslofjorden.jl
#   julia --project -m FjordSim prepare_initial_state  --config examples/oslofjorden/oslofjorden.jl
#   julia --project -m FjordSim add_rivers             --config examples/oslofjorden/oslofjorden.jl
#   julia --project -m FjordSim download_boundary_data --config examples/oslofjorden/oslofjorden.jl
#   julia --project scripts/measure_boundary_tide.jl examples/oslofjorden/oslofjorden.jl \
#       Viker 10.949769 59.036046   # paste its gains into `VIKER_TIDE_GAINS`
#   julia --project -m FjordSim prepare_boundary_data  --config examples/oslofjorden/oslofjorden.jl
#   julia --project -m FjordSim download_atmosphere    --config examples/oslofjorden/oslofjorden.jl
#   julia --project -m FjordSim prepare_atmosphere     --config examples/oslofjorden/oslofjorden.jl
#   julia --project -m FjordSim run_simulation         --config examples/oslofjorden/oslofjorden.jl   # needs a GPU
#   julia --project -m FjordSim validate_simulation    --config examples/oslofjorden/oslofjorden.jl

using FjordSim
using Oceananigans:
    CATKEVerticalDiffusivity,
    HydrostaticSphericalCoriolis,
    SeawaterBuoyancy,
    WENO,
    WENOVectorInvariant
using Oceananigans.TurbulenceClosures: HorizontalScalarBiharmonicDiffusivity
using Oceananigans.Units: day, days, hour, hours, minutes, second
using Dates: DateTime, Dates, Month, Second
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using FjordSim: BoundaryData, Grids, Plotting
using FjordSim.Configs: AbstractBoundaryDataConfig
using FjordSim.Validation: TIDAL_CONSTITUENTS
using LinearAlgebra: qr
using NCDatasets: NCDataset

# --- Station positions ------------------------------------------------------------------------
#
# Grouped by the observation programme they come from. Provenance matters here and is recorded per
# group: a coordinate that came from a published table is exact, and one that came from a description
# in prose is not. `StationWriter` logs how far each station moved when it was snapped to the grid,
# which is the number that says whether an approximate position mattered.

# The three permanent Norwegian Mapping Authority gauges in the fjord, from Kartverket's own station
# register (`vannstand.kartverket.no/tideapi.php?tide_request=stationlist&type=perm`). Exact.
# METreport 11/2017 Table 3 compares simulated and observed tidal constituents at all three, and
# Viker is the station FjordOs' own tidal forcing was calibrated against.
const OSLOFJORD_TIDE_GAUGES = [
    Station(name = "Viker",      longitude = 10.949769, latitude = 59.036046),
    Station(name = "Oscarsborg", longitude = 10.604861, latitude = 59.678073),
    Station(name = "Oslo",       longitude = 10.734510, latitude = 59.908559),
]

# The six Statnett/NIVA/UiO ADCP moorings deployed mid-September to late November 2014, from
# METreport 11/2017 Table 1, kept because a station writer costs next to nothing and a station with
# no observation in the window is a blank row in the scores. Exact, to six decimals. Km1/Kn2 form
# the Filtvedt-Brenntangen transect at the entrance to the Drøbak Sound; Ri1/Rl1/Rm1/Rn1 form the
# Småskjær-Evje transect across the fjord further south. The depths in the table are Statnett's
# terrain model, not the model's. METreport 11/2017 Fig. 11 notes that "the true position of
# Station Km1 is a little to the west" of where the model results were extracted, so its snapped
# offset is worth reading.
const OSLOFJORD_MOORINGS = [
    Station(name = "Km1", longitude = 10.627372, latitude = 59.582064),  # Filtvedt,       153 m
    Station(name = "Kn2", longitude = 10.646087, latitude = 59.581803),  # Brenntangen,     54 m
    Station(name = "Ri1", longitude = 10.497661, latitude = 59.350124),  # Småskjær,        20 m
    Station(name = "Rl1", longitude = 10.581023, latitude = 59.343452),  # Laksetrappa,     75 m
    Station(name = "Rm1", longitude = 10.626822, latitude = 59.352375),  # Botnegrunnen,    96 m
    Station(name = "Rn1", longitude = 10.653576, latitude = 59.363182),  # Evje,            64 m
    # The ExxonMobil bottom-mounted ADCP at the Slagen refinery, measuring at 2.5 m and 10 m since
    # 1997 and the longest current record in the fjord. METreport 11/2017 places it only in prose —
    # "northwest of Turning Dolphin at the Slagen Refinery" (§3.2.2, Fig. 3) — and tabulates no
    # position, so this one is APPROXIMATE and must be replaced with the instrument metadata when
    # the record itself is obtained. It is the comparison most likely to succeed, being in open
    # water rather than in a sill or a narrow sound.
    Station(name = "Slagen", longitude = 10.517, latitude = 59.383),
]

# The ten NIVA CTD stations of the Ytre Oslofjord monitoring programme, from METreport 11/2017
# Table 2. Exact, to four decimals. OF-1, OF-5 and D-3 are the three the report draws Hovmøller
# diagrams for; OF-1, LA-1 and D-2 the three it overlays profiles at. The other nine agree with
# AquaMonitor (`observations/niva_ctd/stations.csv`) to 200 m. OF-1 does not: AquaMonitor puts it
# 5 km west, at 10.6652°E 59.036°N, but Kartverket's soundings there are 76 m within 500 m, and every
# OF-1 cast reaches 440 m. The table's position has a 462 m sounding 228 m away, so it is kept.
const OSLOFJORD_CTD_STATIONS = [
    Station(name = "D-2",  longitude = 10.4210, latitude = 59.6280),  # Inner Drammensfjord
    Station(name = "D-3",  longitude = 10.3140, latitude = 59.7060),  # Solumstrand
    Station(name = "LA-1", longitude = 10.0520, latitude = 59.0190),  # Larviksfjord
    Station(name = "MO-2", longitude = 10.6780, latitude = 59.4840),  # Kippenes
    Station(name = "OF-1", longitude = 10.7540, latitude = 59.0410),  # Torbjørnskjær
    Station(name = "OF-5", longitude = 10.4580, latitude = 59.4870),  # Breiangen
    Station(name = "S-9",  longitude = 11.1620, latitude = 59.1140),  # Haslau, Singlefjord
    Station(name = "SF-1", longitude = 10.2460, latitude = 59.0770),  # Sandefjord
    Station(name = "TØ-1", longitude = 10.3550, latitude = 59.2030),  # Vestfjord
    Station(name = "Ø-1",  longitude = 10.8340, latitude = 59.1370),  # Leira, Vesterelva
]

# The two deep basins of the Inner Oslofjord, monitored for hydrography and oxygen by NIVA on behalf
# of Fagrådet since 1973. Not in either METreport; the positions are the programme's own, from
# AquaMonitor (`observations/niva_ctd/stations.csv`). Exact. They are here because deep-water
# renewal in these two basins is what this domain is judged on, and neither report covers it — this
# is the part of the validation that goes beyond reproducing FjordOs.
const INNER_OSLOFJORD_STATIONS = [
    Station(name = "Dk1", longitude = 10.569384, latitude = 59.814999),  # Steilene,     cast to 99 m
    Station(name = "Ep1", longitude = 10.723783, latitude = 59.786301),  # Bunnefjorden, cast to 149 m
]

# A section along the Drammensfjord's deepest water, from the Drammenselva mouth (the river's own
# cell in `rivers_nve.nc`) through the basin, over the Svelvik sill and down Sandebukta to the head of
# Breiangen. Not an observation programme: it is there to show the basin's vertical mixing and its
# renewal over the sill, the deep freshening `TODO.md` calls the main open problem. The
# positions are cell centres of the deepest path through the 193 m grid's 39-level `bathymetry.nc`,
# found by Dijkstra over wet cells with a step cost of distance / depth², every second cell of its
# 150 plus the last and the sill's shallowest. On the 160 m grid each snaps to the nearest cell, and
# a change to the bathymetry moves the path. The comments give distance along the path and the 193 m model's
# depth.
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

# Fixed-point temperature sites, METreport 11/2017 §3.4 and Figs. 3-4. All APPROXIMATE: the Scanmar
# mooring is placed only as "three kilometres south of Åsgårdstrand" and the three beaches only on a
# map. The beaches measure 40 cm below the surface in a few metres of water, so they snap to the
# model's top cell in whatever column the coastline gives them, and the report already cautions that
# "the model is therefore not expected to capture such detailed effects" — they are a sanity check on
# the surface heat budget, not a skill target.
const OSLOFJORD_TEMPERATURE_STATIONS = [
    Station(name = "Åsgårdstrand", longitude = 10.470, latitude = 59.322),  # Scanmar, 1 m
    Station(name = "Sjøstrand",    longitude = 10.471, latitude = 59.811),
    Station(name = "Hvalstrand",   longitude = 10.483, latitude = 59.837),
    Station(name = "Storøyodden",  longitude = 10.630, latitude = 59.894),
]

# --- The boundary tide ------------------------------------------------------------------------
#
# Norkyst-v3's hindcast sea level along the southern boundary has the wrong tide, and FjordSim
# carries it straight in. Fitted from 2014-05-01 to 2015-12-31, its M2 is 17.7 cm and its S2 1.4 cm,
# against 12.3 cm and 2.9 cm observed at Viker. The model's own tide at Viker is the boundary's (M2
# 17.3 cm), and its amplification from Viker to Oscarsborg is right (M2 x1.17 against x1.16
# observed). The model/observed amplitude ratio is the same at all three gauges, and the boundary
# tide is nearly uniform along the 60 km edge (M2 17.2-17.8 cm, phase within 1.3°). So one
# correction per constituent, applied along the whole edge, fixes the input. FjordOs did the same:
# METreport 4/2016 §4.3 forced its boundary with a tide corrected to Viker.
#
# `TideCorrectedBoundaries` does it after the regrid, inside `prepare_boundary_data`. It holds the
# Norkyst-v3 boundary config and forwards everything else to it. Specific to this dataset and this
# box, which is why it lives here and not in `src/`.
#
# The gains are the operational Norkyst-v3's, which has its own tide. They are measured at the
# boundary rather than at the model's Viker, as the hindcast's were, so they need no run: on the
# hindcast run the model's tide at Viker was the boundary's (M2 17.3 cm against 17.2-17.8 cm along
# the edge). `scripts/measure_boundary_tide.jl` fits the downloaded boundary sea level and Viker's
# record over the same window and prints them; `prepare_boundary_data` refuses to run until they
# are stated.

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

# Observed at Viker over the operational Norkyst-v3's boundary, each phase shift the observed phase
# minus the boundary's, measured with `scripts/measure_boundary_tide.jl` over 2024-10-01 to
# 2026-09-30. For comparison, the hindcast's over 2014-05-01 to 2015-12-31,
# measured against the model at Viker, were M2 x0.713 +2.7°, S2 x2.071 -34.0°, N2 x0.903 -4.1°,
# K2 x1.254 +32.7°, O1 x1.552 +7.7°, K1 x0.353 +17.8°. On that run Oscarsborg and Oslo gave the same
# gains within 0.02 for M2 and 0.15 for S2, and their M2 phase was 11° ahead against Viker's 3°, a
# fault of the interior that a boundary correction cannot reach.
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
        "Viker 10.949769 59.036046` and paste its output into `VIKER_TIDE_GAINS`.",
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


data_root = fjord_data_root("oslofjorden")
start_date = DateTime(2024, 10, 1)
# Two years, so each winter falls whole inside the run: 730 days, as neither 2025 nor 2026 is a leap
# year. Every dataset below is downloaded and prepared for the same period.
stop_date = DateTime(2026, 10, 1)
# The full-field writers start a new file on the 1st of each month, so a finished month can leave
# the VM.
monthly_splits = [(start_date + Month(month) - start_date) / Second(1) for month in 1:23]
# Stated here and passed to the grid and to every model component below, rather than set on
# `Oceananigans.defaults.FloatType` — see `EvenGrid`'s `float_type`.
FT = Float32
# A geometric stretch, ratio 1.25 capped at 33.5 m, from a 1 m surface cell to a deepest face at
# -467 m, with 41-132 m split into eighteen 5.06 m levels. The deepest sounding is 470.0 m, in the
# Hvalerdjupet, which METreport 4/2016 §1.1 describes as "a 400 m deep basin extending northeastward
# from the Skagerrak". The rule above would add a layer to -500.5 m, and measured on the 193 m grid,
# that layer held no water: only eight columns of the trench reached below -467 m, by at most 2.3 m,
# and `snap_partial_bottom_cells` raises such slivers of a 33.5 m layer to the face above, so the
# model is the same either way. `limit_bottom_slope` only ever makes the deepest column shallower;
# a change to the gridding or the box can move it, so re-check the deepest cell in `bathymetry.nc`
# after one.
#
# Both errors cost something. Too shallow and the basin is silently truncated —
# `snap_partial_bottom_cells` skips a sounding below the deepest face and `PartialCellBottom` clips
# it, so nothing complains. Too deep and `grid.Lz` feeds `sqrt(g·Lz)` in `SplitExplicitFreeSurface`,
# buying a shorter barotropic substep for water that is not there.
#
# The stretch alone has five 11-27 m levels between 41 and 132 m, and on them WENO's implicit
# vertical diffusion, ~|w|Δz with the basins' 1-4e-4 m/s grid-scale w, mixed the sill basins at
# 5-10x CATKE's κ. Over the first 60 days on this box at 193 m the split cut the
# Drammensfjord salt-budget diffusivity below 41-105 m from 4.5e-5-3.8e-4 to 1.7e-5-7.0e-5 m²/s, and
# the Inner Oslofjord's 1.4-5x below 41 m, at no cost in wet columns.
z_faces = vcat(
    [-467.0, -433.5, -400.0, -366.5, -333.0, -299.5, -266.0, -232.5, -199.0, -165.5],
    range(-132.0, -41.0; length = 19),
    [-32.0, -25.0, -19.0, -14.5, -10.8, -7.9, -5.5, -3.7, -2.2, -1.0, 0.0],
)

# FjordConfig holds every component side by side. Each step (prepare_forcing(), add_rivers(),
# run_simulation(), ...) picks out the components it reads and hands them to a function that
# dispatches on them, which calls the finer-grained hooks the configs below overload.
FjordConfig(
    # EvenGrid overloads the two grid hooks — domain_grid(config, architecture), which every
    # prepare_* pipeline calls, and simulation_grid(config, bathymetry_file, architecture),
    # which build_simulation reads the model grid through — both dispatching on this config's
    # type. Implemented here by wrapping LatitudeLongitudeGrid / ImmersedBoundaryGrid.
    grid_config = EvenGrid(
        # 160 m both ways: 1.20° is 67.8 km at 59.44°N (157.7-162.3 m from the northern to the
        # southern edge) and 0.98° is 109.0 km, 11.3 M cells. Every coefficient below tied to Δx is
        # scaled to it.
        size      = (424, 681, length(z_faces) - 1),
        halo      = (7, 7, 7),
        longitude = (10.00, 11.20),
        latitude  = (58.95, 59.93),
        # See beside `z_faces` above. `z_faces` is written into `bathymetry.nc` and
        # `simulation_grid` reads the file, so a change here needs `prepare_bathymetry` re-run.
        z_faces   = z_faces,
        float_type = FT,
    ),
    # The basin depths these knobs give, measured on this box at 193 m: the deepest cell within
    # 2 km of each station against the deepest sounding there, in m.
    #
    #                                          Ep1   Fl1   Im2   Hm4   Drøbak   Svelvik
    #   soundings                              154   166   211   119   211      11
    #   inverse-distance gridding, factor 2    114   138   176   104   176       8.8
    #   ... regridded                          111   127   170    88   170       6.6
    #   ... slope limit 0.25 and snap          111   114   166    71   166       5.5
    #   linear, factor 4, 0.5, preserve_sills  153   153   209    99   209      10.8
    #
    # Most of the loss was the gridding, before any smoothing stage ran; the rest was the regrid
    # averaging channels with their flanks, which `preserve_sills` undoes, and the slope limiter. The
    # domain's mean depth went from 62.4 m to 72.9 m. Drøbak's sill stays at 23.5 m against a real
    # ~20 m: no stage that only deepens can correct a sill the regrid made too deep.
    #
    # DybdedataConfig overloads bathymetry_dataset (required), plus regrid_options and
    # smoothing_options (both optional) — hooks prepare_bathymetry calls, dispatching on this
    # config's type.
    bathymetry_config = DybdedataConfig(
        data_root             = data_root,
        output_file           = "bathymetry.nc",
        plot_file             = "bathymetry.png",
        # 4, so the native grid is 40 m. The Svelvik strait is 180 m wide and 11 m deep; a
        # ~97 m native grid holds it at 8.5 m, a ~48 m one at 11.3 m, which `preserve_sills`
        # then carries onto the model grid. It was 2 while the gridding took three hours.
        raw_resolution_factor = 4,
        padding_cells         = 2,
        include_contours      = false,
        contour_stride        = 10,
        interpolation_passes  = 1,
        major_basins          = 1,
        minimum_depth         = 2.0,
        # Flood land within five rows of the open southern edge, so the boundary condition is
        # not reconciling a radiation scheme, a prescribed exterior state and a headland within
        # a couple of cells of each other. On this domain's 193 m grid the band was nearly clear
        # already and moved nine cells, none in the boundary row itself — cheap insurance rather
        # than a correction. The open edges come from `boundary_data_config` below.
        open_boundary_land_cells = 5,
        # Six and no more. On the 193 m grid every remaining interior land component was 7 cells
        # or larger with a minimum width of 2 to 4 cells — compact skerries, not the one-cell
        # ridges this stage exists to remove. Raising it would delete real topography.
        max_island_cells      = 6,
        close_narrow_passages = true,
        spike_ratio           = 0.5,
        minimum_cell_fraction = 0.2,
        # 0.5, back from 0.25. The tighter limit was measured as nearly free on the
        # inverse-distance gridding, which had already averaged the trenches flat; on the
        # linear one it takes the deep cell of every narrow trench up once per shallow
        # neighbour. On this box at 193 m 0.25 cost 34 m in Vestfjorden at Fl1 and 24 m
        # at Hm4, where 0.5 costs 8 m and 4 m. Nothing measured 0.5 as unstable: the earlier
        # runs that died did so at slivers, spikes and one-cell passages, which the stages
        # above now remove.
        max_slope_factor      = 0.5,
        # Regridding averages a channel narrower than a cell with its flanks: Svelvik came out
        # at 5.6 m against 11 m real. This restores each lost sill along one path of cells.
        preserve_sills        = true,
        # The open southern edge crosses the Hvaler trench slope. On this box at 193 m its
        # boundary row held cells walled off from the interior below some level, and the
        # domain minimum sat below −1.5 °C in them on 591 of 640 days.
        raise_boundary_pockets = true,
        geonorge_cache        = true,
        regrid_cache          = false,
        # This box reaches two places the FileGDB has no data for, and the gridding fills both
        # with water extrapolated from soundings tens of kilometres away: the land north of Drammen
        # above 59.76°N west of 10.25°E, a gap in the FileGDB's inland coverage, which became a 7 m
        # sea joined to the Drammensfjord; and Swedish territory south of Strömstad, which became a
        # 100 m one on the eastern wall. Every real wet cell lies within 1.5 km of a sounding, so
        # 2 km removes both and nothing else.
        max_sounding_distance = 2000.0,
    ),
    # NorKystV3InitialStateConfig overloads initial_state_record, initial_state_source_grid and
    # initial_state_variable_names — hooks prepare_initial_state calls — plus
    # download_initial_state, which fetches the instant `start_date` and merges the nest into the
    # 800 m run. Read by `initial_conditions = FromInitialState()` below.
    #
    # The 160 m nest reaches every sill basin and the whole fjord north of about 59.5°N, 49 % of
    # the box; the outer fjord south of about 59.3°N is the 800 m reference run. The nest's own
    # two-way parent is published only from 2024-11-23, after `start_date`. On 2024-10-01T00 the
    # two differ at the seam by a median 1.0 g/kg in surface salinity and 0.15 °C in temperature,
    # against 0.03 g/kg and 0.01 °C between neighbouring cells elsewhere.
    initial_state_config = NorKystV3InitialStateConfig(
        data_root        = data_root,
        output_directory = "norkyst_v3_initial",
        output_file      = "initial_state.nc",
        nest_member      = "m71",
        parent_member    = "m00",
        parameters       = ["temperature", "salinity", "u_eastward", "v_northward"],
        start_date       = start_date,
    ),
    # NVERiversConfig overloads river_locations, river_series, river_discharge and
    # fetch_rivers (all required, hooks add_rivers calls), all dispatching on this
    # config's type. It queries NVE, so the window follows its period and every value has
    # provenance. The HydAPI half needs an API key — `NVE_API_KEY`, free at
    # https://hydapi.nve.no/Users.
    river_config = NVERiversConfig(
        data_root  = data_root,
        output_file = "rivers_nve.nc",
        plot_file   = "rivers_nve.png",
        start_date  = start_date,
        stop_date   = stop_date,
        # Every outlet below is *discovered* from NVE's ELVIS river network and sized by
        # its REGINE catchment, so this setup states no river coordinate at all. It used
        # to state nine, copied from the OF800 dataset's outlet table, and measured
        # against NVE's own mouths they were wrong in ways nothing checked: Mosseleva 7 km
        # and Gjersjøelva 6 km from their rivers, Drammenselva 602 m from the nearest
        # water in a column at the 2 m `minimum_depth` floor, Glomma on a shelf beside
        # both of its beds rather than in either, and Glomma's whole western arm missing.
        #
        # On the smaller 10.2-11.02°E, 59.0-59.93°N box 0.5 m³/s gave 21 mouths carrying 99.6 %
        # of its 1179 m³/s, and lowering it to 0.15 would have added twelve more streams for the
        # remaining 0.3 %.
        minimum_discharge = 0.5,
        # A river's water enters through the surface of its columns, so this is a vertical
        # Courant bound on the 1 m top cell: 1 mm/s is w·Δt/Δz = 0.18 at the 3-minute
        # `max_time_step` below, leaving the time step to the horizontal CFL rather than to
        # the river mouths. `add_rivers` spreads each river over as many of the 160 m
        # columns around its mouth as its peak discharge needs to stay under it — one for a
        # small stream, tens for Glomma — and logs the count.
        maximum_inflow_velocity = 1.0e-3,
        # Overrides, keyed by the `vassdragsnr` of each mouth's terminal ELVIS segment.
        # They never add a river — a mouth discovery finds enters with or without one —
        # and an override matching no discovered mouth is an error. An override's job is
        # attaching a gauge (`discharge_station`, sometimes `temperature_station`) to a
        # mouth that would otherwise run on its catchment normal alone; Glomma's two also
        # split one gauge between its mouths with `discharge_fraction`. The mouths named
        # nowhere below — Numedalslågen, Aulivassdraget, Årosvassdraget, Hølenelva, Alna and
        # the other small ones — run at their catchment normal, which is the whole point of
        # discovery: a river with no gauge still enters at its own size.
        outlets = [
            # Glomma, Norway's largest catchment, reaching the sea through **two** mouths
            # at Fredrikstad: Østerelva east of Kråkerøy and Vesterelva west of it. Both
            # beds are resolved on this grid. ELVIS files
            # Vesterelva under the small Seutelva catchment, so discovery finds it as a
            # 1.7 m³/s stream; these two overrides put Solbergfoss's observed series
            # across both instead, and the 2/3–1/3 split is a stated assumption — NVE's
            # own Nedre Glomma flood report describes the division and publishes no
            # fraction, and none was found elsewhere.
            #
            # Discharge at Solbergfoss and temperature 40 km downstream at Sarpfossen,
            # because no station carries both.
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
            # Drammenselva, the second largest input. Mjøndalen bru carries discharge and
            # temperature on one id at 4 m a.s.l., ~5 km above the fjord.
            NVERiver(
                vassdragsnr = "012.A2", name = "Drammenselva",
                discharge_station = "12.534.0", temperature_station = "12.534.0",
            ),
            # Mossevassdraget at Moss dam, 693 km². No water temperature series.
            NVERiver(vassdragsnr = "003.A4", discharge_station = "3.23.0"),
            # Nordmarkvassdraget is Akerselva's catchment. Temperature deliberately
            # omitted: the 2020 series averages 21.2 °C and peaks at 31.1 °C, which is a
            # dry or sun-exposed sensor, not an Oslo river. The river carries the mean
            # temperature of the gauged ones instead.
            NVERiver(
                vassdragsnr = "006.A10", name = "Akerselva",
                discharge_station = "6.38.0",
            ),
            # Lysakerelva. Temperature omitted for the same reason and worse: mean
            # 19.6 °C, max 37.2 °C in 2020.
            NVERiver(vassdragsnr = "007.A0", discharge_station = "7.29.0"),
            # Sandvikselva, 226 km², and the only small stream here whose temperature
            # survives inspection: mean 8.4 °C, range 1.0-20.9 °C in 2020.
            NVERiver(
                vassdragsnr = "008.A2",
                discharge_station = "8.2.0", temperature_station = "8.2.0",
            ),
            # Lierelva at Oppsal, 223 km² of a 310 km² catchment. No temperature series.
            NVERiver(vassdragsnr = "011.A0", discharge_station = "11.6.0"),
        ],
    ),
    # The exterior state along the open southern edge, from the *hourly* operational Norkyst-v3: a
    # Flather boundary compares the model's own η against the exterior one, and daily means have
    # the tide averaged out of them. Only the boundary row is prepared at that cadence, so the file
    # is a few hundred MB. Its tide is corrected to Viker; see "The boundary tide" above.
    #
    # It is what states the open edge — for the boundary steps and for the boundary condition.
    # NorKystV3BoundaryDataConfig overloads boundary_variable_names and download_boundary_data
    # and shares boundary_time_steps and boundary_source_grid with the other NorKyst configs —
    # hooks prepare_boundary_data and download_boundary_data call.
    #
    # The reference run `m00`, the one member published over the whole run. The source's 15 levels
    # stop at 300 m, so below that the boundary repeats the 300 m value; the edge reaches 317 m.
    boundary_data_config = TideCorrectedBoundaries(
        NorKystV3BoundaryDataConfig(
            data_root        = data_root,
            output_directory = "norkyst_v3_hourly",
            output_file      = "boundaries.nc",
            plot_file        = "boundaries.png",
            # Measured from the 193 m grid's land mask. The southern wall is 311 wet cells of 351,
            # to 317 m; the northern one is dry. The western one is wet for 14 cells, 2.8 km at
            # 58.95-58.97°N off Nevlunghavn, to 41 m — all of it inside the 3.2 km
            # `BoundarySponge` of the southern edge it meets, so the westward outflow METreport
            # 4/2016 §5.1 describes turning "inside of Store Færder" leaves through the southern
            # boundary a few cells from the corner rather than being blocked. The eastern one is wet
            # for 26 cells at 59.07-59.11°N, to 66 m, where the Idefjord approach meets the Swedish
            # border, and the edge of what the Norwegian FileGDB covers: a sill fjord with no large
            # exchange, on data that thins out at the wall.
            #
            # Opening either would cost far more than it buys. `boundary_domain` takes the
            # *bounding box* of the open edges' bands, and a full-width southern band with a
            # full-height western or eastern one spans the whole domain — about 110 GB of hourly
            # download over two years against roughly 5 GB for the southern edge alone.
            open_edges       = :south,
            margin           = 0.05,
            architecture     = :auto,
            member           = "m00",
            # The source has no barotropic pair; the download writes `ubar_eastward` and
            # `vbar_northward` as the depth means of `u_eastward` and `v_northward`.
            parameters       = [
                "temperature", "salinity", "u_eastward", "v_northward", "zeta",
                "ubar_eastward", "vbar_northward",
            ],
            start_date       = start_date,
            stop_date        = stop_date,
        ),
        VIKER_TIDE_GAINS,
    ),
    # NORA3Config overloads atmosphere_time_steps, atmosphere_source_grid and
    # atmosphere_variable_names — hooks prepare_atmosphere calls — plus download_atmosphere
    # (called by the driver-level download_atmosphere) and prescribed_atmosphere /
    # prescribed_radiation (called directly by build_simulation, since the setup is
    # simulated). All dispatch on this config's type.
    atmosphere_config = NORA3Config(
        data_root        = data_root,
        output_directory = "nora3",
        output_file      = "atmosphere.nc",
        plot_file        = "atmosphere.png",
        resolution       = 0.02,
        padding          = 0.1,
        start_date       = start_date,
        stop_date        = stop_date,
        # Snow on the ice, which insulates it far better than ice does.
        snowfall         = true,
    ),
    # CoupledHydrostaticModel overloads coupled_model and model_tracers — the model hooks
    # build_simulation calls, dispatching on this config's type.
    model_config = CoupledHydrostaticModel(
        buoyancy           = SeawaterBuoyancy(FT, equation_of_state = TEOS10EquationOfState(FT)),
        # A BoundarySponge rather than a bare closure tuple: it adds a harmonic viscosity
        # and diffusivity ramping to zero over 20 cells (3.2 km) inward from whichever
        # edges `boundary_data_config` opens, and passes `base` through untouched everywhere
        # else. Without it nothing damped what the open boundary radiates — the 2020 run
        # carried velocities of 0.4 m/s std and 2 m/s peak on the boundary row against
        # 0.13 and 0.43 in the interior, and grid-scale salinity roughness 70x the interior
        # value, which accumulated for 50 days in the near-boundary bottom cells and then
        # ran away. Viscous rather than a tracer relaxation band, so it cannot fight the
        # open boundary's own nudging towards the same data.
        #
        # 20 m² s⁻¹ is set by explicit-diffusion stability, not by taste: Δt ≤ Δx²/4ν is
        # 320 s on this 160 m cell, comfortably clear of `max_time_step` below. Raising it
        # much past 35 would cap the time step instead of the CFL doing it.
        closure            = BoundarySponge(
            # `minimum_tke` is not a safety net on this grid — it *is* the vertical
            # closure over most of the column. CATKE takes `w★ = sqrt(max(minimum_tke, e))`,
            # and measured on the 2020 run at day 20.5 the prognostic `e` sat below a 7e-6
            # floor in 85 % of wet cells: 89 % below 6 m, 95-100 % below 9 m. So κ there is
            # the documented background κ = Cʰⁱᶜ·e_min/N = 0.098·e_min/N with
            # ν = 0.242·e_min/N, which the domain-median κᵤ/κᶜ of 2.38 against
            # Cʰⁱᵤ/Cʰⁱᶜ = 2.47 confirms.
            #
            # 7e-7 rather than ROMS' `GLS_KMIN` default (7.6e-6). Against the NIVA casts of
            # 9 April 2014 on this box, 7e-6 left the top 5 m too salty
            # and too weakly stratified: salinity at 25 m minus 1 m was 7.3 at Ep1 and 6.2 at
            # Dk1, against 11.0 and 9.4 observed and 9.3 and 7.4 at 7e-7. Over 640 days at 7e-7
            # the salinity RMSE against 216 casts was 1.73 g/kg and the temperature RMSE
            # 1.22 °C. The cost was overheating in 2-3 m deep cells.
            #
            # `Cᵇ` stays at Oceananigans' 0.28 rather than NumericalEarth's regional 0.01
            # (`default_ocean_closure`), deliberately: it scales the near-bottom mixing
            # length, and Stigebrandt attributed the high basin-mean diffusivity here to
            # exactly that boundary mixing.
            base = (
                CATKEVerticalDiffusivity(FT; minimum_tke = 7e-7),
                # 9.4e3 m⁴ s⁻¹. Biharmonic damping of the 2Δx mode goes as ν₄·16/Δx⁴, an
                # e-folding of 72 min on this 160 m cell — the 2e4 measured on 193 m cells,
                # scaled by (160/193)⁴. 1e5 on 193 m cells gave 14.5 min. 72 min keeps the
                # scale selectivity that leaves an 8 km baroclinic eddy alone (56 h at 8Δx).
                # For scale, Norkyst-800 — the ROMS system feeding this run — applies no
                # explicit interior viscosity at all.
                #
                # It also matters for a stability limit nothing else watches: the explicit
                # biharmonic needs Δt ≤ Δx⁴/32ν₄, which is 2180 s here, and
                # `AdaptiveTimeStep` measures only the advective CFL
                # (`cell_advection_timescale_coupled_model`).
                #
                # κ is 0. Over the first 60 days on this box, κ = 2e3
                # left the deep-basin salt budgets unchanged and the stations within 0.05 g/kg,
                # but overshot at fronts, leaving 169 wet cells below 0 g/kg, against none at 0.
                HorizontalScalarBiharmonicDiffusivity(FT; ν = 9.4e3, κ = 0),
            ),
            width_cells = 20,
            # 30 and 15 on 193 m cells, scaled by (160/193)².
            viscosity   = 20.0,
            diffusivity = 10.0,
        ),
        # T and S upwinded horizontally around the river mouths: on the Møre og Romsdal setup
        # WENO across a mouth's fresh front made water saltier than any source and a 2 m column
        # −1 °C. Three cells is WENO5's stencil, so no face outside reconstructs across a river
        # column.
        #
        # Salinity's WENO bounded: on this box plain WENO drove the
        # surface cells beside Numedalslågen's mouth below TEOS10's -32 psu floor at day 1.5.
        #
        # Oceananigans gives any tracer a NamedTuple omits the unbounded `Centered()` default,
        # so it names CATKE's `e`, and a setup adding tracers must name those too.
        tracer_advection   = (
            T = RiverMouthAdvection(interior = WENO(FT), width_cells = 3),
            S = RiverMouthAdvection(interior = WENO(FT; bounds = (0, 40)), width_cells = 3),
            e = WENO(FT),
        ),
        momentum_advection = WENOVectorInvariant(FT),
        tracers            = (:T, :S),
        coriolis           = HydrostaticSphericalCoriolis(FT),
        # A slab that grows and melts in place and never drifts, as the landfast ice of the
        # sheltered inner basins does. The ice takes no salt: it forms in water of 1-10 g/kg, and
        # a fixed ice salinity above the water's would drive the water's salinity negative.
        # 2 W m⁻¹ K⁻¹ is ClimaSeaIce's conductivity.
        sea_ice            = SlabSeaIce(ice_salinity = 0, conductivity = 2),
        # Overloads free_surface(config, grid) — its own hook, called from inside
        # coupled_model once the grid exists.
        free_surface       = SplitExplicitFreeSurfaceConfig(cfl = 0.7),
        # Anything else the three constructors coupled_model calls accept, one slot each:
        # :ocean_model, :ocean_simulation, :coupled_model.
        # Nothing extra here, but stated rather than defaulted like every other field.
        extra_kwargs       = (;),
    ),
    # Each piece overloads boundary_condition_sides, and the model config's
    # model_boundary_conditions merges them. Air-sea fluxes and quadratic bottom drag are
    # separate pieces, so either can be swapped alone, plus the genuinely open southern
    # edge — which edge, and whose exterior state, come from `boundary_data_config` above.
    # Dropping it would close the domain; tuple order is merge precedence.
    #
    # The two timescales are Marchesiello et al. (2001): nudge hard towards the data on
    # inflow, let radiation do the work on outflow.
    boundary_condition_configs = (
        AirSeaFluxes(),
        QuadraticBottomDrag(coefficient = 0.003),
        OpenLateralBoundaryFromData(
            # 3 hours, not the 1 day this was. At the time steps this run actually takes
            # (~10 s), a one-day timescale relaxes by 1.2e-4 per step against an advective
            # rate into the boundary cell of ~2.5e-3 s⁻¹ — 200x weaker, which is why the
            # nudging never caught the drift it was there to catch. Both numbers are
            # tuning knobs: too strong an inflow nudge over-constrains an open boundary
            # and reflects.
            inflow_timescale  = 3hours,
            outflow_timescale = 360days,
        ),
    ),
    # How the model is run. Writers, callbacks and time stepping each dispatch on their own
    # supertype; the run itself is read through the run config's accessors.
    run_config = RunConfig(
        results_root       = fjord_results_root("oslofjorden"),
        architecture       = :auto,
        # Both overload attach_writer! — what the run writes. `variables` may name anything
        # `Oceananigans.fields` exposes on the ocean model — add `:w`, `:η` or a
        # biogeochemical tracer by naming it here. Dropping the `CheckpointWriter` turns
        # checkpointing off entirely.
        writers = (
            # Daily. One record of five fields on this grid is 225 MB, so three-hourly over 730 days
            # would be about 1.3 TB against 164 GB daily. The full fields are for maps, sections
            # and the M2 amplitude and phase maps, none of which needs sub-daily sampling —
            # everything that does is the surface file or a station below. A new file every
            # month, ~7 GB.
            #
            # `e` is CATKE's TKE tracer, in for diagnosis: the 2026-08-25 run died at day 11.5 on a
            # TEOS10 `sqrt(Sᴬ + 32)` DomainError — one cell below -32 psu — and the file carried
            # nothing that could show what led there. `e` drives the vertical diffusivity, so it is
            # the field most likely to hold the precursor.
            SnapshotWriter(
                name               = :ocean,
                output_file        = "snapshots_ocean.nc",
                variables          = (:T, :S, :u, :v, :e),
                interval           = 1day,
                split_times        = monthly_splits,
                overwrite_existing = true,
            ),
            # The top level of T, S, u and v, and `η`, hourly, a new file every month. One record
            # is 5.8 MB on this grid, about 4.2 GB a month and 100 GB over 730 days.
            #
            # `η` in JLD2 because it cannot go in the NetCDF file: Oceananigans' NetCDF writer
            # cannot emit a `(Center, Center, Nothing)` user output at all — the free surface
            # asks for a singleton `z_aaf = [0.0]` against the grid's own `z_aaf`, and it fails
            # beside the 3D fields, alone in its own file, and with grid metrics off. See
            # `FieldSnapshotWriter`.
            SurfaceSnapshotWriter(
                name               = :surface,
                output_file        = "snapshots_surface.jld2",
                variables          = (:η, :T, :S, :u, :v),
                interval           = 1hour,
                split_times        = monthly_splits,
                overwrite_existing = true,
            ),
            # Ice and snow thickness and ice concentration, daily, a new file every month.
            SeaIceSnapshotWriter(
                name               = :sea_ice,
                output_file        = "snapshots_sea_ice.jld2",
                variables          = (:h, :ℵ, :hs),
                interval           = 1day,
                split_times        = monthly_splits,
                overwrite_existing = true,
            ),
            # Ten minutes, and the only writer at that cadence. A tide gauge comparison is a
            # harmonic analysis: METreport 11/2017 Table 3 resolves constituents down to MS4 at
            # 6.1033 hours, and the shortest of them needs several samples per period to come back
            # unaliased. Three points at 10 minutes over two years is about 1 MB.
            #
            # JLD2 because η cannot go in a NetCDF writer at all — see `FieldSnapshotWriter`.
            FieldStationWriter(
                name               = :tidegauge,
                output_file        = "stations_tidegauge.jld2",
                variables          = (:η,),
                stations           = OSLOFJORD_TIDE_GAUGES,
                interval           = 10minutes,
                overwrite_existing = true,
            ),
            # Hourly profiles of velocity at the ADCP sites. Hourly is what the observations are
            # and what METreport 11/2017's own analysis assumes: the 49-hour running mean it
            # separates the estuarine circulation with, and the tidal ellipses of its Table 4,
            # both need the tide resolved in the record they are taken from.
            StationWriter(
                name               = :moorings,
                output_file        = "stations_moorings.nc",
                variables          = (:u, :v),
                stations           = OSLOFJORD_MOORINGS,
                interval           = 1hour,
                overwrite_existing = true,
            ),
            # Hourly T/S profiles at the CTD stations and the two inner-basin stations. The
            # observations are a dozen casts a year, so hourly is far more than the comparison at
            # cast dates needs — but it is what resolves a deep-water renewal event, which happens
            # over days and is the thing this domain is judged on. `e` rides along because the
            # standing conclusions about this fjord are both about vertical mixing — too weak in
            # the open fjord, too vigorous inside the sills — and the TKE is what would show it.
            StationWriter(
                name               = :hydrography,
                output_file        = "stations_hydrography.nc",
                variables          = (:T, :S, :e),
                stations           = vcat(OSLOFJORD_CTD_STATIONS, INNER_OSLOFJORD_STATIONS),
                interval           = 1hour,
                overwrite_existing = true,
            ),
            # Hourly whole-column profiles along the Drammensfjord section: T and S for the
            # mixing and renewal, `e` for the TKE, and u, v for the exchange over the sill. One
            # file for all 77 rather than one per station, which cost about 7 minutes of start-up.
            SectionWriter(
                name               = :section,
                output_file        = "stations_section.nc",
                variables          = (:T, :S, :e, :u, :v),
                stations           = DRAMMENSFJORD_SECTION,
                interval           = 1hour,
                overwrite_existing = true,
            ),
            # The fixed-point temperature sites, hourly. The whole column is written rather than
            # the top cell alone, because the explanation METreport 11/2017 gives for its beach
            # discrepancy is a diurnal cycle in the model that is too strong, and that is a
            # statement about the surface layer rather than the surface.
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
        # Overloads attach_callback! — what the run reports while it runs. `report` is the
        # function itself, so a model whose tracers omit :T (which `progress` reads) names its
        # own here instead. An empty tuple runs silently.
        callbacks = (ProgressCallback(name = :progress, interval = 1hour, report = progress),),
        # Overloads attach_time_stepping! and initial_time_step.
        time_stepping = AdaptiveTimeStep(
            initial_time_step    = 1second,
            cfl                  = 0.3,
            max_time_step        = 3minutes,
            max_time_step_change = 1.01,
        ),
        # The Norkyst-v3 state at `start_date`, from `initial_state_config`: T, S, u and v.
        # `(T = 5.0, S = 33.0)` starts from a uniform water column at rest instead, and
        # `FromResults("snapshots_ocean_<tag>.nc")` continues from a previous run.
        initial_conditions = FromInitialState(),
        start_date         = start_date,
        stop_date          = stop_date,
        # One pass. Raise it to spin the deep basins up on the same year, carrying the
        # state over; each repetition writes its own `_loopNN` file.
        loops              = 1,
        pickup             = false,
    ),
)
