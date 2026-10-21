# The Oslofjord setup rebuilt for validation against the FjordOs reports: `oslofjorden.jl`'s physics,
# over the domain and the years those reports cover, with the in-run diagnostics a model-observation
# comparison needs.
#
# It exists to answer one question — how does FjordSim, as it actually stands, score against the same
# observations MET Norway scored their ROMS model FjordOs CL against? — so **nothing about the physics
# is changed or retuned here** beyond the two departures under "What is different". The references
# are Røed et al., *A high-resolution, curvilinear ROMS model for the Oslofjord*, METreport 4/2016,
# which documents FjordOs CL and its 1 April 2014 - 31 December 2015 hindcast, and Hjelmervik
# et al., *Evaluation of the FjordOs-model*, METreport 11/2017, which evaluates it against water
# level, currents, CTD profiles and fixed-point temperature.
#
# Nothing in `src/` is specific to this file. What it uses from the package is general: the
# Norkyst-v3 hindcast configs, the station writers, and `validate_simulation` with its Kartverket
# water-level source.
#
#   julia --project -m FjordSim prepare_bathymetry  --config examples/oslofjorden_validation.jl
#   julia --project -m FjordSim add_rivers          --config examples/oslofjorden_validation.jl
#   julia --project -m FjordSim download_boundaries --config examples/oslofjorden_validation.jl
#   julia --project -m FjordSim prepare_boundaries  --config examples/oslofjorden_validation.jl
#   julia --project -m FjordSim download_atmosphere --config examples/oslofjorden_validation.jl
#   julia --project -m FjordSim prepare_atmosphere  --config examples/oslofjorden_validation.jl
#   julia --project -m FjordSim run_simulation      --config examples/oslofjorden_validation.jl
#   julia --project -m FjordSim validate_simulation --config examples/oslofjorden_validation.jl
#
# `add_rivers` needs `NVE_API_KEY`, exactly as `oslofjorden.jl` does.
#
# # What is reused from `oslofjorden.jl`, and why it is reused rather than copied
#
# The float type, the closure (every coefficient and the `BoundarySponge`, the CATKE `minimum_tke`
# floor excepted, below), the advection schemes (salinity's excepted, below), the boundary conditions and their two Marchesiello timescales, the
# time stepping, every bathymetry smoothing knob and the whole river
# configuration are taken from the config `oslofjorden.jl` returns, not restated. A validation of a model that has drifted from
# the model being validated is worth nothing, and a copy is exactly what drifts. So a change to any of
# them is made in `oslofjorden.jl` and arrives here by itself — but a bathymetry knob is a *data*
# change, and needs this file's `prepare_bathymetry`, `add_rivers` and `prepare_boundaries`
# re-run as much as `oslofjorden.jl`'s.
#
# # What is different, and why
#
# **The dataset.** `oslofjorden.jl` reads the operational NorKyst-800m archive at
# `fou-hi/norkyst800m`, which begins 2017-02-20 and so cannot reach this window at all. This reads the
# Norkyst-v3 hindcast at `romshindcast/norkyst_v3` (`NorKystHindcastBoundariesConfig`), MET's
# continuous free run from 2012, which carries the same variables on the same grid. It also gains
# from the change: the hindcast publishes `ubar_eastward`/`vbar_northward` already rotated to geographic axes, where the operational
# collection publishes ROMS' grid-relative pair and `NorKystBoundariesConfig` has to derotate it.
#
# **The domain.** 10.00-11.20°E by 58.95-59.93°N against `oslofjorden.jl`'s 10.2-11.02 by 59.0-59.93 —
# the box FjordOs CL covers, at the same cell size. The three extensions each buy something the
# validation needs and nothing else does:
#
# - *East to 11.20°E* brings in the Hvaler archipelago, both arms of Glomma's delta and CTD station
#   S-9. Glomma is the largest river in the domain, and its outlet being in the *outer* fjord rather
#   than at the head is what METreport 4/2016 §1.1 calls out as making the Oslofjord's estuarine
#   circulation "deviate considerably from a classical textbook example". Cutting the delta in half is
#   not a small approximation.
# - *West to 10.00°E* brings in Numedalslågen — 6514 km², the third largest catchment in the FjordOs
#   river table — and CTD station LA-1.
# - *South to 58.95°N* is for LA-1 specifically. At 59.019°N it would sit about ten cells inside a
#   sixteen-cell `BoundarySponge` on the old southern edge and be unusable; here it is about
#   thirty-nine cells in. METreport 11/2017 hits the same problem from the other side, noting LA-1 is
#   "close to the southern open boundary of the FjordOs model".
#
# **Salinity advection.** Bounds-preserving WENO for S alone, because Numedalslågen is the river
# that enters among two-level columns, and plain WENO drove salinity beside its mouth below
# TEOS10's -32 psu floor at day 1.5. The measurement is at the statement, beside `model` below.
#
# **The vertical-mixing floor.** CATKE's `minimum_tke` at 7e-7 rather than `oslofjorden.jl`'s 7e-6,
# because the 9 April 2014 NIVA casts are more stratified than the run was at 7e-6. The measurement
# is beside `closure` below.
#
# **The deep levels.** Eighteen 5.06 m levels between 41 and 132 m in place of `oslofjorden.jl`'s five
# 11-27 m ones, because WENO's implicit vertical diffusion on those thick levels, not CATKE, is what
# freshened the sill basins. The measurement is beside `z_faces` below.
#
# **No horizontal tracer diffusion.** The biharmonic κ at 0, ν kept, because κ only made negative
# salinity. The measurement is beside `closure` below.
#
# **The boundary tide.** Norkyst-v3's tide is corrected to Viker, constituent by constituent, because
# its M2 is 40 % too large and its S2 half the observed. The measurement is under "The boundary
# tide" below. It is applied to the prepared file, so it takes effect only once `prepare_boundaries`
# is re-run.
#
# **The initial state.** The previous full run's state of 1 April 2015 rather than Norkyst-v3's of
# 31 March 2014, because Norkyst's deep inner basins are about 2 °C too cold and the sill basins keep
# what they start with. The measurement is beside `initial_conditions` below.
#
# **The writers.** A 640-day run cannot afford `oslofjorden.jl`'s three-hourly full-domain snapshots —
# at 7.5 M cells a record is 155 MB, so that is roughly 790 GB — and a validation does not need them:
# what it needs is high cadence at a few dozen points. So the 3D fields drop to daily (99 GB) for maps
# and sections, and five station writers carry the time resolution at about 1 GB total.
#
# **The data root.** Only the Geonorge FileGDB is shared with `oslofjorden.jl`, by absolute path,
# because the 4.5 GB national download is the same for any Norwegian domain. Everything else is a
# different period and downloads into this file's own root.
#
# # Two gates, settled by measurement
#
# Both were measured on the `bathymetry.nc` `prepare_bathymetry` writes, with `max_sounding_distance`
# in force; the reasoning is beside each value below. A change to the box or to any bathymetry knob
# reopens both.
#
# 1. **`z_faces`.** `oslofjorden.jl`'s own, to -467 m, refined between 41 and 132 m. The deepest
#    sounding is 470.0 m, but one more layer would hold no water: see the reasoning beside the value.
# 2. **`open_edges`.** `:south` alone: the only long wet wall.
#
# # What it can and cannot be expected to show
#
# Stated up front because it decides how the results should be read. At ~193 m this grid is finer
# than NorKyst-800 everywhere and much finer in the vertical, but two to three times coarser than
# FjordOs CL in exactly the places FjordOs was built for. The Drøbak Sound is 1-2 km wide, six to ten
# cells here against fifteen to twenty-five there; the Drøbak Jetty's two openings are ~6 m deep and
# tens of metres wide and are not representable at all; Svelvik is 180 m wide and 11 m deep, one
# cell. So the sill-controlled tidal jets METreport 4/2016 features, and the Drammensfjord exchange
# through Svelvik, are out of reach on this grid whatever is tuned.
#
# What should be competitive or better is everything that is not sill-controlled: tidal elevation,
# open-fjord currents (Slagen above all), the seasonal hydrography at the CTD stations, and deep-water
# renewal in the inner basins. The z-coordinate with `PartialCellBottom` removes the defect both
# reports complain about most — FjordOs' bathymetry had to be smoothed to an rx0 limit until "the
# observed slopes are everywhere steeper" than the model's (METreport 11/2017 Fig. 11), and the
# pressure-gradient error that forces it does not exist here. The only slope limit this bathymetry
# gets is `max_slope_factor = 0.5`: the deepest cell within 2 km of every inner-fjord CTD station is
# within 12 m of the deepest cast there, where the earlier bathymetry was up to 51 m shallower, and
# `preserve_sills` keeps Svelvik at 10.8 m against 11 m real. That is also why `plot_section` draws no
# smoothed-against-real bathymetry comparison, which is the main point of Fig. 11: there is little to
# compare. Its two sections are the Statnett transects — the
# Filtvedt-Brenntangen line at 59.582°N and the Småskjær-Evje line at ~59.35°N, both very nearly
# zonal, so each is one row of this grid. Establishing which half of the expectation above holds,
# quantitatively, is the point of the run.
#
# # Scoring it
#
# `validate_simulation` fetches what it can, scores the run's station files against it and writes
# tables and figures to `<results_root>/validation/`, where the skill and tidal tables are meant to be
# read beside METreport 11/2017's own (its Table 3 for the tidal constituents).
#
# Only **Kartverket water level** is public, and `default_observations` builds it automatically from
# the stations the `:tidegauge` writer records `η` at; its 10-minute, quality-flagged series has been
# checked to cover Viker, Oscarsborg and Oslo across the whole 2014-2015 window. The NIVA CTD series,
# the 2014 Statnett ADCP records, the Fagrådet Inner Oslofjord programme, the Scanmar mooring and the
# beach thermistors are all held by their owners rather than published, so each wants a reader
# subtyping `AbstractObservationConfig` (or a `CsvObservations` over an export) passed to
# `validate_simulation` explicitly — see `docs/adding-a-source.md`. Their station names must be spelt
# as below, since stations are matched between a source and the run by `name`.
#
# Drifter trajectories and the Godafoss Lagrangian statistics (METreport 11/2017 §4.5-4.6) are out of
# scope: FjordSim has no particle tracking.

using FjordSim
using FjordSim: AbstractBoundaryDataConfig, BoundaryData, Grids, Plotting
using FjordSim.Validation: TIDAL_CONSTITUENTS
using LinearAlgebra: qr
using NCDatasets: NCDataset
using Oceananigans: WENO, CATKEVerticalDiffusivity, HorizontalScalarBiharmonicDiffusivity
using Oceananigans.Units
using Dates: DateTime, Dates

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
# METreport 11/2017 Table 1. Exact, to six decimals. Km1/Kn2 form the Filtvedt-Brenntangen transect
# at the entrance to the Drøbak Sound; Ri1/Rl1/Rm1/Rn1 form the Småskjær-Evje transect across the
# fjord further south. The depths in the table are Statnett's terrain model, not the model's.
# METreport 11/2017 Fig. 11 notes that "the true position of Station Km1 is a little to the west" of
# where the model results were extracted, so its snapped offset is worth reading.
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
# renewal in these two basins is what `oslofjorden.jl`'s closure comment names as "what this domain
# is judged on", and neither report covers it — this is the part of the validation that goes beyond
# reproducing FjordOs.
const INNER_OSLOFJORD_STATIONS = [
    Station(name = "Dk1", longitude = 10.569384, latitude = 59.814999),  # Steilene,     cast to 99 m
    Station(name = "Ep1", longitude = 10.723783, latitude = 59.786301),  # Bunnefjorden, cast to 149 m
]

# A section along the Drammensfjord's deepest water, from the Drammenselva mouth (the river's own
# cell in `rivers_nve.nc`) through the basin, over the Svelvik sill and down Sandebukta to the head of
# Breiangen. Not an observation programme: it is there to show the basin's vertical mixing and its
# renewal over the sill, the deep freshening `TODO.md` calls the main open problem. The
# positions are cell centres of the deepest path through this run's 39-level `bathymetry.nc`, found
# by Dijkstra over wet cells with a step cost of distance / depth² (`thalweg.jl` in
# `gs://fjordsim/analysis/oslofjorden_validation_dz5_kappah0/`), every second cell of its 150 plus
# the last and the sill's shallowest. So they snap to themselves, and a change to the bathymetry
# moves the path. The comments give distance along the path and model depth.
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
# `TideCorrectedBoundaries` does it after the regrid, inside `prepare_boundaries`. It holds the
# Norkyst-v3 boundary config and forwards everything else to it. Specific to this dataset and this
# box, which is why it lives here and not in `src/`. The gains are measured for the hindcast's
# 2014-15 tide only. `oslofjorden.jl`'s operational archive would need its own.

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

# Observed over modelled at Viker, from the full 640-day run (stamp 20261014T142318) against
# Kartverket's 10-minute record, 2014-05-01 to 2015-12-31. Each phase shift is the observed phase
# minus the modelled. The constituents are the six above 0.4 cm at Viker that a 20-month record
# separates; standard errors are 0.1 cm. Oscarsborg and Oslo give the same gains within 0.02 for M2
# and 0.15 for S2. Their M2 phase is also 11° ahead, against Viker's 3°, which is a fault of the
# interior that a boundary correction cannot reach.
const VIKER_TIDE_GAINS = [
    TidalGain("M2", 0.713,   2.7),   # 17.27 cm modelled, 12.32 observed
    TidalGain("S2", 2.071, -34.0),   #  1.40,  2.89
    TidalGain("N2", 0.903,  -4.1),   #  3.35,  3.02
    TidalGain("K2", 1.254,  32.7),   #  0.48,  0.61
    TidalGain("O1", 1.552,   7.7),   #  1.19,  1.85
    TidalGain("K1", 0.353,  17.8),   #  1.01,  0.36
]

"""
    TideCorrectedBoundaries(source, gains)

`source`'s open-boundary data with its tide corrected per constituent by `gains`.

`prepare_boundaries` writes the file exactly as `source` would, then corrects it in place
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
BoundaryData.download_boundaries(target_grid, config::TideCorrectedBoundaries) =
    BoundaryData.download_boundaries(target_grid, config.source)

function BoundaryData.prepare_boundaries(
    target_grid,
    config::TideCorrectedBoundaries;
    coverage = nothing,
)
    result = BoundaryData.prepare_boundaries(target_grid, config.source; coverage)
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
            "Re-run `prepare_boundaries` to start from the uncorrected data.",
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

# --- The config -------------------------------------------------------------------------------

# The config `oslofjorden.jl` evaluates to, built from the file beside this one.
base = include("oslofjorden.jl")
data_root = fjord_data_root("oslofjorden_validation")
# `oslofjorden.jl`'s faces with 41-132 m split into eighteen equal levels; see beside `z_faces` below.
z_faces = let faces = base.grid_config.z_faces
    vcat(filter(z -> z <= -132, faces), collect(range(-132.0, -41.0; length = 19))[2:end-1],
         filter(z -> z >= -41, faces))
end
years = [2014, 2015]

# The three configs below are `oslofjorden.jl`'s own instances with only their location and period
# changed, rather than fresh ones restating its fields: every smoothing knob, every river override and
# threshold, and the atmosphere's resolution then carry over by construction, including any field
# added later. `include` evaluates `oslofjorden.jl` afresh, so this mutates nothing shared.

bathymetry = base.bathymetry_config
# Oslofjord's extracted copy of the national FileGDB, by absolute path — the same sharing
# `drammensfjorden.jl` does. Resolved before `data_root` moves, which it is relative to.
bathymetry.geodatabase_file = geodatabase_path(bathymetry)
bathymetry.data_root = data_root
# This box reaches two places the FileGDB has no data for, and the gridding fills both with water
# extrapolated from soundings tens of kilometres away: the land north of Drammen above 59.76°N west
# of 10.25°E, a gap in the FileGDB's inland coverage, which became a 7 m sea joined to the
# Drammensfjord; and Swedish territory south of Strömstad, which became a 100 m one on the eastern
# wall. Every real wet cell lies within 1.5 km of a sounding, so 2 km removes both and nothing else.
# `oslofjorden.jl`'s box touches neither and does not set it.
bathymetry.max_sounding_distance = 2000.0

# Discovery will find more mouths than Oslofjord's 21 because the box is larger — Numedalslågen in
# the west, and Glomma's Singlefjord side in the east — which is the intended effect. The FjordOs
# river table has 37 named rivers over a comparable domain, so the two counts should end up close.
# Every override `oslofjorden.jl` states lies well inside this box, which only grew, so none of them
# matches no mouth.
rivers = base.river_config
rivers.data_root = data_root
rivers.years = years

# Three of `oslofjorden.jl`'s gauges were chosen for 2020 and have no daily values on HydAPI in 2014
# or 2015, so those overrides are restated with gauges that do. Each was checked as daily values
# over both years.
#
# - Drammenselva temperature. Mjøndalen bru (12.534.0) carries discharge throughout but no
#   temperature at all in this window. The nearest main-stem gauge, Døvikfoss (12.298.0), stops on
#   2015-09-26, and `nve_fill_gaps!` would then hold its 11.8 °C to the end of December. Strømstøa
#   (12.15.0), on Ådalselva above Tyrifjorden, covers both years and matches Døvikfoss best of the
#   seven candidates over their 634 common days: bias -0.35 °C, RMSE 0.95 °C. Discharge stays at
#   Mjøndalen.
# - Akerselva discharge. 6.38.0 has nothing in the window; the Maridalsvatn outflow (6.9.0), 209 km²
#   of the catchment, does.
# - Lysakerelva discharge. 7.29.0 has nothing in the window and no gauge in area 7 does, so it is
#   left to `river_discharge`'s fallback, the REGINE catchment normal.
restated = Dict(
    "012.A2" => NVERiver(
        vassdragsnr = "012.A2", name = "Drammenselva",
        discharge_station = "12.534.0", temperature_station = "12.15.0",
    ),
    "006.A10" => NVERiver(vassdragsnr = "006.A10", name = "Akerselva", discharge_station = "6.9.0"),
    "007.A0" => NVERiver(vassdragsnr = "007.A0"),
)
rivers.outlets = [get(restated, outlet.vassdragsnr, outlet) for outlet in rivers.outlets]

atmosphere = base.atmosphere_config
atmosphere.data_root = data_root
atmosphere.years = years

simulation = base.simulation_config

# The first physics departure from `oslofjorden.jl`: salinity is advected by bounds-preserving WENO.
# Plain WENO is not bounded, and this domain adds Numedalslågen at Larvik, a large river entering
# among two-level columns (measured while rivers were still a relaxation to S = 0). The horizontal reconstruction across that front drove the
# neighbouring surface cells to -29 psu over a 110 psu cell beneath by day 1.5, and past -32 a step
# later, where TEOS10's `sqrt(Sᴬ + 32)` threw a DomainError on the GPU. Measured on a Larvik
# sub-domain with no atmosphere: WENO order 3 or 5, `PartialCellBottom` or `GridFittedBottom`,
# CATKE or a constant diffusivity all reproduce it, and switching the river off, upwinding only the
# horizontal, or bounding S each removes it. T and `e` keep the base's scheme. Water below 0 °C is
# real here, and a NamedTuple must name `e`, or Oceananigans gives it the unbounded `Centered()`
# default.
#
# Not changed in `oslofjorden.jl` itself, because a NamedTuple there would leave the NPZD example's
# tracers on `Centered()`.
scheme = simulation.model.tracer_advection

# The other departure: CATKE's `minimum_tke` at 7e-7, a tenfold cut in the background diffusivity
# that sets κ below about 6 m (see `oslofjorden.jl`), with the rest of the closure as it is there.
# `oslofjorden.jl` keeps 7e-6 to fit basin-mean diffusivities from density budgets. Against the
# NIVA casts of 9 April 2014 (day 8), the run with 7e-6 — rivers already a volume source — was too
# salty in the top 5 m and too weakly stratified. Salinity at 25 m minus 1 m was 7.3 at Ep1 and
# 6.2 at Dk1, against 11.0 and 9.4 observed. A 7e-7 run over the same days gave 9.3 and 7.4, with
# the 1 m salinity at Ep1 22.9 against 24.7 at 7e-6 and 20.1 observed. Over 26 days the two differed
# by about 0.1 in the 50 m salinity trend at D-2 and Ep1. See TODO.md › "Vertical mixing is too
# strong".
#
# The biharmonic tracer diffusivity κ is 0, its viscosity ν kept. Over this run's first 60 days
# (`gs://fjordsim/analysis/mixing-tests-2026-10-18/`, on the 26 levels), κ = 2e3 against 0 left the
# deep-basin salt budgets unchanged — Drammensfjord K_S below 41-105 m 4.5e-5 to 3.8e-4 in both —
# and the stations within 0.05 g/kg; what κ did do is overshoot at fronts, leaving 169 wet cells
# below 0 g/kg (minimum -0.43) at day 60, against none at κ = 0. The one systematic difference left
# is 0.04-0.10 °C colder bottom water at D-2, Ep1 and Dk1 at κ = 0, not yet explained.
sponge = simulation.model.closure
closure = BoundarySponge(
    base        = (
        CATKEVerticalDiffusivity(base.grid_config.float_type; minimum_tke = 7e-7),
        HorizontalScalarBiharmonicDiffusivity(base.grid_config.float_type; ν = sponge.base[2].ν, κ = 0),
    ),
    width_cells = sponge.width_cells,
    viscosity   = sponge.viscosity,
    diffusivity = sponge.diffusivity,
)

model = CoupledHydrostaticSimulation(
    buoyancy           = simulation.model.buoyancy,
    closure            = closure,
    tracer_advection   = (T = scheme, S = WENO(base.grid_config.float_type; bounds = (0, 40)), e = scheme),
    momentum_advection = simulation.model.momentum_advection,
    tracers            = simulation.model.tracers,
    coriolis           = simulation.model.coriolis,
    sea_ice            = simulation.model.sea_ice,
    biogeochemistry    = simulation.model.biogeochemistry,
    free_surface       = simulation.model.free_surface,
    extra_kwargs       = simulation.model.extra_kwargs,
)

FjordConfig(
    grid_config = EvenGrid(
        # 351 x 548 cells over 1.20° x 0.98° is 193.8 m in longitude at 59.4°N and 199.1 m in
        # latitude — the same cell size as `oslofjorden.jl`, deliberately. Everything tuned against Δx
        # or Δx⁴ carries over only if Δx does: the biharmonic coefficients are justified by a damping
        # rate that goes as ν₄·16/Δx⁴, and changing the resolution without re-deriving them would
        # silently invalidate the one number in that setup that was measured rather than chosen.
        # 7.50 M cells against Oslofjord's 3.24 M.
        size      = (351, 548, length(z_faces) - 1),
        halo      = base.grid_config.halo,
        longitude = (10.00, 11.20),
        latitude  = (58.95, 59.93),
        # `oslofjorden.jl`'s 26 levels, to -467 m, with the five between 41 and 132 m replaced by
        # eighteen of 5.06 m: 39 in all. Those five were 11-27 m thick, and on them WENO's implicit
        # vertical diffusion, ~|w|Δz with the basins' 1-4e-4 m/s grid-scale w, mixed the sill
        # basins at 5-10x CATKE's κ there. Over this run's first 60 days
        # (`gs://fjordsim/analysis/mixing-tests-2026-10-18/`) the refinement cut the Drammensfjord
        # salt-budget diffusivity below 41, 52, 66, 83 and 105 m from 4.5e-5, 9.7e-5, 2.0e-4, 3.3e-4
        # and 3.8e-4 m²/s to 1.7e-5, 2.0e-5, 3.2e-5, 4.3e-5 and 7.0e-5, and the Inner Oslofjord's
        # 1.4-5x below 41 m. It cost no wet columns, and moved 30 columns domain-wide by more than
        # 6 m.
        #
        # The 26 levels clear that domain's deepest sounding of
        # 463.0 m. This box extends into the Skagerrak and takes in more of the Hvalerdjupet, which
        # METreport 4/2016 §1.1 describes as "a 400 m deep basin extending northeastward from the
        # Skagerrak", and its deepest sounding is 470.0 m. `oslofjorden.jl`'s rule — a geometric
        # stretch of ratio 1.25 from a 1 m surface cell, capped at 33.5 m, with just enough 33.5 m
        # layers to clear the deepest sounding — would add a layer to -500.5 m, and measured, that
        # layer holds no water: only eight columns of the trench reach below -467 m, by at most
        # 2.3 m. Under a -500.5 m face `snap_partial_bottom_cells` raises all eight to -467 m, being
        # slivers of the 33.5 m layer; without it they stay below the face in `bathymetry.nc` and
        # `PartialCellBottom` clamps them to it. So the faces stop at -467 m, and the model is the
        # same either way.
        #
        # Both errors cost something. Too shallow and the basin is silently truncated —
        # `snap_partial_bottom_cells` skips a sounding below the deepest face and `PartialCellBottom`
        # clips it, so nothing complains. Too deep and `grid.Lz` feeds `sqrt(g·Lz)` in
        # `SplitExplicitFreeSurface`, buying a shorter barotropic substep for water that is not
        # there. `z_faces` is written into `bathymetry.nc` and `simulation_grid` reads the file, not
        # this config, so a change here needs `prepare_bathymetry` re-run.
        z_faces   = z_faces,
        # `oslofjorden.jl`'s, since the model below is built in it.
        float_type = base.grid_config.float_type,
    ),
    bathymetry_config = bathymetry,
    river_config = rivers,
    # The tide in it corrected to Viker; see "The boundary tide" above.
    boundary_config = TideCorrectedBoundaries(
        NorKystHindcastBoundariesConfig(
            data_root        = data_root,
            output_directory = "norkyst_v3_hourly",
            output_file      = base.boundary_config.output_file,
            plot_file        = base.boundary_config.plot_file,
            # Measured from the land mask. The southern wall is 311 wet cells of 351, to 317 m; the
            # northern one is dry. The western one is wet for 14 cells, 2.8 km at 58.95-58.97°N off
            # Nevlunghavn, to 41 m — all of it inside the sixteen-cell `BoundarySponge` of the
            # southern edge it meets, so the westward outflow METreport 4/2016 §5.1 describes
            # turning "inside of Store Færder" leaves through the southern boundary a few cells from
            # the corner rather than being blocked. The eastern one is wet for 26 cells at
            # 59.07-59.11°N, to 66 m, where the Idefjord approach meets the Swedish border, and the
            # edge of what the Norwegian FileGDB covers: a sill fjord with no large exchange, on
            # data that thins out at the wall.
            #
            # Opening either would cost far more than it buys. `boundary_domain` takes the *bounding
            # box* of the open edges' bands, and a full-width southern band with a full-height
            # western or eastern one spans the whole domain — about 110 GB of hourly download over
            # two years against roughly 5 GB for the southern edge alone.
            open_edges       = :south,
            margin           = base.boundary_config.margin,
            architecture     = base.boundary_config.architecture,
            # `ubar_eastward`/`vbar_northward`, not the operational collection's `ubar`/`vbar`:
            # these come from the `sdepth` files and are already rotated to geographic axes.
            # `boundary_variable_names` maps them onto the `ubar`/`vbar` the read side expects.
            parameters       = [
                "temperature", "salinity", "u_eastward", "v_northward", "zeta",
                "ubar_eastward", "vbar_northward",
            ],
            years            = years,
        ),
        VIKER_TIDE_GAINS,
    ),
    atmosphere_config = atmosphere,
    simulation_config = SimulationConfig(
        results_root        = fjord_results_root("oslofjorden_validation"),
        architecture        = simulation.architecture,
        model               = model,
        boundary_conditions = simulation.boundary_conditions,
        writers = (
            # Daily, not `oslofjorden.jl`'s three-hourly. One record of five fields on this grid is
            # 155 MB, so three-hourly over 640 days is about 790 GB against 99 GB daily. The full
            # fields are here for maps, for the two Statnett transect sections and for the M2
            # amplitude and phase maps, none of which needs sub-daily sampling — everything that
            # does is a station below.
            SnapshotWriter(
                name               = :ocean,
                output_file        = "snapshots_ocean.nc",
                variables          = (:T, :S, :u, :v, :e),
                interval           = 1day,
                overwrite_existing = true,
            ),
            FieldSnapshotWriter(
                name               = :surface,
                output_file        = "snapshots_surface.jld2",
                variables          = (:η,),
                interval           = 1day,
                overwrite_existing = true,
            ),
            # Ten minutes, and the only writer at that cadence. A tide gauge comparison is a
            # harmonic analysis: METreport 11/2017 Table 3 resolves constituents down to MS4 at
            # 6.1033 hours, and the shortest of them needs several samples per period to come back
            # unaliased. Three points at 10 minutes over two years is about 1 MB, so the cadence
            # costs nothing; at full-domain resolution it would be unthinkable.
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
            # Hourly profiles of velocity at the seven ADCP sites. Hourly is what the observations
            # are and what the report's own analysis assumes: the 49-hour running mean it separates
            # the estuarine circulation with, and the tidal ellipses of its Table 4, both need the
            # tide resolved in the record they are taken from.
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
            # cast dates needs — but it is also what the Hovmøller diagrams of METreport 11/2017
            # Figs. 25-28 are drawn from, and what resolves a deep-water renewal event, which
            # happens over days and is the thing this domain is judged on. `e` rides along because
            # the reports' two standing conclusions are both about vertical mixing — too weak in the
            # open fjord, too vigorous inside the sills — and the TKE is what would show it.
            StationWriter(
                name               = :hydrography,
                output_file        = "stations_hydrography.nc",
                variables          = (:T, :S, :e),
                stations           = vcat(OSLOFJORD_CTD_STATIONS, INNER_OSLOFJORD_STATIONS),
                interval           = 1hour,
                overwrite_existing = true,
            ),
            # Hourly whole-column profiles along the Drammensfjord section: T and S for the
            # mixing and renewal, `e` for the TKE, and u, v for the exchange over the sill. About
            # 12 MB a station over 640 days, 0.9 GB for the 77.
            StationWriter(
                name               = :section,
                output_file        = "stations_section.nc",
                variables          = (:T, :S, :e, :u, :v),
                stations           = DRAMMENSFJORD_SECTION,
                interval           = 1hour,
                overwrite_existing = true,
            ),
            # The fixed-point temperature sites. Hourly matches the Scanmar mooring's own cadence,
            # and the beaches' three-hourly daytime sampling is a subset of it. The whole column is
            # written rather than the top cell alone, because the report's explanation for the beach
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
        callbacks           = simulation.callbacks,
        time_stepping       = simulation.time_stepping,
        # The previous full run's own state on 1 April 2015 (stamp 20261014T142318, archived as
        # `oslofjorden_validation_2026-10-16`), a year of this model's spin-up, rather than the
        # Norkyst-v3 state of 31 March 2014. Against the 9 April 2014 Inner Oslofjord casts it halves
        # the temperature error below 30 m — RMSE 1.52 -> 0.77 °C at 30-60 m and 2.20 -> 1.18 °C
        # below, Norkyst's deep inner water being about 2 °C too cold — and below 30 m is where the
        # start matters, since the sill basins are not renewed and keep what they start with.
        # Salinity below 60 m is a tie (1.17 against 1.14 g/kg), its bias turned from +0.7 to -1.0
        # by that run's own freshening, which also makes 10-30 m worse (0.86 -> 1.98); the upper
        # layer forgets its start within weeks. The outer fjord's first casts, ten weeks in, tie
        # below 60 m.
        #
        # So this run no longer starts the way FjordOs did, from the "semi-hot start" of a coarser
        # NorKyst (METreport 4/2016 §5): the comparison with it is of two runs, not of one start.
        #
        # The file is that run's snapshot record cut out by `scripts/extract_initial_state.jl`, so a
        # GCP run stages it like any prepared input; with no date, `FromResults` reads its only
        # record. The free surface and the turbulence field start from their defaults.
        initial_conditions  = FromResults(joinpath(data_root, "initial_state_20150401.nc")),
        # 1 April 2014 to the end of 31 December 2015, the FjordOs CL hindcast window exactly: 275
        # days of 2014 plus the whole of 2015, so the last day is included rather than ending at
        # midnight on it. The observation campaigns sit well inside it: the Statnett moorings are
        # mid-September to late November 2014, five months in; the CTD casts and the Slagen current
        # year are 2015, nine months in and later. That is the same spin-up margin FjordOs had,
        # which matters because spin-up is one of the things being compared.
        start_date          = DateTime(2014, 4, 1),
        stop_time           = 640days,
        loops               = 1,
        pickup              = false,
    ),
)
