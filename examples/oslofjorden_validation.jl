# The Oslofjord setup rebuilt for validation against the FjordOs reports: `oslofjorden()`'s physics,
# over the domain and the years those reports cover, with the in-run diagnostics a model-observation
# comparison needs.
#
# It exists to answer one question — how does FjordSim, as it actually stands, score against the same
# observations MET Norway scored their ROMS model FjordOs CL against? — so **nothing about the physics
# is changed or retuned here**. The references are Røed et al., *A high-resolution, curvilinear ROMS
# model for the Oslofjord*, METreport 4/2016, which documents FjordOs CL and its 1 April 2014 -
# 31 December 2015 hindcast, and Hjelmervik et al., *Evaluation of the FjordOs-model*, METreport
# 11/2017, which evaluates it against water level, currents, CTD profiles and fixed-point temperature.
#
# Nothing in `src/` is specific to this file. What it uses from the package is general: the
# Norkyst-v3 hindcast configs, the station writers, and `validate_simulation` with its Kartverket
# water-level source.
#
#   julia --project -m FjordSim prepare_bathymetry  --config examples/oslofjorden_validation.jl
#   julia --project -m FjordSim download_forcing    --config examples/oslofjorden_validation.jl
#   julia --project -m FjordSim prepare_forcing     --config examples/oslofjorden_validation.jl
#   julia --project -m FjordSim add_rivers          --config examples/oslofjorden_validation.jl
#   julia --project -m FjordSim download_boundaries --config examples/oslofjorden_validation.jl
#   julia --project -m FjordSim prepare_boundaries  --config examples/oslofjorden_validation.jl
#   julia --project -m FjordSim download_atmosphere --config examples/oslofjorden_validation.jl
#   julia --project -m FjordSim prepare_atmosphere  --config examples/oslofjorden_validation.jl
#   julia --project -m FjordSim run_simulation      --config examples/oslofjorden_validation.jl
#   julia --project -m FjordSim validate_simulation --config examples/oslofjorden_validation.jl
#
# `add_rivers` needs `NVE_API_KEY`, exactly as `oslofjorden()` does.
#
# # What is reused from `oslofjorden()`, and why it is reused rather than copied
#
# The float type, the closure (every coefficient, the `BoundarySponge`, the CATKE `minimum_tke`
# floor), the advection schemes (salinity's excepted, below), the boundary conditions and their two Marchesiello timescales, the
# time stepping, the initial-condition rule, every bathymetry smoothing knob and the whole river
# configuration are taken from the config `oslofjorden()` returns, not restated. A validation of a model that has drifted from
# the model being validated is worth nothing, and a copy is exactly what drifts. So a change to any of
# them is made in `oslofjorden()` and arrives here by itself — but a bathymetry knob is a *data*
# change, and needs this file's `prepare_bathymetry`, `prepare_forcing`, `add_rivers` and
# `prepare_boundaries` re-run as much as `oslofjorden()`'s.
#
# # What is different, and why
#
# **The dataset.** `oslofjorden()` reads the operational NorKyst-800m archive at
# `fou-hi/norkyst800m`, which begins 2017-02-20 and so cannot reach this window at all. This reads the
# Norkyst-v3 hindcast at `romshindcast/norkyst_v3` (`NorKystHindcastConfig`,
# `NorKystHindcastBoundariesConfig`), MET's continuous free run from 2012, which carries the same
# variables on the same grid. The boundary half gains from the change: the hindcast publishes
# `ubar_eastward`/`vbar_northward` already rotated to geographic axes, where the operational
# collection publishes ROMS' grid-relative pair and `NorKystBoundariesConfig` has to derotate it.
#
# **The domain.** 10.00-11.20°E by 58.95-59.93°N against `oslofjorden()`'s 10.2-11.02 by 59.0-59.93 —
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
# **The writers.** A 640-day run cannot afford `oslofjorden()`'s three-hourly full-domain snapshots —
# at 4.4 M cells a record is 89 MB, so that is roughly 450 GB — and a validation does not need them:
# what it needs is high cadence at a few dozen points. So the 3D fields drop to daily (57 GB) for maps
# and sections, and four station writers carry the time resolution at about 150 MB total.
#
# **The data root.** Only the Geonorge FileGDB is shared with `oslofjorden()`, by absolute path,
# because the 4.5 GB national download is the same for any Norwegian domain. Everything else is a
# different period and downloads into this file's own root.
#
# # Two gates, settled by measurement
#
# Both were measured on the `bathymetry.nc` `prepare_bathymetry` writes, with `max_sounding_distance`
# in force; the reasoning is beside each value below. A change to the box or to any bathymetry knob
# reopens both.
#
# 1. **`z_faces`.** The deepest sounding is 417.9 m, so `oslofjorden()`'s rule gives its own faces
#    plus one 33.5 m layer, to -433.5 m.
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
# pressure-gradient error that forces it does not exist here. With `max_slope_factor = 0.25`, the
# only slope limit this bathymetry gets costs 0.6 m at the Drøbak sill against a 395.1 m basin. That
# is also why `plot_section` draws no smoothed-against-real bathymetry comparison, which is the main
# point of Fig. 11: there is nothing to compare. Its two sections are the Statnett transects — the
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
using Oceananigans: WENO
using Oceananigans.Units
using Dates: DateTime

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
# diagrams for; OF-1, LA-1 and D-2 the three it overlays profiles at.
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
# of Fagrådet since 1973. Not in either METreport, so these are APPROXIMATE and are to be replaced
# from the programme's own station register. They are here because deep-water renewal in these two
# basins is what `docs/setups.md` names as "what this domain is judged on", and neither report
# covers it — this is the part of the validation that goes beyond reproducing FjordOs.
const INNER_OSLOFJORD_STATIONS = [
    Station(name = "Dk1", longitude = 10.5717, latitude = 59.7217),  # Vestfjorden,  ~105 m
    Station(name = "Ep1", longitude = 10.7217, latitude = 59.8117),  # Bunnefjorden, ~155 m
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

# --- The config -------------------------------------------------------------------------------

base = oslofjorden()
data_root = fjord_data_root("oslofjorden_validation")
years = [2014, 2015]

# The three configs below are `oslofjorden()`'s own instances with only their location and period
# changed, rather than fresh ones restating its fields: every smoothing knob, every river override and
# threshold, and the atmosphere's resolution then carry over by construction, including any field
# added later. `oslofjorden()` builds new instances on every call, so this mutates nothing shared.

bathymetry = base.bathymetry_config
# Oslofjord's extracted copy of the national FileGDB, by absolute path — the same sharing
# `drammensfjorden()` does. Resolved before `data_root` moves, which it is relative to.
bathymetry.geodatabase_file = geodatabase_path(bathymetry)
bathymetry.data_root = data_root
# This box reaches two places the FileGDB has no data for, and the gridding fills both with water
# extrapolated from soundings tens of kilometres away: the land north of Drammen above 59.76°N west
# of 10.25°E, a gap in the FileGDB's inland coverage, which became a 7 m sea joined to the
# Drammensfjord; and Swedish territory south of Strömstad, which became a 100 m one on the eastern
# wall. Every real wet cell lies within 1.5 km of a sounding, so 2 km removes both and nothing else.
# `oslofjorden()`'s box touches neither and does not set it.
bathymetry.max_sounding_distance = 2000.0

# Discovery will find more mouths than Oslofjord's 21 because the box is larger — Numedalslågen in
# the west, and Glomma's Singlefjord side in the east — which is the intended effect. The FjordOs
# river table has 37 named rivers over a comparable domain, so the two counts should end up close.
# Every override `oslofjorden()` states lies well inside this box, which only grew, so none of them
# matches no mouth.
rivers = base.forcing_config.rivers
rivers.data_root = data_root
rivers.years = years

# Three of `oslofjorden()`'s gauges were chosen for 2020 and have no daily values on HydAPI in 2014
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
#   left to `river_lambdas`' fallback, the REGINE catchment normal.
restated = Dict(
    "012.A2" => NVERiver(
        vassdragsnr = "012.A2", name = "Drammenselva",
        discharge_station = "12.534.0", temperature_station = "12.15.0",
        plume_depth = Inf,
    ),
    "006.A10" => NVERiver(vassdragsnr = "006.A10", name = "Akerselva", discharge_station = "6.9.0"),
    "007.A0" => NVERiver(vassdragsnr = "007.A0"),
)
rivers.outlets = [get(restated, outlet.vassdragsnr, outlet) for outlet in rivers.outlets]

atmosphere = base.atmosphere_config
atmosphere.data_root = data_root
atmosphere.years = years

simulation = base.simulation_config

# The one physics departure from `oslofjorden()`: salinity is advected by bounds-preserving WENO.
# Plain WENO is not bounded, and this domain adds Numedalslågen at Larvik, a large river relaxed to
# S = 0 among two-level columns. The horizontal reconstruction across that front drove the
# neighbouring surface cells to -29 psu over a 110 psu cell beneath by day 1.5, and past -32 a step
# later, where TEOS10's `sqrt(Sᴬ + 32)` threw a DomainError on the GPU. Measured on a Larvik
# sub-domain with no atmosphere: WENO order 3 or 5, `PartialCellBottom` or `GridFittedBottom`,
# CATKE or a constant diffusivity all reproduce it, and switching the river off, upwinding only the
# horizontal, or bounding S each removes it. T and `e` keep the base's scheme. Water below 0 °C is
# real here, and a NamedTuple must name `e`, or Oceananigans gives it the unbounded `Centered()`
# default.
#
# Not changed in `oslofjorden()` itself, because a NamedTuple there would leave the NPZD example's
# tracers on `Centered()`.
scheme = simulation.model.tracer_advection
model = CoupledHydrostaticSimulation(
    buoyancy           = simulation.model.buoyancy,
    closure            = simulation.model.closure,
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
        # latitude — the same cell size as `oslofjorden()`, deliberately. Everything tuned against Δx
        # or Δx⁴ carries over only if Δx does: the biharmonic coefficients are justified by a damping
        # rate that goes as ν₄·16/Δx⁴, and changing the resolution without re-deriving them would
        # silently invalidate the one number in that setup that was measured rather than chosen.
        # 4.43 M cells against Oslofjord's 3.00 M.
        size      = (351, 548, 25),
        halo      = base.grid_config.halo,
        longitude = (10.00, 11.20),
        latitude  = (58.95, 59.93),
        # `oslofjorden()`'s 24 levels reach -400 m, which clears that domain's deepest sounding of
        # 395.1 m. This box extends into the Skagerrak and takes in the Hvalerdjupet, which
        # METreport 4/2016 §1.1 describes as "a 400 m deep basin extending northeastward from the
        # Skagerrak", and its deepest sounding is 417.9 m. `oslofjorden()`'s rule — a geometric
        # stretch of ratio 1.25 from a 1 m surface cell, capped at 33.5 m, with just enough 33.5 m
        # layers to clear the deepest sounding — therefore adds exactly one layer.
        #
        # Both errors cost something. Too shallow and the basin is silently truncated —
        # `snap_partial_bottom_cells` skips a sounding below the deepest face and `PartialCellBottom`
        # clips it, so nothing complains. Too deep and `grid.Lz` feeds `sqrt(g·Lz)` in
        # `SplitExplicitFreeSurface`, buying a shorter barotropic substep for water that is not
        # there. `z_faces` is written into `bathymetry.nc` and `simulation_grid` reads the file, not
        # this config, so a change here needs `prepare_bathymetry` re-run.
        z_faces   = [
            -433.5, -400.0, -366.5, -333.0, -299.5, -266.0, -232.5, -199.0, -165.5,
            -132.0, -105.0, -83.0, -66.0, -52.0, -41.0, -32.0, -25.0,
            -19.0, -14.5, -10.8, -7.9, -5.5, -3.7, -2.2, -1.0, 0.0,
        ],
        # `oslofjorden()`'s, since the model below is built in it.
        float_type = base.grid_config.float_type,
    ),
    bathymetry_config = bathymetry,
    forcing_config = NorKystHindcastConfig(
        data_root        = data_root,
        output_directory = "norkyst_v3",
        output_file      = base.forcing_config.output_file,
        plot_file        = base.forcing_config.plot_file,
        architecture     = base.forcing_config.architecture,
        parameters       = base.forcing_config.parameters,
        years            = years,
        rivers           = rivers,
    ),
    boundary_config = NorKystHindcastBoundariesConfig(
        data_root        = data_root,
        output_directory = "norkyst_v3_hourly",
        output_file      = base.boundary_config.output_file,
        plot_file        = base.boundary_config.plot_file,
        # Measured from the land mask. The southern wall is 311 wet cells of 351, to 317 m; the
        # northern one is dry. The western one is wet for 14 cells, 2.8 km at 58.95-58.97°N off
        # Nevlunghavn, to 41 m — all of it inside the sixteen-cell `BoundarySponge` of the southern
        # edge it meets, so the westward outflow METreport 4/2016 §5.1 describes turning "inside of
        # Store Færder" leaves through the southern boundary a few cells from the corner rather than
        # being blocked. The eastern one is wet for 26 cells at 59.07-59.11°N, to 66 m, where the
        # Idefjord approach meets the Swedish border, and the edge of what the Norwegian FileGDB
        # covers: a sill fjord with no large exchange, on data that thins out at the wall.
        #
        # Opening either would cost far more than it buys. `boundary_domain` takes the *bounding
        # box* of the open edges' bands, and a full-width southern band with a full-height western
        # or eastern one spans the whole domain — about 110 GB of hourly download over two years
        # against roughly 5 GB for the southern edge alone.
        open_edges       = :south,
        margin           = base.boundary_config.margin,
        architecture     = base.boundary_config.architecture,
        # `ubar_eastward`/`vbar_northward`, not the operational collection's `ubar`/`vbar`: these
        # come from the `sdepth` files and are already rotated to geographic axes.
        # `boundary_variable_names` maps them onto the `ubar`/`vbar` the read side expects.
        parameters       = [
            "temperature", "salinity", "u_eastward", "v_northward", "zeta",
            "ubar_eastward", "vbar_northward",
        ],
        years            = years,
    ),
    atmosphere_config = atmosphere,
    simulation_config = SimulationConfig(
        results_root        = fjord_results_root("oslofjorden_validation"),
        architecture        = simulation.architecture,
        model               = model,
        boundary_conditions = simulation.boundary_conditions,
        writers = (
            # Daily, not `oslofjorden()`'s three-hourly. One record of five fields on this grid is
            # 89 MB, so three-hourly over 640 days is about 450 GB against 57 GB daily. The full
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
            # Figs. 25-28 are drawn from, and what resolves a deep-water renewal event, which happens
            # over days and is the thing this domain is judged on. `e` rides along because the
            # reports' two standing conclusions are both about vertical mixing — too weak in the open
            # fjord, too vigorous inside the sills — and the TKE is what would show it.
            StationWriter(
                name               = :hydrography,
                output_file        = "stations_hydrography.nc",
                variables          = (:T, :S, :e),
                stations           = vcat(OSLOFJORD_CTD_STATIONS, INNER_OSLOFJORD_STATIONS),
                interval           = 1hour,
                overwrite_existing = true,
            ),
            # The fixed-point temperature sites. Hourly matches the Scanmar mooring's own cadence,
            # and the beaches' three-hourly daytime sampling is a subset of it. The whole column is
            # written rather than the top cell alone, because the report's explanation for the beach
            # discrepancy is a diurnal cycle in the model that is too strong, and that is a statement
            # about the surface layer rather than the surface.
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
        # The Norkyst-v3 state at the start, which is the same "semi-hot start" from a coarser
        # NorKyst that METreport 4/2016 §5 used — and, per METreport 11/2017 §5, one of the two
        # things its authors blamed for FjordOs' stratification errors. Starting the same way is what
        # makes the comparison a comparison.
        #
        # From the record twelve hours before `start_date` rather than at it, because the hindcast
        # forcing keeps one record a day at 12:00 and `FromForcing()` needs one exactly at midnight.
        # Moving `start_date` to noon instead would shift the FjordOs window; half a day of state
        # under five months of spin-up is the smaller departure.
        initial_conditions  = FromForcing(DateTime(2014, 3, 31, 12)),
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
