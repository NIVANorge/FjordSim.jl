# Station positions the validation writes model output at, grouped by the observation programme they
# come from. Kept at module scope rather than inside `oslofjorden_validation()` because the analysis
# side reads the same lists to know which observations to fetch.
#
# Provenance matters here and is recorded per group: a coordinate that came from a published table
# is exact, and one that came from a description in prose is not. `StationWriter` logs how far each
# station moved when it was snapped to the grid, which is the number that says whether an
# approximate position mattered.

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

"""
    oslofjorden_validation()

The Oslofjord setup rebuilt for validation against the FjordOs reports: the same physics as
[`oslofjorden`](@ref), over the domain and the years those reports cover, with the in-run
diagnostics a model-observation comparison needs.

It exists to answer one question — how does FjordSim, as it actually stands, score against the same
observations MET Norway scored their ROMS model FjordOs CL against? — so **nothing about the physics
is changed or retuned here**. Every closure coefficient, the `BoundarySponge`, the CATKE
`minimum_tke` floor, the scalar `tracer_advection`, the two Marchesiello timescales and the whole
river configuration are carried over from `oslofjorden()` unchanged. What differs is the domain, the
period, the dataset the two come from, and what the run writes.

The two reports are Røed et al., *A high-resolution, curvilinear ROMS model for the Oslofjord*,
METreport 4/2016, which documents FjordOs CL and its 1 April 2014 - 31 December 2015 hindcast, and
Hjelmervik et al., *Evaluation of the FjordOs-model*, METreport 11/2017, which evaluates it against
water level, currents, CTD profiles and fixed-point temperature. The station lists above are their
Tables 1 and 2.

# What is different from `oslofjorden()`, and why

**The dataset.** `oslofjorden()` reads the operational NorKyst-800m archive, which begins
2017-02-20 and so cannot reach this window at all. This setup reads the Norkyst-v3 hindcast
(`NorKystHindcastConfig`, `NorKystHindcastBoundariesConfig`), MET's continuous free run from 2012,
which carries the same variables on the same grid. The boundary half gains from the change: the
hindcast publishes `ubar_eastward`/`vbar_northward` already rotated to geographic axes, where the
operational collection publishes ROMS' grid-relative pair and `NorKystBoundariesConfig` has to
derotate it.

**The domain.** 10.00-11.20°E by 58.95-59.93°N against `oslofjorden()`'s 10.2-11.02 by 59.0-59.93 —
the box FjordOs CL covers, at the same cell size. The three extensions each buy something the
validation needs and nothing else does:

- *East to 11.20°E* brings in the Hvaler archipelago, both arms of Glomma's delta and CTD station
  S-9. Glomma is the largest river in the domain, and its outlet being in the *outer* fjord rather
  than at the head is what METreport 4/2016 §1.1 calls out as making the Oslofjord's estuarine
  circulation "deviate considerably from a classical textbook example". Cutting the delta in half is
  not a small approximation.
- *West to 10.00°E* brings in Numedalslågen — 6514 km², the third largest catchment in the FjordOs
  river table — and CTD station LA-1.
- *South to 58.95°N* is for LA-1 specifically. At 59.019°N it would sit about ten cells inside a
  sixteen-cell `BoundarySponge` on the old southern edge and be unusable; here it is about thirty-nine
  cells in. METreport 11/2017 hits the same problem from the other side, noting LA-1 is "close to
  the southern open boundary of the FjordOs model".

**The writers.** A 640-day run cannot afford `oslofjorden()`'s three-hourly full-domain snapshots —
at 4.6 M cells that is roughly 470 GB — and a validation does not need them: what it needs is high
cadence at a few dozen points. So the 3D fields drop to daily for maps and sections, and the station
writers carry the time resolution at about 150 MB total. See the `writers` tuple below.

# What this setup can and cannot be expected to show

Stated up front because it decides how the results should be read. At ~193 m this grid is finer than
NorKyst-800 everywhere and much finer in the vertical, but it is two to three times coarser than
FjordOs CL in exactly the places FjordOs was built for. The Drøbak Sound is 1-2 km wide, six to ten
cells here against fifteen to twenty-five there; the Drøbak Jetty's two openings are ~6 m deep and
tens of metres wide and are not representable at all; Svelvik is 180 m wide and 11 m deep, one cell.
So the sill-controlled tidal jets METreport 4/2016 features, and the Drammensfjord exchange through
Svelvik, are out of reach on this grid whatever is tuned.

What should be competitive or better is everything that is not sill-controlled: tidal elevation,
open-fjord currents (Slagen above all), the seasonal hydrography at the CTD stations, and deep-water
renewal in the inner basins. The z-coordinate with `PartialCellBottom` removes the defect both
reports complain about most — FjordOs' bathymetry had to be smoothed to an rx0 limit until "the
observed slopes are everywhere steeper" than the model's (METreport 11/2017 Fig. 11), and the
pressure-gradient error that forces it does not exist here. Establishing which half of that
expectation holds, quantitatively, is the point of the run.

# Preparation

Shares only the Geonorge FileGDB with `oslofjorden()`, by absolute path, because the 4.5 GB national
download is the same for any Norwegian domain. Everything else is a different period and downloads
into this setup's own root, which is what the per-setup `data_root` convention requires.

`prepare_bathymetry` must run before the two gates below are settled — see `z_faces` and
`open_edges` in the config.
"""
function oslofjorden_validation()
    data_root = joinpath(homedir(), "FjordSim_data", "oslofjorden_validation")
    oslofjorden_data_root = joinpath(homedir(), "FjordSim_data", "oslofjorden")
    FT = Oceananigans.defaults.FloatType

    return FjordConfig(
        grid_config = EvenGrid(
            # 351 x 548 cells over 1.20° x 0.98° is 193.8 m in longitude at 59.4°N and 199.1 m in
            # latitude — the same cell size as `oslofjorden()`, deliberately. Everything tuned
            # against Δx or Δx⁴ carries over only if Δx does: the biharmonic coefficients below are
            # justified by a damping rate that goes as ν₄·16/Δx⁴, and changing the resolution
            # without re-deriving them would silently invalidate the one number in this file that
            # was measured rather than chosen. 4.61 M cells against Oslofjord's 3.00 M.
            size      = (351, 548, 26),
            halo      = (7, 7, 7),
            longitude = (10.00, 11.20),
            latitude  = (58.95, 59.93),
            # PROVISIONAL, and the first of two gates this setup has to pass before it is run.
            #
            # `oslofjorden()`'s 24 levels reach -400 m, which clears that domain's deepest sounding
            # of 395.1 m. This box extends into the Skagerrak and takes in the Hvalerdjupet, which
            # METreport 4/2016 §1.1 describes as "a 400 m deep basin extending northeastward from
            # the Skagerrak", so the floor here is deeper and by an unknown amount. Two 33.5 m
            # layers are added at the bottom as headroom while that is measured.
            #
            # Settle it by running `prepare_bathymetry`, reading the deepest sounding it reports,
            # and regenerating this list by the rule `oslofjorden()` used: a geometric stretch of
            # ratio 1.25 from a 1 m surface cell, capped at 33.5 m, with just enough 33.5 m layers
            # to clear that sounding. Both errors cost something. Too shallow and the basin is
            # silently truncated — `snap_partial_bottom_cells` skips a sounding below the deepest
            # face and `PartialCellBottom` clips it, so nothing complains. Too deep and `grid.Lz`
            # feeds `sqrt(g·Lz)` in `SplitExplicitFreeSurface`, buying a shorter barotropic substep
            # for water that is not there. Then re-run `prepare_bathymetry`, since `z_faces` is
            # written into `bathymetry.nc` and `simulation_grid` reads the file, not this config.
            z_faces   = [
                -467.0, -433.5, -400.0, -366.5, -333.0, -299.5, -266.0, -232.5, -199.0, -165.5,
                -132.0, -105.0, -83.0, -66.0, -52.0, -41.0, -32.0, -25.0,
                -19.0, -14.5, -10.8, -7.9, -5.5, -3.7, -2.2, -1.0, 0.0,
            ],
        ),
        bathymetry_config = DybdedataConfig(
            data_root             = data_root,
            output_file           = "bathymetry.nc",
            plot_file             = "bathymetry.png",
            # Oslofjord's extracted copy of the national FileGDB, by absolute path — the same
            # sharing `drammensfjorden()` does, and the only thing shared with another setup here.
            geodatabase_file      = joinpath(oslofjorden_data_root, GEONORGE_DYBDEDATA_GDB),
            raw_resolution_factor = 2,
            padding_cells         = 2,
            include_contours      = false,
            contour_stride        = 10,
            interpolation_passes  = 1,
            major_basins          = 1,
            minimum_depth         = 2.0,
            open_boundary_land_cells = 5,
            max_island_cells      = 6,
            close_narrow_passages = true,
            spike_ratio           = 0.5,
            minimum_cell_fraction = 0.2,
            max_slope_factor      = 0.25,
            geonorge_cache        = true,
            regrid_cache          = false,
        ),
        forcing_config = NorKystHindcastConfig(
            data_root        = data_root,
            output_directory = "norkyst_v3",
            output_file      = "forcing.nc",
            plot_file        = "forcing.png",
            architecture     = :auto,
            parameters       = ["temperature", "salinity", "u_eastward", "v_northward"],
            years            = [2014, 2015],
            rivers           = NVERiversConfig(
                data_root  = data_root,
                output_file = "forcing_rivers_nve.nc",
                plot_file   = "forcing_rivers_nve.png",
                years       = [2014, 2015],
                # Unchanged from `oslofjorden()`, including the threshold. Discovery will find more
                # mouths than Oslofjord's 21 because the box is larger — Numedalslågen in the west,
                # and Glomma's Singlefjord side in the east — which is the intended effect. The
                # FjordOs river table has 37 named rivers over a comparable domain, so the two
                # counts should end up close.
                minimum_discharge = 0.5,
                default_plume_depth = 5.0,
                minimum_levels = 0,
                minimum_relaxation_timescale = 600.0,
                # The same eight overrides, verbatim from `oslofjorden()`. Each attaches an observed
                # NVE gauge to a mouth discovery finds anyway; three of them also state geometry.
                # An override matching no discovered mouth is an error, and all eight lie well
                # inside the box, which only grew.
                outlets = [
                    NVERiver(
                        vassdragsnr = "002.A21", name = "Glomma (Osterelva)",
                        discharge_station = "2.605.0", temperature_station = "2.1087.0",
                        discharge_fraction = 2 // 3, plume_depth = Inf,
                    ),
                    NVERiver(
                        vassdragsnr = "002.2A", name = "Glomma (Vesterelva)",
                        discharge_station = "2.605.0", temperature_station = "2.1087.0",
                        discharge_fraction = 1 // 3, plume_depth = Inf,
                    ),
                    NVERiver(
                        vassdragsnr = "012.A2", name = "Drammenselva",
                        discharge_station = "12.534.0", temperature_station = "12.534.0",
                        plume_depth = Inf,
                    ),
                    NVERiver(vassdragsnr = "003.A4", discharge_station = "3.23.0"),
                    NVERiver(
                        vassdragsnr = "006.A10", name = "Akerselva",
                        discharge_station = "6.38.0",
                    ),
                    NVERiver(vassdragsnr = "007.A0", discharge_station = "7.29.0"),
                    NVERiver(
                        vassdragsnr = "008.A2",
                        discharge_station = "8.2.0", temperature_station = "8.2.0",
                    ),
                    NVERiver(vassdragsnr = "011.A0", discharge_station = "11.6.0"),
                ],
            ),
        ),
        boundary_config = NorKystHindcastBoundariesConfig(
            data_root        = data_root,
            output_directory = "norkyst_v3_hourly",
            output_file      = "boundaries.nc",
            plot_file        = "boundaries.png",
            # PROVISIONAL, and the second gate. `oslofjorden()`'s box has land on all three of its
            # other walls; this one does not obviously. At 10.00°E the coast near Nevlunghavn leaves
            # water on the western wall somewhere south of about 59.05°N, and at 11.20°E the Hvaler
            # and Singlefjord approach leaves water on the eastern one. A closed wall standing in
            # water is a real error rather than a cosmetic one — the southwest corner is where
            # METreport 4/2016 §5.1 has the fjord's outflow turning west "inside of Store Færder to
            # join the westward flowing current in the Skagerrak", and a wall there blocks it.
            #
            # Settle it from the land mask `prepare_bathymetry` writes, not by assumption: measure
            # the wet run along each wall and open the edges where it is long. The cost is why this
            # is not simply set to all three now. `boundary_domain` takes the *bounding box* of the
            # open edges, so adjacent edges (south and east) give a corner strip while opposite ones
            # (west and east) give the whole domain — about 110 GB of hourly download over two
            # years against roughly 5 GB for a single edge.
            open_edges       = :south,
            margin           = 0.05,
            architecture     = :auto,
            # `ubar_eastward`/`vbar_northward`, not the operational collection's `ubar`/`vbar`:
            # these come from the `sdepth` files and are already rotated to geographic axes.
            # `boundary_variable_names` maps them onto the `ubar`/`vbar` the read side expects.
            parameters       = [
                "temperature", "salinity", "u_eastward", "v_northward", "zeta",
                "ubar_eastward", "vbar_northward",
            ],
            years            = [2014, 2015],
        ),
        atmosphere_config = NORA3Config(
            data_root        = data_root,
            output_directory = "nora3",
            output_file      = "atmosphere.nc",
            plot_file        = "atmosphere.png",
            resolution       = 0.02,
            padding          = 0.1,
            years            = [2014, 2015],
        ),
        simulation_config = SimulationConfig(
            results_root       = joinpath(homedir(), "FjordSim_results", "oslofjorden_validation"),
            architecture       = :auto,
            # Verbatim from `oslofjorden()`. The justifications for every number here are in that
            # setup and in `docs/setups.md` under "The 2020 diagnosis"; they are not repeated,
            # because repeating them would invite editing one copy. If any of this changes, it
            # changes there first and this setup follows — a validation of a model that has drifted
            # from the model being validated is worth nothing.
            model              = CoupledHydrostaticSimulation(
                buoyancy           = SeawaterBuoyancy(FT, equation_of_state = TEOS10EquationOfState(FT)),
                closure            = BoundarySponge(
                    base = (
                        CATKEVerticalDiffusivity(minimum_tke = 7e-6),
                        HorizontalScalarBiharmonicDiffusivity(ν = 2e4, κ = 2e3),
                    ),
                    width_cells = 16,
                    viscosity   = 30.0,
                    diffusivity = 15.0,
                ),
                tracer_advection   = WENO(),
                momentum_advection = WENOVectorInvariant(FT),
                tracers            = (:T, :S),
                coriolis           = HydrostaticSphericalCoriolis(FT),
                sea_ice            = FreezingLimitedOceanTemperature(),
                biogeochemistry    = nothing,
                free_surface       = SplitExplicitFreeSurfaceConfig(cfl = 0.7),
                extra_kwargs       = (;),
            ),
            boundary_conditions = MergedBoundaryConditions(
                AirSeaFluxes(),
                QuadraticBottomDrag(coefficient = 0.003),
                OpenLateralBoundaryFromData(
                    inflow_timescale  = 3hours,
                    outflow_timescale = 360days,
                ),
            ),
            writers = (
                # Daily, not `oslofjorden()`'s three-hourly. One record of five fields on this grid
                # is 92 MB, so three-hourly over 640 days is about 470 GB against 59 GB daily. The
                # full fields are here for maps, for the two Statnett transect sections and for the
                # M2 amplitude and phase maps, none of which needs sub-daily sampling — everything
                # that does is a station below.
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
                # 6.1033 hours, and the shortest of them needs several samples per period to come
                # back unaliased. Three points at 10 minutes over two years is about 1 MB, so the
                # cadence costs nothing; at full-domain resolution it would be unthinkable.
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
                # Hourly profiles of velocity at the seven ADCP sites. Hourly is what the
                # observations are and what the report's own analysis assumes: the 49-hour running
                # mean it separates the estuarine circulation with, and the tidal ellipses of its
                # Table 4, both need the tide resolved in the record they are taken from.
                StationWriter(
                    name               = :moorings,
                    output_file        = "stations_moorings.nc",
                    variables          = (:u, :v),
                    stations           = OSLOFJORD_MOORINGS,
                    interval           = 1hour,
                    overwrite_existing = true,
                ),
                # Hourly T/S profiles at the CTD stations and the two inner-basin stations. The
                # observations are a dozen casts a year, so hourly is far more than the comparison
                # at cast dates needs — but it is also what the Hovmøller diagrams of METreport
                # 11/2017 Figs. 25-28 are drawn from, and what resolves a deep-water renewal event,
                # which happens over days and is the thing this domain is judged on. `e` rides
                # along because the reports' two standing conclusions are both about vertical
                # mixing — too weak in the open fjord, too vigorous inside the sills — and the TKE
                # is what would show it.
                StationWriter(
                    name               = :hydrography,
                    output_file        = "stations_hydrography.nc",
                    variables          = (:T, :S, :e),
                    stations           = vcat(OSLOFJORD_CTD_STATIONS, INNER_OSLOFJORD_STATIONS),
                    interval           = 1hour,
                    overwrite_existing = true,
                ),
                # The fixed-point temperature sites. Hourly matches the Scanmar mooring's own
                # cadence, and the beaches' three-hourly daytime sampling is a subset of it. The
                # whole column is written rather than the top cell alone, because the report's
                # explanation for the beach discrepancy is a diurnal cycle in the model that is too
                # strong, and that is a statement about the surface layer rather than the surface.
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
            callbacks = (ProgressCallback(name = :progress, interval = 1hour, report = progress),),
            time_stepping = AdaptiveTimeStep(
                initial_time_step    = 1second,
                cfl                  = 0.3,
                max_time_step        = 3minutes,
                max_time_step_change = 1.01,
            ),
            # The Norkyst-v3 state at `start_date`, which is the same "semi-hot start" from a
            # coarser NorKyst that METreport 4/2016 §5 used — and, per METreport 11/2017 §5, one of
            # the two things its authors blamed for FjordOs' stratification errors. Starting the
            # same way is what makes the comparison a comparison.
            initial_conditions = FromForcing(),
            # 1 April 2014 to the end of 31 December 2015, the FjordOs CL hindcast window exactly:
            # 275 days of 2014 plus the whole of 2015, so the last day is included rather than
            # ending at midnight on it. The
            # observation campaigns sit well inside it: the Statnett moorings are mid-September to
            # late November 2014, five months in; the CTD casts and the Slagen current year are
            # 2015, nine months in and later. That is the same spin-up margin FjordOs had, which
            # matters because spin-up is one of the things being compared.
            start_date         = DateTime(2014, 4, 1),
            stop_time          = 640days,
            loops              = 1,
            pickup             = false,
        ),
    )
end
