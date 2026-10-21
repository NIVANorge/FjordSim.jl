# Lista: the open coast of southernmost Norway, from Sokndal past Lindesnes (6.22-7.25°E,
# 57.89-58.34°N), with Fedafjorden, the Flekkefjord sounds, Rosfjorden and the northern slope of the
# Norwegian Trench. A NorKyst-800m initial state, hourly NorKyst-800m data on an open boundary on
# three sides, NVE rivers, NORA3 atmosphere, run from 1 January 2024 to 1 October 2025.
#
# # Preparing and running it locally
#
#   julia --project -m FjordSim prepare_bathymetry     --config examples/lista.jl
#   julia --project -m FjordSim download_initial_state --config examples/lista.jl
#   julia --project -m FjordSim prepare_initial_state  --config examples/lista.jl
#   julia --project -m FjordSim add_rivers             --config examples/lista.jl
#   julia --project -m FjordSim download_boundaries    --config examples/lista.jl
#   julia --project -m FjordSim prepare_boundaries     --config examples/lista.jl
#   julia --project -m FjordSim download_atmosphere    --config examples/lista.jl
#   julia --project -m FjordSim prepare_atmosphere     --config examples/lista.jl
#   julia --project -m FjordSim run_simulation         --config examples/lista.jl   # needs a GPU
#
# `add_rivers` needs `NVE_API_KEY` (free at https://hydapi.nve.no/Users). `prepare_bathymetry`
# downloads the national Geonorge FileGDB (4.5 GB) if it is missing. Data goes to
# `~/FjordSim_data/lista/`, and results to `~/FjordSim_results/lista/`.
#
# # Settle these after `prepare_bathymetry`, before the rest
#
# The first was estimated from the raw Geonorge soundings, the second from EMODnet's bathymetry
# (`rest.emodnet-bathymetry.eu`). Neither from the `bathymetry.nc` the model reads. If either
# changes, re-run every prepare step after `prepare_bathymetry`.
#
# 1. **`z_faces`.** The deepest Geonorge sounding in the box is 442.1 m, at 7.010°E 57.893°N on the
#    southern edge. The deepest face, -433.5 m, deliberately cuts the few cells below it; the
#    reasoning is beside `z_faces` below. Check that no more than a handful of cells, all on the
#    southern edge, reach below -433.5 m in `bathymetry.nc`.
# 2. **`open_edges`.** EMODnet has the southern edge wet all the way across (300-444 m), the western
#    edge wet for 93 % of its length (to 352 m) and the eastern edge for its southern 31 % (to 357 m),
#    and the northern edge dry. Check them against the land mask in `bathymetry.png`.
#
# # What follows `examples/oslofjorden_validation.jl`, and what does not
#
# The physics is `oslofjorden.jl`'s, scaled to this grid, except where the Oslofjord validation
# measured a departure that is not specific to its own box. Four carry over: bounded WENO for
# salinity, CATKE's `minimum_tke` at 7e-7, the biharmonic tracer diffusivity κ at 0, and 5.06 m
# levels between 41 and 132 m. Each is argued beside its value. Two do not. Its boundary tide
# correction was fitted to the Norkyst-v3 hindcast at Viker, while this box reads the operational
# NorKyst-800m and holds no tide gauge. Its initial state is a previous run's, which a first run of
# this box cannot have.

using FjordSim
using Oceananigans
using Oceananigans.Units
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using NumericalEarth: FreezingLimitedOceanTemperature
using Dates: DateTime

data_root = fjord_data_root("lista")
FT = Float32
# The hourly NorKyst-800m collection ends on 2025-10-05 on thredds.met.no, so 2025 is the last year
# the boundary data can reach.
years = [2024, 2025]
# `oslofjorden.jl`'s levels from a 1 m surface cell stretching by 1.25, with 41-132 m split into
# eighteen 5.06 m levels as in `examples/oslofjorden_validation.jl`. There `oslofjorden.jl`'s five
# 11-27 m levels let WENO's implicit vertical diffusion, ~|w|Δz, mix the sill basins at 5-10 times
# CATKE's κ, and the split cut the Drammensfjord's salt-budget diffusivity below 41 m by 2.6-7.7
# times in 60 days. The mechanism is not specific to that grid. Here it acts wherever the levels are
# thick, so on the sill basins of Fedafjorden and Rosfjorden where they reach below 41 m, and on the
# coastal current's lower interface. Unmeasured on this grid.
#
# Below 132 m, `oslofjorden.jl`'s 33.5 m layers rather than thicker ones, for the same reason: the
# trench slope is where the coastal current meets the Atlantic water under it.
# They stop at -433.5 m, short of the deepest sounding (gate 1 in the header). Binned per model
# cell, the raw Geonorge soundings reach below -433.5 m in only a few cells, all on the southern
# open edge, and by at most about 9 m. A layer to -467 m would leave those as thin, isolated pits
# beside the open boundary, the kind that sat at the freezing point for most of the Oslofjord
# validation run (TODO.md › "Fill single-cell pits on open edges"). Cut at -433.5 m,
# `PartialCellBottom` clamps them to the face, inside the boundary sponge. It also keeps `grid.Lz`,
# which sets the barotropic substep through `sqrt(g·Lz)`, to water that exists.
#
# 38 levels on 123 x 100 columns is 0.47 M cells.
z_faces = vcat(
    [-433.5, -400.0, -366.5, -333.0, -299.5, -266.0, -232.5, -199.0, -165.5],
    range(-132.0, -41.0; length = 19),
    [-32.0, -25.0, -19.0, -14.5, -10.8, -7.9, -5.5, -3.7, -2.2, -1.0, 0.0],
)

FjordConfig(
    grid_config = EvenGrid(
        # 0.00837° by 0.0045°: about 490 m by 500 m at 58.1°N, 2.5 times Oslofjord's 193 m. Every
        # coefficient below that depends on grid spacing (Δx or Δx⁴) is scaled from `oslofjorden.jl`
        # by that factor, not copied.
        size      = (123, 100, length(z_faces) - 1),
        halo      = (7, 7, 7),
        longitude = (6.22, 7.25),
        latitude  = (57.89, 58.34),
        # See beside `z_faces` above. `z_faces` is written into `bathymetry.nc` and
        # `simulation_grid` reads the file, so a change here needs `prepare_bathymetry` re-run.
        z_faces   = z_faces,
        float_type = FT,
    ),
    bathymetry_config = DybdedataConfig(
        data_root             = data_root,
        output_file           = "bathymetry.nc",
        plot_file             = "bathymetry.png",
        # 4 rather than Oslofjord's 2. The model cells are 2.5 times larger, so a 2x native grid
        # would average each one from only four ~245 m samples. At 4x the native grid is still
        # only 492 x 400.
        raw_resolution_factor = 4,
        padding_cells         = 2,
        include_contours      = false,
        contour_stride        = 10,
        interpolation_passes  = 1,
        major_basins          = 1,
        minimum_depth         = 2.0,
        # Off. The coast runs into the eastern open edge, which is land for its northern two thirds,
        # and the stage floods outwards from the water along the whole edge band. It would advance
        # the coastline several kilometres up that band.
        open_boundary_land_cells = 0,
        max_island_cells      = 6,
        close_narrow_passages = true,
        spike_ratio           = 0.5,
        # Must equal the `minimum_fractional_cell_height` the grid gives `PartialCellBottom`.
        minimum_cell_fraction = 0.2,
        # Oslofjord's limit and sill restoration. Here the limit acts on the trench slope, which
        # drops from about 100 m to 300 m within a few cells. Compare `bathymetry.png` with the raw
        # soundings before trusting it.
        max_slope_factor      = 0.5,
        preserve_sills        = true,
        geonorge_cache        = true,
        regrid_cache          = false,
    ),
    # One NorKyst daily mean, the day `start_date` falls in, for `FromInitialState()` below.
    initial_state_config = NorKystInitialStateConfig(
        data_root        = data_root,
        output_directory = "norkyst_initial",
        output_file      = "initial_state.nc",
        parameters       = ["temperature", "salinity", "u_eastward", "v_northward"],
    ),
    river_config = NVERiversConfig(
        data_root   = data_root,
        output_file = "rivers_nve.nc",
        plot_file   = "rivers_nve.png",
        years       = years,
        # Mouths are discovered from NVE's ELVIS network and REGINE catchments. At 0.5 m³/s
        # that finds 20 of them, carrying 310.8 of the 315.0 m³/s that the 35 mouths above
        # 0.15 m³/s carry. Sira (130.5), Kvina (89.5) and Lygna (39.6) dominate.
        minimum_discharge = 0.5,
        # The vertical Courant bound on the 1 m top cell, as on `oslofjorden`: 1 mm/s is
        # w·Δt/Δz = 0.18 at the 3-minute `max_time_step`. On 490 m cells most rivers need one
        # column; the log of `add_rivers` states how many each got.
        maximum_inflow_velocity = 1.0e-3,
        # Overrides, keyed by each mouth's terminal ELVIS `vassdragsnr` (from the discovery).
        # A gauge's daily series is the river's inflow, so its floods reach the fjord. Every
        # series named here was checked for daily values over the run window.
        outlets = [
            # Røynestad (25.85.0), 1139 km². Kvina is regulated by Sira-Kvina, which diverts
            # its upper catchment to Sira: the gauge's 2024-2025 mean is 20.7 m³/s against
            # 72.6 natural. It is used unscaled, so it does not count the residual catchment
            # below the gauge. Rafossen (25.35.0) has temperature for every day of the window.
            NVERiver(
                vassdragsnr = "025.A0", name = "Kvina",
                discharge_station = "25.85.0", temperature_station = "25.35.0",
            ),
            # Tingvatn (24.9.0) drains only 272 km² of Lygna. Its 2024-2025 mean of 17.4 m³/s
            # is scaled by 2.38, the mouth's REGINE normal (39.64) over the gauge's own
            # (16.66), which gives about 41 m³/s.
            #
            # Kvina's temperature at Rafossen, because both of Lygna's own temperature gauges
            # start on 2024-07-31: Tingvatn (24.9.0) and Lygna ndf. Lygne (24.5.0), which also
            # stops on 2025-07-28. The gap fill copies the nearest recorded day, so January to July
            # 2024 would all have been forced at the temperature of 31 July.
            NVERiver(
                vassdragsnr = "024.A0", name = "Lygna",
                discharge_station = "24.9.0", temperature_station = "25.35.0",
                discharge_fraction = 2.38,
            ),
            # HydAPI has no public discharge for Sira or Sokno in this window, so both run on
            # their REGINE natural normals. Sira's normal understates it by the water
            # Sira-Kvina brings in from Kvina. Temperature is from Sira ndf. Sirdalsvatnet
            # (26.47.0) and Sokno ovf. Litlå (26.53.0), both complete.
            NVERiver(vassdragsnr = "026.A", name = "Sira", temperature_station = "26.47.0"),
            NVERiver(vassdragsnr = "026.4A0", name = "Sokno", temperature_station = "26.53.0"),
        ],
    ),
    boundary_config = NorKystBoundariesConfig(
        data_root        = data_root,
        output_directory = "norkyst_hourly",
        output_file      = "boundaries.nc",
        plot_file        = "boundaries.png",
        # Three sides, unlike a fjord. The Norwegian Coastal Current comes in through the eastern
        # edge and leaves through the western one, so closing either would block it. The
        # bounding box of three edge bands is the whole domain, so the hourly download covers all
        # of it: tens of GB for this window.
        open_edges       = (:south, :west, :east),
        margin           = 0.05,
        architecture     = :auto,
        parameters       = [
            "temperature", "salinity", "u_eastward", "v_northward", "zeta", "ubar", "vbar",
        ],
        years            = years,
    ),
    atmosphere_config = NORA3Config(
        data_root        = data_root,
        output_directory = "nora3",
        output_file      = "atmosphere.nc",
        plot_file        = "atmosphere.png",
        resolution       = 0.02,
        padding          = 0.1,
        years            = years,
    ),
    simulation_config = SimulationConfig(
        results_root       = fjord_results_root("lista"),
        architecture       = :auto,
        model              = CoupledHydrostaticSimulation(
            buoyancy           = SeawaterBuoyancy(FT, equation_of_state = TEOS10EquationOfState(FT)),
            closure            = BoundarySponge(
                base = (
                    # The floor sets a background κ = 0.098·e_min/N wherever CATKE's own `e` is
                    # weaker. On Oslofjord that was 85 % of wet cells. With N measured from NorKyst
                    # columns in the trench here (February and August 2024), 7e-7 gives 4e-6 m² s⁻¹
                    # in the summer pycnocline and 1-3e-5 below 100 m. Oceananigans' default of
                    # 1e-9 gives 1e-8, below molecular heat diffusion, so it means no background
                    # mixing at all.
                    # The open part of the box is flushed by the coastal current within days, so
                    # there the floor hardly matters. It matters in the sill basins of Fedafjorden
                    # and Rosfjorden.
                    #
                    # 7e-7 rather than `oslofjorden.jl`'s 7e-6, as in
                    # `examples/oslofjorden_validation.jl`. There 7e-6 mixed the water columns
                    # towards uniform and freshened the sill basins. Over 640 days at 7e-7 the
                    # salinity RMSE against 216 casts was 1.73 g/kg and the temperature RMSE
                    # 1.22 °C, and summer heat no longer reached too deep. The cost was overheating
                    # in 2-3 m deep cells, which this box has among the skerries at
                    # `minimum_depth = 2.0`.
                    # Unmeasured here: check `e` against the floor in the snapshots.
                    CATKEVerticalDiffusivity(FT; minimum_tke = 7e-7),
                    # Oslofjord's ν = 2e4 m⁴ s⁻¹ scaled by (490/193)⁴. That keeps the same damping
                    # times: about 74 min for the 2Δx mode and about 52 h for an 8Δx (4 km) eddy. The
                    # explicit limit Δt ≤ Δx⁴/32ν₄ is 2200 s.
                    #
                    # κ is 0, as in `examples/oslofjorden_validation.jl`. There κ = 2e3 left the
                    # deep-basin salt budgets and the stations unchanged over 60 days, but
                    # overshot at fronts and left 169 wet cells below 0 g/kg, against none at κ = 0.
                    HorizontalScalarBiharmonicDiffusivity(FT; ν = 8e5, κ = 0),
                ),
                # 8 cells, about 4 km, against Oslofjord's 16 cells (3.2 km). Sixteen would cover
                # a quarter of the domain's width on each side wall.
                width_cells = 8,
                # Explicit-diffusion stability Δt ≤ Δx²/4ν is about 600 s at 100 m² s⁻¹ on this
                # grid, comfortably above `max_time_step`.
                viscosity   = 100.0,
                diffusivity = 50.0,
            ),
            # Salinity bounded, as in `examples/oslofjorden_validation.jl`. Plain WENO drove
            # salinity below TEOS10's -32 psu floor beside a large river entering shallow
            # columns, and Kvina and Lygna enter at the heads of narrow fjords, Sira into a small bay. A
            # NamedTuple has to name `e`, or CATKE's TKE falls back to `Centered()`.
            tracer_advection   = (T = WENO(FT), S = WENO(FT; bounds = (0, 40)), e = WENO(FT)),
            momentum_advection = WENOVectorInvariant(FT),
            tracers            = (:T, :S),
            coriolis           = HydrostaticSphericalCoriolis(FT),
            sea_ice            = FreezingLimitedOceanTemperature(FT),
            biogeochemistry    = nothing,
            free_surface       = SplitExplicitFreeSurfaceConfig(cfl = 0.7),
            extra_kwargs       = (;),
        ),
        boundary_conditions = MergedBoundaryConditions(
            AirSeaFluxes(),
            QuadraticBottomDrag(coefficient = 0.003),
            # Marchesiello et al. (2001): nudge hard on inflow, radiate on outflow. Oslofjord's
            # timescales, applied to every open edge.
            OpenLateralBoundaryFromData(
                inflow_timescale  = 3hours,
                outflow_timescale = 360days,
            ),
        ),
        writers = (
            # Daily. One record of these five fields is 5.9 MB, so 639 days is 3.8 GB, against
            # 30 GB at three-hourly.
            SnapshotWriter(
                name               = :ocean,
                output_file        = "snapshots_ocean.nc",
                variables          = (:T, :S, :u, :v, :e),
                interval           = 1day,
                overwrite_existing = true,
            ),
            # η goes to JLD2 because the NetCDF writer cannot write a (Center, Center, Nothing)
            # output (see `FieldSnapshotWriter`). Hourly to resolve the tide, and only 2D, so it
            # stays under 1 GB.
            FieldSnapshotWriter(
                name               = :surface,
                output_file        = "snapshots_surface.jld2",
                variables          = (:η,),
                interval           = 1hour,
                overwrite_existing = true,
            ),
            CheckpointWriter(interval = 12hours, cleanup = true),
        ),
        callbacks = (ProgressCallback(name = :progress, interval = 1hour, report = progress),),
        # The advective CFL limits Δt here, not `max_time_step`: Oslofjord's validation run took
        # 25-32 s steps at 193 m, so this grid should take roughly 2.5 times that.
        time_stepping = AdaptiveTimeStep(
            initial_time_step    = 1second,
            cfl                  = 0.3,
            max_time_step        = 3minutes,
            max_time_step_change = 1.01,
        ),
        # The NorKyst daily mean of `start_date`: T, S, u and v.
        initial_conditions = FromInitialState(),
        # 1 January 2024 to 1 October 2025: 366 + 273 days. NorKyst ends on 2025-10-05, so this is
        # the last whole month the boundary data covers.
        start_date         = DateTime(2024, 1, 1),
        stop_time          = 639days,
        loops              = 1,
        pickup             = false,
    ),
)
