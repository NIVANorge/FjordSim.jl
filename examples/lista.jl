# Lista: the open coast of southernmost Norway, from Sokndal past Lindesnes (6.22-7.25°E,
# 57.89-58.34°N), with Fedafjorden, the Flekkefjord sounds, Rosfjorden and the northern slope of the
# Norwegian Trench. NorKyst-800m forcing and an open boundary on three sides, NVE rivers, NORA3
# atmosphere, run from 1 January 2024 to 1 October 2025.
#
# # Preparing and running it locally
#
#   julia --project -m FjordSim prepare_bathymetry  --config examples/lista.jl
#   julia --project -m FjordSim download_forcing    --config examples/lista.jl
#   julia --project -m FjordSim prepare_forcing     --config examples/lista.jl
#   julia --project -m FjordSim add_rivers          --config examples/lista.jl
#   julia --project -m FjordSim download_boundaries --config examples/lista.jl
#   julia --project -m FjordSim prepare_boundaries  --config examples/lista.jl
#   julia --project -m FjordSim download_atmosphere --config examples/lista.jl
#   julia --project -m FjordSim prepare_atmosphere  --config examples/lista.jl
#   julia --project -m FjordSim run_simulation      --config examples/lista.jl   # needs a GPU
#
# `add_rivers` needs `NVE_API_KEY` (free at https://hydapi.nve.no/Users). `prepare_bathymetry` reads
# Oslofjord's extracted copy of the national Geonorge FileGDB, `~/FjordSim_data/oslofjorden/`, by
# absolute path, and downloads it there (4.5 GB) if it is missing. Everything else goes to
# `~/FjordSim_data/lista/`, and results to `~/FjordSim_results/lista/`.
#
# # On GCP
#
# Bathymetry locally, because the FileGDB is already on this machine and is not in the bucket. The
# downloads and regrids then run on the VM, which opens OPeNDAP about five times faster:
#
#   julia --project -m FjordSim prepare_bathymetry --config examples/lista.jl
#   gcp/fjordsim-gcp push-data lista
#   gcp/fjordsim-gcp up && gcp/fjordsim-gcp bootstrap
#   gcp/fjordsim-gcp run --config examples/lista.jl --steps download_forcing,prepare_forcing,add_rivers,download_boundaries,prepare_boundaries,download_atmosphere,prepare_atmosphere
#   gcp/fjordsim-gcp run --config examples/lista.jl --steps run_simulation --gpu --run-id lista-2024
#   gcp/fjordsim-gcp logs lista-2024
#   gcp/fjordsim-gcp pull-results lista
#   gcp/fjordsim-gcp down
#
# `run` takes `NVE_API_KEY` from your environment. The VM is Spot and is deleted after 24 hours, so
# a long run usually has to be resumed from its last checkpoint:
#
#   gcp/fjordsim-gcp up && gcp/fjordsim-gcp bootstrap
#   gcp/fjordsim-gcp run --config examples/lista.jl --steps run_simulation --gpu --run-id lista-2024 --resume
#
# The file's basename, `lista`, is what names the staged directories in the bucket, so keep it.
#
# # Settle these after `prepare_bathymetry`, before the rest
#
# Both were estimated from EMODnet's bathymetry (`rest.emodnet-bathymetry.eu`), not from the
# Geonorge soundings the model uses. If either changes, re-run every prepare step after
# `prepare_bathymetry`.
#
# 1. **`z_faces`.** EMODnet's deepest point in the box is 444 m, on the southern edge. The deepest
#    face, -459.5 m, has to clear the deepest depth in `bathymetry.nc`.
# 2. **`open_edges`.** EMODnet has the southern edge wet all the way across (300-444 m), the western
#    edge wet for 93 % of its length (to 352 m) and the eastern edge for its southern 31 % (to 357 m),
#    and the northern edge dry. Check them against the land mask in `bathymetry.png`.

using FjordSim
using Oceananigans
using Oceananigans.Units
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using NumericalEarth: FreezingLimitedOceanTemperature
using Dates: DateTime

data_root = fjord_data_root("lista")
FT = Float32
# NorKyst-800m, daily and hourly alike, ends on 2025-10-05 on thredds.met.no, so 2025 is the last
# year either forcing source can reach.
years = [2024, 2025]

FjordConfig(
    grid_config = EvenGrid(
        # 0.00837° by 0.0045°: about 490 m by 500 m at 58.1°N, 2.5 times Oslofjord's 193 m. Every
        # coefficient below that depends on grid spacing (Δx or Δx⁴) is scaled from `oslofjorden()`
        # by that factor, not copied.
        size      = (123, 100, 24),
        halo      = (7, 7, 7),
        longitude = (6.22, 7.25),
        latitude  = (57.89, 58.34),
        # `oslofjorden()`'s upper 16 layers unchanged, from a 1 m surface cell stretching by 1.25,
        # down to -132 m. Below that, one 33.5 m layer and then seven of 42 m: each step up in
        # thickness is still a ratio of about 1.25, and the count stays at 24. That reaches -459.5 m, deep enough for the trench
        # (gate 1 in the header).
        z_faces   = [
            -459.5, -417.5, -375.5, -333.5, -291.5, -249.5, -207.5, -165.5,
            -132.0, -105.0, -83.0, -66.0, -52.0, -41.0, -32.0, -25.0,
            -19.0, -14.5, -10.8, -7.9, -5.5, -3.7, -2.2, -1.0, 0.0,
        ],
        float_type = FT,
    ),
    bathymetry_config = DybdedataConfig(
        data_root             = data_root,
        output_file           = "bathymetry.nc",
        plot_file             = "bathymetry.png",
        # Shared with `oslofjorden()`: it is one national file, the same for any Norwegian domain.
        geodatabase_file      = geodatabase_path(oslofjorden().bathymetry_config),
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
        # Oslofjord's limit, where it was measured to cost almost nothing. Here it acts on the
        # trench slope, which drops from about 100 m to 300 m within a few cells. Compare
        # `bathymetry.png` with the raw soundings before trusting it.
        max_slope_factor      = 0.25,
        geonorge_cache        = true,
        regrid_cache          = false,
    ),
    forcing_config = NorKystConfig(
        data_root        = data_root,
        output_directory = "norkyst",
        output_file      = "forcing.nc",
        plot_file        = "forcing.png",
        architecture     = :auto,
        parameters       = ["temperature", "salinity", "u_eastward", "v_northward"],
        years            = years,
        rivers           = NVERiversConfig(
            data_root   = data_root,
            output_file = "forcing_rivers_nve.nc",
            plot_file   = "forcing_rivers_nve.png",
            years       = years,
            # Mouths are discovered from NVE's ELVIS network and REGINE catchments. At 0.5 m³/s
            # that finds 20 of them, carrying 310.8 of the 315.0 m³/s that the 35 mouths above
            # 0.15 m³/s carry. Sira (130.5), Kvina (89.5) and Lygna (39.6) dominate.
            minimum_discharge = 0.5,
            # A 5 m surface plume: four cells, to -5.5 m, on this grid.
            default_plume_depth = 5.0,
            minimum_levels = 0,
            # A floor on 1/λ. λΔt stays at or below 0.3 even at the 3-minute `max_time_step`.
            minimum_relaxation_timescale = 600.0,
            # Overrides, keyed by each mouth's terminal ELVIS `vassdragsnr` (from the discovery).
            # Only the *mean* of a discharge series is used (λ = Q̄/V), so a gauge's job here is
            # Q̄. Every series named here was checked for daily values over the run window.
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
                    # weaker. On Oslofjord that was 85 % of wet cells. Measured from NorKyst columns
                    # in the trench here (February and August 2024), 7e-6 gives 4e-5 m² s⁻¹ in the
                    # summer pycnocline and 1-3e-4 below 100 m. Oceananigans' default of 1e-9 gives
                    # 1e-8, below molecular heat diffusion, so it means no background mixing at all.
                    # The open part of the box is flushed by the coastal current within days, so
                    # there the floor hardly matters. It matters in the sill basins of Fedafjorden
                    # and Rosfjorden, which is why Oslofjord's sill-fjord value is kept. Unmeasured
                    # here: check `e` against the floor in the snapshots.
                    CATKEVerticalDiffusivity(FT; minimum_tke = 7e-6),
                    # Oslofjord's 2e4 m⁴ s⁻¹ scaled by (490/193)⁴. That keeps the same damping
                    # times: about 74 min for the 2Δx mode and about 52 h for an 8Δx (4 km) eddy. The
                    # explicit limit Δt ≤ Δx⁴/32ν₄ is 2200 s.
                    HorizontalScalarBiharmonicDiffusivity(FT; ν = 8e5, κ = 8e4),
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
        # The NorKyst state at `start_date`. `prepare_forcing` pads the daily 12:00 records so
        # that one lands exactly at midnight on 1 January.
        initial_conditions = FromForcing(),
        # 1 January 2024 to 1 October 2025: 366 + 273 days. NorKyst ends on 2025-10-05, so this is
        # the last whole month either forcing source covers.
        start_date         = DateTime(2024, 1, 1),
        stop_time          = 639days,
        loops              = 1,
        pickup             = false,
    ),
)
