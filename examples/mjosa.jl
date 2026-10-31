# The Mjøsa setup: a closed freshwater lake, with NVE lake depth-map bathymetry, NVE rivers — the
# Vorma outflow closing the water budget on the observed lake level — and NORA3 atmosphere, on a
# 266 x 456 x 38 grid of 200 m cells covering 10.36-11.34°E, 60.36-61.18°N. It runs the calendar
# year 2025 from a lake at rest.
#
# z = 0 is the lake surface, 123 m above sea level and regulated between about 120 and 123 m; the
# depth map's depths are below it. No edge is open, and the water is fresh (S = 0 throughout).
# There is no ice model: `FreezingLimitedOceanTemperature` keeps the water at or above 0 °C but does
# not insulate it, so winter heat loss and wind mixing are too strong wherever the lake freezes.
#
# # Preparing and running it
#
#   julia --project -m FjordSim prepare_bathymetry  --config examples/mjosa.jl
#   julia --project -m FjordSim add_rivers          --config examples/mjosa.jl   # needs NVE_API_KEY
#   julia --project -m FjordSim download_atmosphere --config examples/mjosa.jl
#   julia --project -m FjordSim prepare_atmosphere  --config examples/mjosa.jl
#   julia --project -m FjordSim run_simulation      --config examples/mjosa.jl   # needs a GPU

using FjordSim
using Oceananigans:
    CATKEVerticalDiffusivity,
    HydrostaticSphericalCoriolis,
    SeawaterBuoyancy,
    WENO,
    WENOVectorInvariant
using Oceananigans.TurbulenceClosures: HorizontalScalarBiharmonicDiffusivity
using Oceananigans.Units: days, hour, hours, minutes, second
using Dates: DateTime
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using NumericalEarth: FreezingLimitedOceanTemperature

data_root = fjord_data_root("mjosa")
FT = Float32

# Every gauge below is scaled by one factor, so the gauged rivers carry the whole 16 186 km² of
# land draining to the lake (16 555 km² catchment less the 369 km² lake) rather than the 13 079 km²
# their six gauges measure: Losna 11 205, Aulestad 870, Moelva 170, Vismunda 192, Lena 184 and
# Svartelva 458. The ungauged 19 % then enters at the gauged mouths. The outlet closes the budget
# whatever the inflows are, so this sets the lake's through-flow, not its level: scaled, the 2025
# inflow is ~364 m³/s against ~344 m³/s at the Vorma gauge (Ertesekken, scaled to the lake's
# catchment).
inflow_scale = 16_186 / 13_079

FjordConfig(
    grid_config = EvenGrid(
        # 200 m in both directions. The box holds the lake's whole outline (10.39-11.31°E,
        # 60.37-61.17°N) with a margin of land on every side, so every edge is closed.
        size      = (266, 456, 38),
        halo      = (7, 7, 7),
        longitude = (10.36, 11.34),
        latitude  = (60.36, 61.18),
        # A geometric stretch from a 1 m surface cell, ratio 1.15, to 18.8 m at 137.6 m: 12
        # levels in the top 20 m and 2-4 m cells through the summer thermocline at 10-25 m.
        # Below, 16 even 19.1 m layers to the deepest cell `prepare_bathymetry` writes, 443.3 m.
        z_faces   = [
            -443.3, -424.2, -405.1, -386.0, -366.9, -347.8, -328.7, -309.6, -290.5, -271.3,
            -252.2, -233.1, -214.0, -194.9, -175.8, -156.7, -137.6, -118.8, -102.4, -88.2,
            -75.8, -65.1, -55.7, -47.6, -40.5, -34.4, -29.0, -24.3, -20.3, -16.8, -13.7,
            -11.1, -8.8, -6.7, -5.0, -3.5, -2.1, -1.0, 0.0,
        ],
        float_type = FT,
    ),
    # NVE's depth map of Mjøsa, lake number 118: contours every 20 m to 440 m and the 453 m
    # deepest point. Geonorge's sea-chart soundings have none in the lake.
    bathymetry_config = NVELakeBathymetryConfig(
        data_root             = data_root,
        output_file           = "bathymetry.nc",
        plot_file             = "bathymetry.png",
        lake_number           = 118,
        raw_resolution_factor = 2,
        major_basins          = 1,
        # The Oslofjord smoothing set, less its open-edge stages.
        minimum_depth         = 2.0,
        # Islets up to 0.24 km²; Helgøya, at 18 km² or ~450 cells, stays.
        max_island_cells      = 6,
        close_narrow_passages = true,
        spike_ratio           = 0.5,
        # Must equal the `minimum_fractional_cell_height` `Grids.jl` gives `PartialCellBottom`.
        minimum_cell_fraction = 0.2,
        max_slope_factor      = 0.5,
        regrid_cache          = false,
    ),
    # Stated outlets: NVE's outlet discovery finds only mouths draining to the sea. The mouths are
    # where each river's ELVIS segment meets the lake's own segments; the gauges carry 2025 daily
    # discharge.
    river_config = NVERiversConfig(
        data_root               = data_root,
        output_file             = "rivers_nve.nc",
        plot_file               = "rivers_nve.png",
        years                   = [2025],
        # 1 mm/s through the 1 m top cell is a vertical Courant number of 0.18 at the
        # 3-minute `max_time_step`; Lågen and Gausa's ~1100 m³/s scaled peak spreads over ~27
        # columns.
        maximum_inflow_velocity = 1.0e-3,
        outlets = [
            # Gudbrandsdalslågen at Lillehammer, with Gausa, which joins it above the mouth. Lågen
            # has no water temperature for 2025 after August, so both take Glomma's at
            # Glomstadfoss, the nearest large inland river with a full year.
            NVERiver(
                id = 1, name = "Gudbrandsdalslågen", longitude = 10.40035, latitude = 61.13979,
                discharge_station = "2.145.0", temperature_station = "2.43.0",
                discharge_fraction = inflow_scale,
            ),
            NVERiver(
                id = 2, name = "Gausa", longitude = 10.40035, latitude = 61.13979,
                discharge_station = "2.28.0", temperature_station = "2.43.0",
                discharge_fraction = inflow_scale,
            ),
            NVERiver(
                id = 3, name = "Moelva", longitude = 10.69257, latitude = 60.92273,
                discharge_station = "2.1147.0", discharge_fraction = inflow_scale,
            ),
            NVERiver(
                id = 4, name = "Vismunda", longitude = 10.6063, latitude = 60.96825,
                discharge_station = "2.463.0", temperature_station = "2.463.0",
                discharge_fraction = inflow_scale,
            ),
            NVERiver(
                id = 5, name = "Lenaelva", longitude = 10.9646, latitude = 60.65964,
                discharge_station = "2.634.0", discharge_fraction = inflow_scale,
            ),
            NVERiver(
                id = 6, name = "Svartelva", longitude = 11.14499, latitude = 60.78803,
                discharge_station = "2.464.0", temperature_station = "2.464.0",
                discharge_fraction = inflow_scale,
            ),
            # Vorma at Minnesund removes what the rivers bring in less what the lake stores, on
            # the level at Hamar and the lake's 369.28 km², so the model level follows the
            # observed one (120.3-123.0 m in 2025).
            NVERiver(
                id = 7, name = "Vorma", longitude = 11.23633, latitude = 60.39769,
                level_station = "2.101.0", level_area = 369.28e6,
            ),
        ],
    ),
    atmosphere_config = NORA3Config(
        data_root        = data_root,
        output_directory = "nora3",
        output_file      = "atmosphere.nc",
        plot_file        = "atmosphere.png",
        resolution       = 0.02,
        padding          = 0.1,
        years            = [2025],
    ),
    model_config = CoupledHydrostaticModel(
        # TEOS-10 at S = 0 has fresh water's density maximum near 4 °C.
        buoyancy           = SeawaterBuoyancy(FT, equation_of_state = TEOS10EquationOfState(FT)),
        # No `BoundarySponge`: no edge is open.
        closure            = (
            # Oslofjorden's measured background-mixing floor (`oslofjorden.jl`).
            CATKEVerticalDiffusivity(FT; minimum_tke = 7e-7),
            # Oslofjorden's 2e4 on 193 m cells scaled by (200/193)⁴: the 2Δx mode e-folds in
            # Δx⁴/16ν₄ = 72 min, and the explicit limit Δx⁴/32ν₄ = 2170 s is far above
            # `max_time_step`.
            HorizontalScalarBiharmonicDiffusivity(FT; ν = 2.3e4, κ = 0),
        ),
        # Upwinded around the river mouths, as in `oslofjorden.jl`; Oceananigans gives a tracer
        # the NamedTuple omits an unbounded `Centered()`, so `e` is named.
        tracer_advection   = (
            T = RiverMouthAdvection(interior = WENO(FT), width_cells = 3),
            S = RiverMouthAdvection(interior = WENO(FT; bounds = (0, 40)), width_cells = 3),
            e = WENO(FT),
        ),
        momentum_advection = WENOVectorInvariant(FT),
        tracers            = (:T, :S),
        coriolis           = HydrostaticSphericalCoriolis(FT),
        sea_ice            = FreezingLimitedOceanTemperature(FT),
        free_surface       = SplitExplicitFreeSurfaceConfig(cfl = 0.7),
        extra_kwargs       = (;),
    ),
    boundary_condition_configs = (
        AirSeaFluxes(),
        QuadraticBottomDrag(coefficient = 0.003),
    ),
    run_config = RunConfig(
        results_root       = fjord_results_root("mjosa"),
        architecture       = :auto,
        writers = (
            SnapshotWriter(
                name               = :ocean,
                output_file        = "snapshots_ocean.nc",
                variables          = (:T, :S, :u, :v),
                interval           = 3hours,
                overwrite_existing = true,
            ),
            FieldSnapshotWriter(
                name               = :surface,
                output_file        = "snapshots_surface.jld2",
                variables          = (:η,),
                interval           = 3hours,
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
        # At rest and uniform: on 1 January a deep dimictic lake is near-isothermal, its deep
        # water at 3.5-4 °C, and no open source has a lake-wide profile for the day.
        initial_conditions = (T = 3.5, S = 0.0),
        start_date         = DateTime(2025, 1, 1),
        stop_time          = 365days,
        loops              = 1,
        pickup             = false,
    ),
)
