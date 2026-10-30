# Møre og Romsdal: the coast from Stad to Hitra (5.12-9.00°E, 62.00-63.43°N), with Storfjorden,
# Romsdalsfjorden, the Kristiansund fjords, the Møre shelf and the continental slope off it in the
# north-western corner. A NorKyst-800m initial state, hourly NorKyst-800m data on the open western
# and northern edges, NVE rivers and NORA3 atmosphere, run over 2023 with T and S only.
#
# # Preparing and running it
#
#   julia --project -m FjordSim prepare_bathymetry     --config examples/moro.jl
#   julia --project -m FjordSim download_initial_state --config examples/moro.jl
#   julia --project -m FjordSim prepare_initial_state  --config examples/moro.jl
#   julia --project -m FjordSim add_rivers             --config examples/moro.jl
#   julia --project -m FjordSim download_boundary_data --config examples/moro.jl
#   julia --project -m FjordSim prepare_boundary_data  --config examples/moro.jl
#   julia --project -m FjordSim download_atmosphere    --config examples/moro.jl
#   julia --project -m FjordSim prepare_atmosphere     --config examples/moro.jl
#   julia --project -m FjordSim run_simulation         --config examples/moro.jl   # needs a GPU
#
# `add_rivers` needs `NVE_API_KEY` (free at https://hydapi.nve.no/Users). Data goes to
# `~/FjordSim_data/moro/`, the Geonorge FileGDB to its own directory there, `geonorge/`, downloaded
# (4.5 GB) by `prepare_bathymetry` if it is missing. Results go to `~/FjordSim_results/moro/`.
#
# # What follows `examples/oslofjorden_validation.jl`, and what does not
#
# The physics is the validation run's, with every coefficient that depends on grid spacing scaled
# from 193 m to 500 m: bounded WENO for salinity, CATKE's `minimum_tke` at 7e-7, the biharmonic
# tracer diffusivity κ at 0, `oslofjorden.jl`'s top levels, its bathymetry knobs, sponge, bottom
# drag and open-boundary timescales. Each is argued beside its value, as are the three it does not
# take: 10 m levels between 41 and 131 m, a stronger slope limit and a CFL of 0.4. Three things do
# not carry over. The validation run reads the Norkyst-v3 hindcast because its window predates
# 2017, while 2023 is in the operational archive. Its tide correction was fitted at Viker. Its
# initial state is a previous run's.

using FjordSim
using FjordSim.Bathymetry: GEONORGE_DYBDEDATA_GDB
using Oceananigans:
    CATKEVerticalDiffusivity,
    HydrostaticSphericalCoriolis,
    SeawaterBuoyancy,
    WENO,
    WENOVectorInvariant
using Oceananigans.TurbulenceClosures: HorizontalScalarBiharmonicDiffusivity
using Oceananigans.Units: day, days, hour, hours, minutes, second
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using NumericalEarth: FreezingLimitedOceanTemperature
using Dates: DateTime

data_root = fjord_data_root("moro")
FT = Float32
years = [2023]

# `oslofjorden.jl`'s top levels down to 41 m, then 10 m levels to 131 m and five levels thickening
# from 13 to 28.5 m down to 232.5 m. `examples/oslofjorden_validation.jl` splits 41-132 m into
# eighteen 5.06 m levels, which cut the Drammensfjord's salt-budget diffusivity below 41 m by
# 2.6-7.7 times in 60 days (WENO's implicit vertical diffusion goes as |w|Δz). Here those levels set
# the time step instead: flow up the steep slopes drives a single-column w of 2-5 cm/s in the bottom
# cell, nearly unchanged up through the thin levels, and over the first 41 days of 2023 its |w|/Δz
# held Δt at 67-186 s at CFL 0.3. 10 m levels halve that term. The graded levels below 131 m replace
# a jump from 5.06 to 33.5 m at 132 m.
#
# Below 232.5 m, `oslofjorden.jl`'s 33.5 m layers, continued down the continental slope to -1137 m,
# against a deepest cell of 1153 m in `bathymetry.nc`. The 30 columns below the face, by at most
# 16 m, are on the slope within 11 cells of the open edges, inside the sponge, and `PartialCellBottom`
# clamps them to it. Thicker deep layers would cost the fjord basins, which reach 688 m, their
# vertical structure. `grid.Lz` of 1137 m gives a barotropic wave speed of 106 m/s, so the
# split-explicit substep is about 3 s.
#
# 52 levels on 396 x 315 columns is 6.5 M cells.
z_faces = vcat(
    collect(range(-1137.0, -232.5; step = 33.5)),
    [-204.0, -180.0, -160.0, -144.0, -131.0],
    collect(range(-121.0, -41.0; step = 10.0)),
    [-32.0, -25.0, -19.0, -14.5, -10.8, -7.9, -5.5, -3.7, -2.2, -1.0, 0.0],
)

FjordConfig(
    grid_config = EvenGrid(
        # 0.0098° by 0.00454°: 512 m by 505 m at 62.0°N and 488 m by 505 m at 63.43°N, 2.6 times
        # Oslofjord's 193 m. Every coefficient below that depends on grid spacing is scaled from
        # `oslofjorden_validation.jl` by Δx = 500 m, not copied.
        size       = (396, 315, length(z_faces) - 1),
        halo       = (7, 7, 7),
        longitude  = (5.12, 9.0),
        latitude   = (62.00, 63.43),
        # `z_faces` is written into `bathymetry.nc` and `simulation_grid` reads the file, so a
        # change here needs `prepare_bathymetry` re-run.
        z_faces    = z_faces,
        float_type = FT,
    ),
    # `oslofjorden.jl`'s knobs, which the validation run takes unchanged, but for `max_slope_factor`.
    bathymetry_config = DybdedataConfig(
        data_root             = data_root,
        output_file           = "bathymetry.nc",
        plot_file             = "bathymetry.png",
        # The national FileGDB, in its own directory under this setup's data root.
        geodatabase_file      = joinpath("geonorge", GEONORGE_DYBDEDATA_GDB),
        raw_resolution_factor = 4,
        padding_cells         = 2,
        include_contours      = false,
        contour_stride        = 10,
        interpolation_passes  = 1,
        major_basins          = 1,
        minimum_depth         = 2.0,
        # Off, unlike Oslofjord's 5: the northern open edge runs through the Smøla and Hitra
        # skerries, and the stage would flood the islands within the edge band.
        open_boundary_land_cells = 0,
        # 6 cells, 1.5 km² at this resolution. Check that the Smøla and Hitra skerries keep their
        # sounds open.
        max_island_cells      = 6,
        close_narrow_passages = true,
        spike_ratio           = 0.5,
        # Must equal the `minimum_fractional_cell_height` the grid gives `PartialCellBottom`.
        minimum_cell_fraction = 0.2,
        # 0.4 rather than 0.5: the upslope w above the steepest slopes limits the time step. On
        # `oslofjorden_validation` 0.25 took 34 m off a trench where 0.5 took 8 m (FjordSim's
        # `DybdedataConfig`); `preserve_sills` restores the sill depths.
        max_slope_factor      = 0.4,
        preserve_sills        = true,
        # The western and northern open edges cut across the shelf slope, where an edge column can
        # reach a level the column inside it does not, and the radiation condition stops acting on
        # such a face. Run without it, the domain's fastest current reached 2.9-3.8 m/s.
        raise_boundary_pockets = true,
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
    # Salinity 0 (`NVERiversConfig`'s default `constants`) and NVE's temperature.
    river_config = NVERiversConfig(
        data_root   = data_root,
        output_file = "rivers_nve.nc",
        plot_file   = "rivers_nve.png",
        years       = years,
        # Mouths are discovered from NVE's ELVIS network and REGINE catchments at 0.5 m³/s, as on
        # `oslofjorden.jl`, where that kept 99.6 % of the discharge.
        minimum_discharge = 0.5,
        # The vertical Courant bound on the 1 m top cell: 1 mm/s is w·Δt/Δz = 0.18 at the 3-minute
        # `max_time_step`. On 500 m cells most rivers need one column; the log of `add_rivers`
        # states how many each got.
        maximum_inflow_velocity = 1.0e-3,
        # Water temperature only: every river with a HydAPI daily temperature station on its own
        # watercourse, near the mouth, recorded on all 365 days of 2023. The other mouths carry
        # these rivers' daily mean (`nve_fill_temperatures!`); with none, `add_rivers` fails.
        # Surna's is the one above Trollheim kraftstasjon: the plant's outlet station (112.33.0)
        # lacks 18 days of 2023. Discharge is still every river's REGINE normal; attach
        # `discharge_station`s to the large rivers (Driva, Rauma, Surna, Eira, Aura) as
        # `examples/lista.jl` does.
        outlets = [
            NVERiver(vassdragsnr = "109.A0", name = "Driva", temperature_station = "109.44.0"),
            NVERiver(vassdragsnr = "112.A0", name = "Surna", temperature_station = "112.35.0"),
            NVERiver(vassdragsnr = "103.A20", name = "Rauma", temperature_station = "103.42.0"),
            NVERiver(vassdragsnr = "104.A", name = "Eira", temperature_station = "104.27.0"),
            NVERiver(vassdragsnr = "104.2A0", name = "Visa", temperature_station = "104.23.0"),
            NVERiver(vassdragsnr = "105.A", name = "Osenvassdraget", temperature_station = "105.1.0"),
            NVERiver(vassdragsnr = "096.1A", name = "Hareidselva", temperature_station = "96.3.0"),
            NVERiver(vassdragsnr = "094.4A", name = "Austefjordvassdraget", temperature_station = "94.13.0"),
            NVERiver(vassdragsnr = "102.11Z", name = "Hildreelva", temperature_station = "102.1.0"),
            NVERiver(vassdragsnr = "107.61Z", name = "Sagelva", temperature_station = "107.8.0"),
        ],
    ),
    boundary_data_config = NorKystBoundaryDataConfig(
        data_root        = data_root,
        output_directory = "norkyst_hourly",
        output_file      = "boundaries.nc",
        plot_file        = "boundaries.png",
        # The Norwegian Coastal Current enters through the western edge and leaves through the
        # northern one. In `bathymetry.nc` the western edge is 86 % water and the northern one
        # 85 %; the southern and eastern edges are land but for one fjord-head cell each, 7 m
        # and 5 m deep. The bounding box of the two edge bands is the whole domain, so the
        # hourly download is the whole box for a year: roughly six times Lista's.
        open_edges       = (:west, :north),
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
    model_config = CoupledHydrostaticModel(
        buoyancy           = SeawaterBuoyancy(FT, equation_of_state = TEOS10EquationOfState(FT)),
        closure            = BoundarySponge(
            base = (
                # 7e-7 rather than `oslofjorden.jl`'s 7e-6, as in `oslofjorden_validation.jl`,
                # where 7e-6 mixed the columns towards uniform and freshened the sill basins.
                # Over 640 days at 7e-7 the salinity RMSE against 216 casts was 1.73 g/kg and
                # the temperature RMSE 1.22 °C. Unmeasured here: check `e` against the floor.
                CATKEVerticalDiffusivity(FT; minimum_tke = 7e-7),
                # The validation run's ν = 2e4 m⁴ s⁻¹ scaled by (500/193)⁴ = 45, which keeps its
                # damping times: about 72 min for the 2Δx mode and 56 h at 8Δx. The explicit
                # limit Δt ≤ Δx⁴/32ν₄ is 2000 s at the narrowest 488 m cells.
                #
                # κ is 0, as there: κ = 2e3 left the deep-basin salt budgets unchanged over
                # 60 days but overshot at fronts, leaving 169 wet cells below 0 g/kg.
                HorizontalScalarBiharmonicDiffusivity(FT; ν = 9e5, κ = 0),
            ),
            # The validation run's 16 cells, 8 km here against its 3.2 km, which is under 5 % of
            # this box on each side.
            width_cells = 16,
            # Its 30 and 15 m² s⁻¹ scaled by (500/193)² = 6.7. Δt ≤ Δx²/4ν is 300 s at the
            # narrowest cells, above `max_time_step`.
            viscosity   = 200.0,
            diffusivity = 100.0,
        ),
        # Salinity bounded, as in `oslofjorden_validation.jl`, where plain WENO drove salinity
        # below TEOS10's -32 psu floor beside a large river entering shallow columns. A
        # NamedTuple leaves any tracer it omits on `Centered()`, so it names CATKE's `e` as well.
        #
        # T and S upwinded horizontally around the river mouths: WENO across a mouth's fresh front
        # made water saltier than any source (35.43 beside the Driva; none exceeds 35.05) and a
        # 2 m Surna column −1 °C. Three cells is WENO5's stencil, so no face outside reconstructs
        # across a river column; 4.1 % of the wet cells.
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
        # Marchesiello et al. (2001): nudge hard on inflow, radiate on outflow. Oslofjord's
        # timescales, applied to every open edge.
        OpenLateralBoundaryFromData(
            inflow_timescale  = 3hours,
            outflow_timescale = 360days,
        ),
    ),
    run_config = RunConfig(
        results_root       = fjord_results_root("moro"),
        architecture       = :auto,
        writers = (
            # Daily, as in the validation run. One record of these five fields is 147 MB, so the
            # year is 54 GB.
            SnapshotWriter(
                name               = :ocean,
                output_file        = "snapshots_ocean.nc",
                variables          = (:T, :S, :u, :v, :e),
                interval           = 1day,
                overwrite_existing = true,
            ),
            # `η` in JLD2, which the NetCDF writer cannot take (`FieldSnapshotWriter`). Hourly, for
            # the tide.
            FieldSnapshotWriter(
                name               = :surface,
                output_file        = "snapshots_surface.jld2",
                variables          = (:η,),
                interval           = 1hour,
                overwrite_existing = true,
            ),
            # The two permanent Kartverket gauges in the box, from Kartverket's station register
            # (`vannstand.kartverket.no/tideapi.php?tide_request=stationlist&type=perm`). Ten
            # minutes, as in the validation run, for the harmonic analysis `validate_simulation`
            # does against Kartverket's own record.
            FieldStationWriter(
                name               = :tidegauge,
                output_file        = "stations_tidegauge.jld2",
                variables          = (:η,),
                stations           = [
                    Station(name = "Ålesund",      longitude = 6.151946, latitude = 62.469414),
                    Station(name = "Kristiansund", longitude = 7.736140, latitude = 63.113920),
                ],
                interval           = 10minutes,
                overwrite_existing = true,
            ),
            CheckpointWriter(interval = 12hours, cleanup = true),
        ),
        callbacks = (ProgressCallback(name = :progress, interval = 1hour, report = progress),),
        # The advective CFL limits Δt, not `max_time_step`. At `oslofjorden.jl`'s 0.3 the first
        # 41 days of 2023 took 67-186 s steps. At 0.5, with Δt near 90 s, the top cells on the open
        # western edge went unstable on 12 February 2023 at day 42.3, under a 0.8-0.9 m/s
        # NorKyst surface inflow at low tide.
        time_stepping = AdaptiveTimeStep(
            initial_time_step    = 1second,
            cfl                  = 0.4,
            max_time_step        = 3minutes,
            max_time_step_change = 1.01,
        ),
        # The NorKyst daily mean of `start_date` for T, S, u and v.
        initial_conditions = FromInitialState(),
        # The calendar year 2023.
        start_date         = DateTime(2023, 1, 1),
        stop_time          = 365days,
        loops              = 1,
        pickup             = false,
    ),
)
