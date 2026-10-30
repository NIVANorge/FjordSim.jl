# The Drammensfjord setup: Geonorge Sjøkart Dybdedata bathymetry, a NorKyst-800m initial state, NVE
# rivers, hourly NorKyst-800m data on the open southern edge and NORA3 atmosphere, on a 170 x 256 x 18 grid covering 10.20-10.50°E, 59.53-59.76°N — 99 m
# cells, and the whole fjord from the Drammenselva and Lierelva mouths down past the Svelvik sill to
# the last narrows before Breiangen. It runs the calendar year 2020, as `oslofjorden.jl` does.
#
# The southern edge is at 59.53°N and the domain ends at 10.50°E because that is the southernmost cut
# whose *eastern* wall is still land: below 59.52°N Breiangen reaches past 10.50°E, so a closed
# eastern wall would sit in open water across the mouth.
#
# Every file goes to this setup's own data root, `~/FjordSim_data/drammensfjorden/`, the Geonorge
# FileGDB (4.5 GB, downloaded by `prepare_bathymetry` if it is missing) included.
#
# # Preparing and running it
#
#   julia --project -m FjordSim prepare_bathymetry     --config examples/drammensfjorden.jl
#   julia --project -m FjordSim download_initial_state --config examples/drammensfjorden.jl
#   julia --project -m FjordSim prepare_initial_state  --config examples/drammensfjorden.jl
#   julia --project -m FjordSim add_rivers             --config examples/drammensfjorden.jl
#   julia --project -m FjordSim download_boundary_data --config examples/drammensfjorden.jl
#   julia --project -m FjordSim prepare_boundary_data  --config examples/drammensfjorden.jl
#   julia --project -m FjordSim download_atmosphere    --config examples/drammensfjorden.jl
#   julia --project -m FjordSim prepare_atmosphere     --config examples/drammensfjorden.jl
#   julia --project -m FjordSim run_simulation         --config examples/drammensfjorden.jl   # needs a GPU

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

data_root = fjord_data_root("drammensfjorden")
# Stated here and passed to the grid and to every model component below, rather than set on
# `Oceananigans.defaults.FloatType` — see `EvenGrid`'s `float_type`.
FT = Float32

FjordConfig(
    # Overloads LatitudeLongitudeGrid(architecture, config) — the one grid hook.
    grid_config = EvenGrid(
        # 99.2 m in longitude, 100.0 m in latitude. The northern edge is 59.76° rather than
        # 59.75° so Lierelva's mouth at 59.7502°N is inside the domain and the northern wall is
        # entirely land.
        size      = (170, 256, 18),
        halo      = (7, 7, 7),
        longitude = (10.20, 10.50),
        latitude  = (59.53, 59.76),
        # A geometric stretch, ratio 1.33 tapering to ~13 m, from a 1 m surface cell to a
        # deepest face at -130 m — the inner basin reaches 124 m. The 11 levels ending at
        # -100 m in 25 m steps that this replaced left the basin floor outside the grid and,
        # where it did not, held it in one 25 m cell with no vertical structure.
        z_faces   = [
            -130.0, -116.6, -103.4, -90.4, -77.6, -65.0, -53.0, -42.2, -33.2,
            -25.8, -19.8, -15.0, -11.2, -8.2, -5.8, -3.8, -2.2, -1.0, 0.0,
        ],
        float_type = FT,
    ),
    # Overloads bathymetry_dataset (required), regrid_options and smoothing_options (both
    # optional) — the hooks prepare_bathymetry dispatches on.
    bathymetry_config = DybdedataConfig(
        data_root             = data_root,
        output_file           = "bathymetry.nc",
        plot_file             = "bathymetry.png",
        raw_resolution_factor = 2,
        padding_cells         = 2,
        include_contours      = false,
        contour_stride        = 10,
        interpolation_passes  = 1,
        # One basin, which is what drops the lakes: Geonorge's soundings include freshwater
        # bodies, and this domain holds several. It also drops Sandebukta, which enters the
        # fjord from the south and is therefore disconnected at this cut.
        major_basins          = 1,
        # The Oslofjord smoothing set, which this setup carried none of; `DybdedataConfig`
        # documents the failure each stage prevents. `open_boundary_land_cells` is 2 rather
        # than Oslofjord's 5 because half of this southern row is land, so the stage has far
        # more to flood here, and every cell it floods lacks an exterior profile of its own.
        minimum_depth         = 2.0,
        open_boundary_land_cells = 2,
        max_island_cells      = 6,
        close_narrow_passages = true,
        spike_ratio           = 0.5,
        # Must equal the `minimum_fractional_cell_height` `Grids.jl` gives `PartialCellBottom`.
        minimum_cell_fraction = 0.2,
        # Oslofjord's limit and sill restoration; see there for what each costs and buys. Here
        # the native grid already resolves Svelvik, so `preserve_sills` holds it at 11.2 m.
        max_slope_factor      = 0.5,
        preserve_sills        = true,
        # Where an open-edge column reaches a level the column inside it does not, the radiation
        # condition stops acting on that face.
        raise_boundary_pockets = true,
        geonorge_cache        = true,
        regrid_cache          = false,
    ),
    # One NorKyst daily mean, the day `start_date` falls in, for `FromInitialState()` below. Its
    # own directory: the download is subset to this domain's box.
    initial_state_config = NorKystInitialStateConfig(
        data_root        = data_root,
        output_directory = "norkyst_initial",
        output_file      = "initial_state.nc",
        parameters       = ["temperature", "salinity", "u_eastward", "v_northward"],
    ),
    # Overloads river_locations, river_series, river_discharge and fetch_rivers (all
    # required). Its own `data_root`: the ELVIS and REGINE caches record the domain they
    # were fetched for under a fixed name, so sharing Oslofjord's `nve` directory would
    # overwrite them.
    river_config = NVERiversConfig(
        data_root         = data_root,
        output_file       = "rivers_nve.nc",
        plot_file         = "rivers_nve.png",
        years             = [2020],
        # Discovery at Oslofjord's threshold: every mouth NVE has in this domain carrying
        # 0.5 m³/s or more, sized by its own REGINE catchment. It replaced the OF800
        # dataset's outlet table, which landed one river in the domain, 602 m from the
        # nearest water. No river coordinate appears in this file.
        minimum_discharge = 0.5,
        # The vertical Courant bound on the 1 m top cell, as on `oslofjorden`: 1 mm/s is
        # w·Δt/Δz = 0.18 at the 3-minute `max_time_step` below. On 99 m cells that spreads
        # Drammenselva's flood over many more columns than Oslofjord's grid needs; the log
        # of `add_rivers` states how many.
        maximum_inflow_velocity = 1.0e-3,
        # Overrides keyed by the mouth's terminal ELVIS `vassdragsnr`, both attaching an
        # observed gauge to a mouth discovery finds anyway. They need `NVE_API_KEY` (free
        # at https://hydapi.nve.no/Users); the map services do not.
        outlets = [
            # Mjøndalen bru carries discharge and temperature on one id, ~5 km above the
            # fjord.
            NVERiver(
                vassdragsnr = "012.A2", name = "Drammenselva",
                discharge_station = "12.534.0", temperature_station = "12.534.0",
            ),
            # Lierelva at Oppsal, 223 km² of a 310 km² catchment. No temperature series,
            # so its water enters at the temperature of the fjord it enters.
            NVERiver(vassdragsnr = "011.A0", discharge_station = "11.6.0"),
        ],
    ),
    # The hourly exterior state along the open southern edge, and the one place that edge is
    # stated. Its own download directory: the band is a thin strip along *this* grid's southern
    # edge, 60 km north of Oslofjord's, so it is different data rather than a subset.
    boundary_data_config = NorKystBoundaryDataConfig(
        data_root        = data_root,
        output_directory = "norkyst_hourly",
        output_file      = "boundaries.nc",
        plot_file        = "boundaries.png",
        open_edges       = :south,
        margin           = 0.05,
        architecture     = :auto,
        parameters       = [
            "temperature", "salinity", "u_eastward", "v_northward", "zeta", "ubar", "vbar",
        ],
        years            = [2020],
    ),
    # Overloads atmosphere_time_steps, atmosphere_source_grid, atmosphere_variable_names,
    # download_atmosphere, prescribed_atmosphere and prescribed_radiation.
    atmosphere_config = NORA3Config(
        data_root        = data_root,
        output_directory = "nora3",
        output_file      = "atmosphere.nc",
        plot_file        = "atmosphere.png",
        resolution       = 0.02,
        padding          = 0.1,
        years            = [2020],
    ),
    # Overloads coupled_model and model_tracers.
    model_config = CoupledHydrostaticModel(
        buoyancy           = SeawaterBuoyancy(FT, equation_of_state = TEOS10EquationOfState(FT)),
        # A `BoundarySponge` for the same reason `oslofjorden.jl` has one: nothing else damps
        # what the open boundary radiates. 13 m² s⁻¹ rather than Oslofjord's 30, because
        # explicit `Δt ≤ Δx²/4ν` is 82 s at 30 on a 99 m cell — below `max_time_step`, so
        # the sponge and not the CFL would silently set the time step. At 13 it is 189 s.
        # `width_cells` stays 16: what the sponge damps is grid-scale noise, so its width
        # is a count of wavelengths rather than a distance.
        closure            = BoundarySponge(
            base = (
                # CATKE takes `w★ = sqrt(max(minimum_tke, e))`, so this floor, not the
                # prognostic TKE, sets κ = 0.098·e_min/N over most of the column.
                # `oslofjorden.jl` has moved to 7e-7, measured there; whether this setup follows
                # is open (TODO.md, "Vertical mixing is too strong").
                CATKEVerticalDiffusivity(FT; minimum_tke = 7e-6),
                # A biharmonic coefficient is only meaningful against Δx⁴, and this grid's
                # 99 m cell is a factor 14 smaller in Δx⁴ than Oslofjord's 193 m one. 1e3
                # damps the 2Δx mode with a 101 min e-folding (`Δx⁴/16ν₄`) against
                # Oslofjord's 72 min at 2e4, and leaves a 3026 s explicit-stability limit
                # `Δt ≤ Δx⁴/32ν₄` that `AdaptiveTimeStep` never measures. The 1e5 this file
                # once copied put that limit at 24 s, below the steps the run takes.
                HorizontalScalarBiharmonicDiffusivity(FT; ν = 1e3, κ = 1e2),
            ),
            width_cells = 16,
            viscosity   = 13.0,
            diffusivity = 6.5,
        ),
        # T and S upwinded horizontally around the river mouths, and salinity's WENO bounded, as in
        # `oslofjorden.jl`, where both are argued. Oceananigans gives any tracer a NamedTuple omits
        # the unbounded `Centered()` default, so it names CATKE's `e`.
        tracer_advection   = (
            T = RiverMouthAdvection(interior = WENO(FT), width_cells = 3),
            S = RiverMouthAdvection(interior = WENO(FT; bounds = (0, 40)), width_cells = 3),
            e = WENO(FT),
        ),
        momentum_advection = WENOVectorInvariant(FT),
        tracers            = (:T, :S),
        coriolis           = HydrostaticSphericalCoriolis(FT),
        sea_ice            = FreezingLimitedOceanTemperature(FT),
        # Overloads free_surface(config, grid), called from inside coupled_model once the
        # grid exists.
        free_surface       = SplitExplicitFreeSurfaceConfig(cfl = 0.7),
        # Anything else the three constructors coupled_model calls accept, one slot each:
        # :ocean_model, :ocean_simulation, :coupled_model.
        extra_kwargs       = (;),
    ),
    # Each piece overloads boundary_condition_sides, so any one can be swapped alone; the
    # model config merges them. Dropping the last would close the domain; tuple order is merge
    # precedence. The two timescales are
    # Marchesiello et al. (2001), and 3 hours rather than a day because at the ~10 s steps
    # these runs take a one-day timescale is ~200x weaker than advection into the boundary
    # cell.
    boundary_condition_configs = (
        AirSeaFluxes(),
        QuadraticBottomDrag(coefficient = 0.003),
        OpenLateralBoundaryFromData(
            inflow_timescale  = 3hours,
            outflow_timescale = 360days,
        ),
    ),
    # How the model is run. Writers, callbacks and time stepping each dispatch on their own
    # supertype; the run itself is read through the run config's accessors.
    run_config = RunConfig(
        results_root       = fjord_results_root("drammensfjorden"),
        architecture       = :auto,
        # Both overload attach_writer!. `variables` may name anything `Oceananigans.fields`
        # exposes on the ocean model; dropping the `CheckpointWriter` turns checkpointing off.
        writers = (
            SnapshotWriter(
                name               = :ocean,
                output_file        = "snapshots_ocean.nc",
                variables          = (:T, :S, :u, :v),
                interval           = 3hours,
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
        # The NorKyst daily mean of `start_date`: T, S, u and v.
        initial_conditions = FromInitialState(),
        # The whole of 2020, a leap year, so 366 days from midnight lands on midnight. The
        # prepare steps pad their time axes to reach the ends of this window.
        start_date         = DateTime(2020, 1, 1),
        stop_time          = 366days,
        # One pass. Raise it to spin the basin up on the same year, carrying the state
        # over; each repetition writes its own `_loopNN` file.
        loops              = 1,
        pickup             = false,
    ),
)
