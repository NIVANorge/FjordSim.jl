module Configs

export FjordConfig

"""
    AbstractGridConfig

Supertype for grid configurations. A concrete subtype describes how to build the simulation
grid.

# Methods a subtype provides
- `domain_grid(config, architecture)`: the bare grid, with no bathymetry. Required. Every
  `prepare_*` pipeline builds its target geometry with this.
- `simulation_grid(config, bathymetry_file, architecture)`: the grid the simulation runs on, with
  the processed bathymetry immersed into it. Required. `build_simulation` and every step that
  regrids onto the model grid read it back through this.

Both are FjordSim generics rather than methods on `Oceananigans.LatitudeLongitudeGrid` and
`Oceananigans.ImmersedBoundaries.ImmersedBoundaryGrid`, which is what the grid hook used to be.
Those names state a *result* type, so a config describing a rectilinear, single-column or
curvilinear grid could not implement the hook under either name without lying about what it
returns. These two name the *role* instead, leaving the return type to the subtype.

`FjordSim.Grids.EvenGrid` is the built-in implementation, and returns a `LatitudeLongitudeGrid`
and an `ImmersedBoundaryGrid` wrapping one.
"""
abstract type AbstractGridConfig end

"""
    AbstractBathymetryConfig

Supertype for bathymetry configurations. A concrete subtype describes one bathymetry source
and the pipeline that turns it into a processed FjordSim bathymetry NetCDF.

# Fields a subtype provides
- `data_root`, `output_file`, `plot_file`: resolved by `bathymetry_path` and `plot_path`.

# Methods a subtype provides
- `bathymetry_dataset(target_grid, config)`: the NumericalEarth dataset to regrid from.
  Required; `prepare_bathymetry` calls it and everything after it is source-agnostic.
- `regrid_options(config)`: named tuple forwarded to `NumericalEarth.regrid_bathymetry`.
  Optional, defaults to `(;)`.
- `smoothing_options(config)`: named tuple forwarded to `smooth_bathymetry_gaps!`, which is how a
  source opts into the six post-regrid stages — `max_island_cells`, `close_narrow_passages`,
  `spike_ratio`, `minimum_cell_fraction`, `max_slope_factor` (plus `minimum_depth`, which the
  slope limiter must not undo) and `preserve_sills`. Optional, defaults to `(;)`, leaving only the
  topological cleanup every source gets.

`FjordSim.Bathymetry.DybdedataConfig` is the built-in implementation, for Geonorge Sjøkart
Dybdedata; `src/Bathymetry/geonorge.jl` is the template to copy for a new source.
"""
abstract type AbstractBathymetryConfig end

"""
    AbstractForcingConfig

Supertype for forcing configurations. A concrete subtype describes one forcing dataset and
how to download and subset it.

# Fields a subtype provides
- `data_root`, `output_file`, `plot_file`: resolved by `forcing_path` and `plot_path`.
- `output_directory`: resolved by `forcing_directory`; only needed by a dataset that downloads.
- `parameters`: source variable names to prepare.
- `architecture`: `:auto`, `:cpu` or `:gpu`, resolved by `interpolation_architecture` to decide
  where `prepare_forcing` runs its interpolation kernel.

The open edge is *not* a field here. It belongs to the open-boundary data config
(`AbstractBoundaryDataConfig`), which hangs off `FjordConfig` rather than off this config, and
reaches `regrid_forcing` as its `edges` keyword — a forcing dataset says nothing about which edge of
the domain is open. It was a field here while the edge carried an interior relaxation band this
config's `relaxation_timescale` set the rate of; that band is gone.

# Methods a subtype provides
- `forcing_time_steps(config)`: the downloaded time records, as `SourceRecord`s. Required.
- `forcing_source_grid(config, filepath)`: geometry of the source data, e.g. a
  `ProjectedSourceGrid`. Required.
- `forcing_variable_names(config)`: source variable name => FjordSim forcing name. Required.
- `fetch_forcing(target_grid, config)`: fetch the source data. Only if it downloads.
- `simulation_forcing(config, grid, filepath, tracers, reference_date)`: the forcing term object
  `coupled_model` consumes, read from the resolved forcing file at `filepath`, so this takes
  the path rather than resolving it. Optional; defaults to `forcing_from_file`, the FjordSim NetCDF layout
  every built-in source already writes. A source whose prepared files are not that contract
  overrides it.
- `forcing_date_range(config, filepath)`: first and last date the prepared file covers, as a
  `(first, last)` tuple, or `nothing` to skip the check. Optional; defaults to reading `ds["time"]`
  from the prepared NetCDF. The mirror of `atmosphere_date_range`, and a source that overrode
  `simulation_forcing` to read a different layout overrides this too — otherwise
  `validate_time_coverage` would check its coverage with a reader it does not use.

`FjordSim.Forcing.NorKystConfig` is the built-in implementation, for NorKyst-800m;
`src/Forcing/norkyst.jl` is the template to copy for a new dataset.
"""
abstract type AbstractForcingConfig end

"""
    AbstractInitialStateConfig

Supertype for initial-state configurations. A concrete subtype describes one dataset the ocean state
at `start_date` is read from, and how to download it. `prepare_initial_state` regrids that one record
onto the simulation grid, and `FromInitialState` initial conditions read it back.

Independent of the forcing config: the forcing file has to span the whole run window, because the
run reads it as forcing terms at every step, while this file is a single record read once. A setup
wanting a realistic start without interior forcing names this and no forcing config.

The record is the one `start_date` falls in, so changing `start_date` means re-running both steps;
the prepared file records the `start_date` it was prepared for, and `build_simulation` refuses a
mismatch.

# Fields a subtype provides
- `data_root`, `output_file`: resolved by `initial_state_path`, the prepared file.
- `output_directory`: resolved by `initial_state_directory`; the directory the download writes into.
- `parameters`: source variable names to prepare.

# Methods a subtype provides
- `initial_state_record(config, date)`: the downloaded `SourceRecord` holding the state at `date`.
  Required.
- `initial_state_source_grid(config, filepath)`: geometry of the source data, e.g. a
  `ProjectedSourceGrid`. Required.
- `initial_state_variable_names(config)`: source variable name => FjordSim state name. Required.
- `fetch_initial_state(target_grid, config, date)`: fetch the record for `date`. Only if it
  downloads.

`FjordSim.InitialState.NorKystInitialStateConfig` is the built-in implementation, for one day of
NorKyst-800m; `src/InitialState/norkyst_initial_state.jl` is the template to copy for a new dataset.
"""
abstract type AbstractInitialStateConfig end

"""
    AbstractRiverConfig

Supertype for river configurations. A concrete subtype describes one river dataset: where its
rivers reach the sea, how much water each one carries, and what that water is. `add_rivers` turns
it into a river file of its own, which the run reads as a volume inflow through the surface of a
few columns at every mouth.

Rivers are optional: a `FjordConfig` whose `river_config` is `nothing` skips the step entirely.
They are independent of the forcing config — the river file stands alone, so a setup can name rivers
with no interior forcing.

# Fields a subtype provides
- `data_root`, `output_file`, `plot_file`: resolved by `river_forcing_path`, the river file, and by
  `plot_path`.
- `search_radius`: read by the default `river_search_radius`.
- `maximum_inflow_velocity`: m/s, read by the default `river_maximum_inflow_velocity` — how fast a
  river's water may cross the surface of the columns it enters through, which decides how many
  columns a large river is spread over.

# Methods a subtype provides
- `river_locations(config)`: the river outlets, as `RiverLocation`s. Required.
- `river_discharge(config, times)`: a `(river, time)` matrix of discharge in m³/s, one row per
  `river_locations` entry. Required.
- `river_series(config, times)`: FjordSim tracer name => a `(river, time)` matrix of the values the
  river water carries, one row per `river_locations` entry. Every entry of a river the grid carries
  has to be finite — there is no "carries the water it enters", which is unstable — and a model
  tracer it names no matrix for enters at zero. Required.
- `fetch_rivers(target_grid, config)`: fetch the source data. Only if it downloads. Takes the
  bare domain grid, like `fetch_forcing` and `fetch_boundary_data`, because a dataset that has
  to *find* its rivers rather than being told where they are needs the domain bounds; one whose
  outlets are stated in the config ignores it.
- `river_search_radius(config)`: how far to look for a coastal cell, in cells. Optional,
  defaults to `config.search_radius`.
- `river_maximum_inflow_velocity(config)`: in m/s. Optional, defaults to
  `config.maximum_inflow_velocity`.

`FjordSim.Rivers.NVERiversConfig` is the built-in implementation: it discovers river mouths from
NVE's map services and reads observed discharge and water temperature from NVE's HydAPI, which needs
an API key. `src/Rivers/nve_rivers.jl` is the template to copy for a new dataset.
"""
abstract type AbstractRiverConfig end

"""
    AbstractBoundaryDataConfig

Supertype for open-boundary *data* configurations. A concrete subtype describes one dataset of
exterior state along a lateral boundary, and how to download and regrid it into the boundary NetCDF
the open lateral boundary condition reads.

Not to be confused with `AbstractBoundaryConditionConfig`, one word away: that one describes the
boundary *condition* — which scheme acts where, and how fast it nudges — while this one describes
the exterior *values* the condition relaxes towards. One is physics, the other is data.

Open-boundary data is optional: a `FjordConfig` whose `boundary_data_config` is `nothing` skips both
boundary steps entirely, and its lateral boundaries are closed walls — `water_mask` unmasks no
velocity face row, and no boundary condition has an exterior state to nudge towards.

It hangs off `FjordConfig` rather than off the forcing config, and states the open edge itself. The
exterior values were always a separate file read by a separate pipeline; being configured *inside*
the forcing config said otherwise, and meant a setup running with no interior forcing could not have
a data-driven open boundary at all.

The boundary file is separate from the prepared forcing file rather than a part of it, because the
exterior state a Flather/Orlanski boundary needs is **hourly** — a daily mean has the tide averaged
out of it, and elevation with no tide is not worth prescribing. Only the boundary row needs that
cadence, so it is a few hundred MB rather than the hundreds of GB an hourly interior forcing would
be. Every variable in the file is named for its side (`south_T`, `south_eta`, …), so a second open
edge is a set of new variables in the same file rather than a new file or a format change.

# Fields a subtype provides
- `data_root`, `output_file`, `plot_file`: resolved by `boundary_data_path` and `plot_path`.
- `output_directory`: resolved by `boundary_data_directory`; the directory the download writes its
  intermediate files into.
- `open_edges`: which lateral boundaries the domain is open on, each one of `:south`, `:north`,
  `:west` or `:east`, read by `open_edges` as a `Vector{Symbol}`. A domain in the open ocean names all
  four; a fjord names the one it connects to the sea through. Stated once, here:
  `download_boundary_data` and `prepare_boundary_data` read it off this config rather than taking it beside
  it, `water_mask` unmasks each of those edges' velocity face rows for both this pipeline and the
  interior forcing's, and the open lateral boundary condition puts its schemes on every one of them.
- `parameters`: source variable names to prepare.
- `architecture`: `:auto`, `:cpu` or `:gpu`, resolved by `interpolation_architecture`.

# Methods a subtype provides
- `boundary_time_steps(config)`: the downloaded time records, as `SourceRecord`s. Required.
- `boundary_source_grid(config, filepath)`: geometry of the source data, e.g. a
  `ProjectedSourceGrid`. Required.
- `boundary_variable_names(config)`: source variable name => FjordSim boundary name. Required.
- `fetch_boundary_data(target_grid, config)`: fetch the source data, reading the edges it subsets
  along with `open_edges(config)`. Only if it downloads.
- `boundary_date_range(config)`: first and last date the prepared file covers, as a
  `(first, last)` tuple, or `nothing` to skip the check. Optional; defaults to reading `ds["time"]`
  from the prepared NetCDF, so `validate_time_coverage` also guards the `Cyclical()` wrap here.
- `boundary_source_slab(config, reader, step, source_name)`: one source slab as the writer
  interpolates from. Optional; defaults to `blended_slab`, a plain read. Override it only for a
  variable that must be *derived* from more than one of the source's own — which a velocity component
  is, when the source states it along its own grid axes rather than eastward/northward.

`FjordSim.BoundaryData.NorKystBoundaryDataConfig` is the built-in implementation, for MET Norway's hourly
NorKyst-800m collection; `src/Boundaries/BoundaryData/norkyst_boundary_data.jl` is the template to copy
for a new dataset.
"""
abstract type AbstractBoundaryDataConfig end

"""
    AbstractAtmosphereConfig

Supertype for atmosphere configurations. A concrete subtype describes one atmospheric
reanalysis dataset, how to download it and how to regrid it onto the regular longitude/latitude
grid the simulation-time reader consumes.

The prepared file is what `FjordSim.Atmospheres.NORA3.MultiYearNORA3` reads, so its layout is a
fixed contract rather than a dataset detail: variables of shape `(longitude, latitude, time)`,
uniformly spaced 1D `lon` and `lat` cell centers, and a CF-encoded `time`. `prepare_atmosphere`
writes it; a subtype only supplies the dataset-specific parts.

Unlike the forcing grid, the atmosphere grid is independent of the ocean grid — NumericalEarth
interpolates between them — so it only has to cover the ocean domain with a margin.

# Fields a subtype provides
- `data_root`, `output_file`, `plot_file`: resolved by `atmosphere_path` and `plot_path`.
- `output_directory`: resolved by `atmosphere_directory`; the directory the download writes its
  intermediate files into. Only needed by a dataset that downloads.
- `resolution`: target grid spacing in degrees.
- `padding`: degrees of margin added around the ocean domain, read by
  `atmosphere_target_axes`.

# Methods a subtype provides
- `atmosphere_time_steps(config)`: the downloaded time records, as `AtmosphereRecord`s.
  Required.
- `atmosphere_source_grid(config, filepath)`: geometry of the downloaded data, e.g. a
  `ProjectedAtmosphereGrid`. Required.
- `atmosphere_variable_names(config)`: source variable name => prepared variable name.
  Required.
- `fetch_atmosphere(target_grid, config)`: fetch the source data. Only if it downloads.
- `prescribed_atmosphere(config, architecture; reference_date)`,
  `prescribed_radiation(config, architecture; reference_date)`: read the prepared file back at
  simulation time, as the NumericalEarth objects `coupled_model` consumes.
  Only if the setup is simulated.
- `atmosphere_date_range(config)`: first and last date of the prepared file, for
  `build_simulation`'s coverage check. Optional; the default is `nothing`, which skips the check.

`FjordSim.Atmospheres.NORA3Config` is the built-in implementation, for the MET Norway NORA3
reanalysis; `src/Atmospheres/nora3_source.jl` is the template to copy for a new dataset.
"""
abstract type AbstractAtmosphereConfig end

"""
    AbstractRunConfig

Supertype for run configurations. A concrete subtype describes how a built model is *run*:
where it writes, what it reports, how it steps, where its state starts and over which window. What
the model *is* lives in the setup's `AbstractModelConfig`, not here.

The pipeline reads a run config only through the methods below, never through its fields, so a
subtype stores or derives each value as it likes. The writers, callbacks, time stepping and initial
conditions it returns each dispatch on their own supertype. A different way of running — an
ensemble, a nested run — is a setup type of its own overloading `run_simulation`.

# Methods a subtype provides
- `results_root(config)`: the directory the run writes into, and what a writer's `output_file` and
  `FjordSim.CLI.log_path` resolve against. The one config rooted at a results directory rather
  than a `data_root`, being the only one that writes rather than reads.
- `FjordSim.Simulations.run_architecture(config)`: the live `CPU()` or `GPU()` the run uses.
- `writers(config)`: a tuple of `AbstractWriterConfig`s, attached by `attach_writer!`.
- `callbacks(config)`: a tuple of `AbstractCallbackConfig`s, attached by `attach_callback!`.
- `time_stepping(config)`: an `AbstractTimeSteppingConfig`, attached by `attach_time_stepping!`.
- `initial_conditions(config)`: an `AbstractInitialConditionsConfig`, or a `NamedTuple` `set!`
  accepts, resolved by `FjordSim.Simulations.resolve_initial_conditions`.
- `start_date(config)`: the calendar instant model time zero stands for. With `stop_time` it is the
  run window `coverage_window` reports, which every prepare step spans.
- `stop_time(config)`, `loops(config)`, `pickup(config)`: run control, read by `build_simulation`
  and `run_simulation`.

`FjordSim.Simulations.RunConfig` is the built-in implementation, returning its fields of the same
names.
"""
abstract type AbstractRunConfig end

function results_root end
function writers end
function callbacks end
function time_stepping end
function initial_conditions end
function start_date end
function stop_time end
function loops end
function pickup end

"""
    AbstractInitialConditionsConfig

Supertype for where a run's ocean state starts from, as a run config's `initial_conditions`. A plain
`NamedTuple` of constants, fields or functions is accepted beside it and passes through unchanged.

# Methods a subtype provides
- `FjordSim.Simulations.resolve_initial_conditions(initial_conditions, grid, config, tracers)`:
  the `NamedTuple` `set!(model; ...)` consumes. Required. `config` is the whole setup, so a source
  reads whichever component holds its data; `tracers` are the model config's `model_tracers`.

`FjordSim.Simulations.FromForcing`, `FromInitialState` and `FromResults` are the built-in
implementations.
"""
abstract type AbstractInitialConditionsConfig end

"""
    AbstractModelConfig

Supertype for model configurations. A concrete subtype names the model components that do not
depend on the grid, and its `coupled_model` method assembles them, on a grid, into the model a
`Simulation` steps.

`build_simulation` derives the grid, the forcing, the boundary conditions, the initial state and the
atmosphere and hands them to this one call as `model_inputs`, so a different model is a new subtype
rather than an edit to an existing method body.

# Methods a subtype provides
- `coupled_model(model, grid; forcing, boundary_conditions, initial_conditions, biogeochemistry,
  atmosphere, radiation, boundary_data_config, river_columns, stop_time, initial_time_step)`: build
  and return the model. Required. `boundary_data_config` is the setup's `AbstractBoundaryDataConfig`
  (or `nothing`), so a component acting on the open edges reads them from the one place they are
  stated. `river_columns` is every column a river enters through, or `nothing`; `model_closure`
  takes it. `stop_time` and `initial_time_step` are there because NumericalEarth's
  `OceanSeaIceModel` wraps the ocean as a `Simulation` of its own.
- `model_tracers(model)`: the tracer names the model carries, as a tuple of `Symbol`s. Required —
  `build_simulation` needs them before the model exists, to pick which forcing terms to build,
  which lateral tracer boundaries to open and which state variables to read back.
- `model_boundary_conditions(model, pieces, grid, forcing, boundary_data_config, tracers,
  boundaries)`: turn the setup's `boundary_condition_configs` tuple into the objects the model consumes.
  Optional; the default merges the pieces with `FjordSim.Utils.recursive_merge` and wraps each field
  in Oceananigans `FieldBoundaryConditions`.

`FjordSim.Simulations.CoupledHydrostaticModel` is the built-in implementation.
"""
abstract type AbstractModelConfig end

"""
    AbstractFreeSurfaceConfig

Supertype for free-surface configurations. A concrete subtype names a free-surface solver's knobs,
and its `free_surface(config, grid)` method builds the object `coupled_model` passes to
`HydrostaticFreeSurfaceModel`.

Nested inside a model config as its own field rather than a plain scalar, for the same reason the
model itself is a config rather than the built object: `SplitExplicitFreeSurface(grid, cfl = ...)`
needs the grid, which does not exist until `coupled_model` is called.

# Methods a subtype provides
- `free_surface(config, grid)`: build the free-surface object. Required.

`FjordSim.Simulations.SplitExplicitFreeSurfaceConfig` is the built-in implementation.
"""
abstract type AbstractFreeSurfaceConfig end

"""
    AbstractClosureConfig

Supertype for *closure* configurations. A concrete subtype names a turbulence closure's knobs, and
its `model_closure(config, grid, boundary_data_config, river_columns)` method builds the closure — or
closure tuple — that `coupled_model` passes to `HydrostaticFreeSurfaceModel`.

The `closure` field of a model config holds either a pre-built Oceananigans closure
(the usual case, and what every setup wrote before this supertype existed) or one of these. Which
one it is is resolved by dispatch on `model_closure`, the way `initial_conditions` is resolved by
`FjordSim.Simulations.resolve_initial_conditions` — three methods rather than a branch, and an
identity fallback so a pre-built closure passes through untouched and no existing setup changes.

It exists because a closure can need what only exists at build time. `FjordSim.Simulations.BoundarySponge`
needs two such things: the grid, for the domain extent and cell size its ramp is expressed in, and
the *open edges*, which are stated once on the `AbstractBoundaryDataConfig` and read back with
`open_edges` — so a setup that opens a second edge cannot end up sponging only the first.
`FjordSim.Simulations.RiverMouthMixing` needs the columns the rivers enter through, which only the
river file `add_rivers` wrote knows.

# Methods a subtype provides
- `model_closure(config, grid, boundary_data_config, river_columns)`: build the closure. Required.
"""
abstract type AbstractClosureConfig end

"""
    AbstractBiogeochemistryConfig

Supertype for *biogeochemistry* configurations. A concrete subtype names a biogeochemical model's
knobs, and its `model_biogeochemistry(config, grid)` method builds the object `coupled_model`
passes to `HydrostaticFreeSurfaceModel`.

A `FjordConfig` field of its own, `biogeochemistry_config`, beside the model config rather than
inside it: a biogeochemical model is swapped independently of the physics it is coupled to. The
model config still lists its tracers in `model_tracers`, because `build_simulation` needs them before
anything is built. It exists because biogeochemical models such as OceanBioME's `NPZD(grid; ...)`
need the grid, which does not exist when a setup file is evaluated.

# Methods a subtype provides
- `model_biogeochemistry(config, grid)`: build the biogeochemistry. Required.
"""
abstract type AbstractBiogeochemistryConfig end

"""
    AbstractBoundaryConditionConfig

Supertype for boundary-condition configurations. A concrete subtype describes one *contribution* to
the model's boundary conditions, and its `boundary_condition_sides` method builds that contribution.

A setup composes its boundaries out of several of these — air-sea fluxes, bottom drag and an open
lateral edge are three independent statements — as the `boundary_condition_configs` tuple of its
`FjordConfig`, and opts out of a piece by leaving it out. The model config's
`model_boundary_conditions` combines them, by default merging with `FjordSim.Utils.recursive_merge`
in tuple order (last wins). So each method here returns the *nested named tuple* shape
(`(u = (top = …, bottom = …), T = (top = …,))`), not `FieldBoundaryConditions`.

# Methods a subtype provides
- `boundary_condition_sides(config, grid, forcing, boundary_data_config, tracers[, boundaries])`: the
  nested named tuple this piece contributes, keyed by field and then by side. Required. Implement the
  six-argument form only if the piece reads the prepared open-boundary state (`boundaries`, from
  `FjordSim.BoundaryData.boundary_series`, or `nothing`); the five-argument form is reached through
  a forwarding fallback, which is what lets a piece that carries its own exterior state ignore it.

`FjordSim.BoundaryConditions.AirSeaFluxes`, `FjordSim.BoundaryConditions.QuadraticBottomDrag` and
`FjordSim.BoundaryConditions.OpenLateralBoundaryFromData` are the built-in implementations. Not to
be confused with `AbstractBoundaryDataConfig`, which is the exterior *data* rather than the scheme —
and which is what the fourth argument holds, so a piece can read `open_edges` from it.
"""
abstract type AbstractBoundaryConditionConfig end

"""
    AbstractWriterConfig

Supertype for output-writer configurations. A concrete subtype describes one thing a run writes, and
its `attach_writer!` method attaches it to the simulation.

A tuple of these replaces the hardcoded snapshot writer and the `checkpoint_interval` threshold that
used to decide whether a `Checkpointer` existed: a setup that wants no checkpoints names no
checkpointing writer, and a setup that wants two snapshot files at different intervals names two
writers. Which simulation a writer attaches to — the coupled one or the ocean sub-simulation — is
the method's business, not the caller's.

# Fields a subtype provides
- `output_file`: resolved by `results_path`, **only** for a writer whose `output_path_trait` is
  `FjordSim.Simulations.NamesOutputFile`. A writer that names no product file (a checkpointer,
  whose files carry no run tag and are pruned) needs no such field, and `results_path` reports that
  as an error rather than reaching for a field that is not there.

# Methods a subtype provides
- `attach_writer!(simulation, writer, config, loop)`: attach this writer for repetition `loop`.
  Required.
- `FjordSim.Simulations.checkpoint_trait(writer)`: whether this writer contributes a `Checkpointer`.
  Defaults to `NotCheckpointing()`.
- `FjordSim.Simulations.output_path_trait(writer)`: whether this writer names a file the run should
  report. Defaults to `NamesNoOutputFile()`.

`FjordSim.Simulations.SnapshotWriter` and `FjordSim.Simulations.CheckpointWriter` are the built-in
implementations.
"""
abstract type AbstractWriterConfig end

"""
    AbstractCallbackConfig

Supertype for callback configurations. A concrete subtype describes one thing a run does *while* it
runs — report progress, accumulate a diagnostic, watch for a stopping condition — and its
`attach_callback!` method installs it on the simulation.

The same shape as `AbstractWriterConfig`, and for the same reason: what a run reports used to be a
single hardcoded `Callback(progress, TimeInterval(config.progress_interval))` under the fixed key
`:progress`, so changing the report, the schedule type, or adding a second diagnostic all meant
rewriting `build_simulation`. A tuple of these makes each of those a config choice, and a setup that
wants a silent run names an empty tuple.

# Fields a subtype provides
- `name`: the key it occupies in `simulation.callbacks`. Two callbacks under one name replace each
  other, which `FjordSim.Simulations.validate_callbacks` rejects.

# Methods a subtype provides
- `attach_callback!(simulation, callback, config)`: attach this callback. Required. Which simulation
  it attaches to — the coupled one or the ocean sub-simulation — is the method's business, exactly
  as it is for a writer.

`FjordSim.Simulations.ProgressCallback` is the built-in implementation.
"""
abstract type AbstractCallbackConfig end

"""
    AbstractTimeSteppingConfig

Supertype for time-stepping configurations. A concrete subtype describes how the model's time step
is chosen, and its `attach_time_stepping!` method installs that choice on the simulation.

# Methods a subtype provides
- `attach_time_stepping!(simulation, config)`: install the time-step policy. Required. It must not
  set `simulation.Δt` — the step a loop starts from is `initial_time_step`, applied once when the
  simulation is built, and `FjordSim.Simulations.restart_loop!` deliberately leaves `Δt` alone so a
  repetition keeps the step the previous one converged on.
- `initial_time_step(config)`: the `Δt` the simulation starts at, in seconds. Required.

`FjordSim.Simulations.AdaptiveTimeStep` is the built-in implementation.
"""
abstract type AbstractTimeSteppingConfig end

"""
    AbstractFjordConfig

Supertype for a whole setup: what a config file evaluates to, and what every pipeline step takes
from the command line.

Every pipeline step (`prepare_bathymetry`, …, `run_simulation`, `validate_simulation`) and
`build_simulation` has exactly one method, on this type, and holds no logic of its own: it reads the
components it needs through the accessor functions below and passes them to a function of another
name that does the work, dispatching on those components alone — `prepare_bathymetry(config)` is
`process_bathymetry(grid_config(config), bathymetry_config(config), boundary_data_config(config))`.
So a subtype reuses every step by implementing the accessors, and replaces one step, or adds a kind
of input, by overloading that one function on itself rather than copying a body.

# Methods a subtype provides
- The accessors `grid_config`, `bathymetry_config`, `forcing_config`, `initial_state_config`,
  `river_config`, `boundary_data_config`, `atmosphere_config`, `model_config`,
  `biogeochemistry_config`, `boundary_condition_configs` and `run_config`, each returning that
  component or `nothing`. Required.
- `FjordSim.Simulations.model_inputs(config, grid, run)`: the keywords
  `coupled_model` receives. Optional; overload it to feed the model an input of a new kind.
- `FjordSim.CLI.pipeline_steps(config)`: the steps the command line offers, as `name => function`
  pairs in run order. Optional; overload it to add a step for an input of a new kind.

`FjordConfig` is the built-in implementation.
"""
abstract type AbstractFjordConfig end

"""
    FjordConfig

Complete setup configuration for one fjord, as returned by a config file loaded with
`fjord_config`.

The fields are parametric over the abstract config supertypes rather than fixed to the
built-in types, so a new grid, bathymetry source or forcing dataset only needs a struct
subtyping the matching supertype plus methods overloaded on it — `FjordConfig` and the
functions taking it are untouched. Each instantiation still has concrete field types.

Every component is a field of its own, side by side: the data pipelines, the model, its
biogeochemistry, its boundary conditions and the run. Each step dispatches on the few it reads.

# Fields
- `grid_config`: any `AbstractGridConfig`.
- `bathymetry_config`: any `AbstractBathymetryConfig`.
- `forcing_config`: any `AbstractForcingConfig`, or `nothing` for a setup that runs with no
  forcing at all. Defaults to `nothing`, so a setup opts in by naming one.
- `initial_state_config`: any `AbstractInitialStateConfig`, or `nothing` for a setup whose initial
  conditions come from elsewhere. Defaults to `nothing`. Independent of `forcing_config`: the
  initial state is one record at `start_date`, in a file of its own.
- `river_config`: any `AbstractRiverConfig`, or `nothing` for a setup with no rivers. Defaults to
  `nothing`, so a setup opts in by naming one. Independent of `forcing_config`: the river file is
  its own, written by `add_rivers` from the bathymetry alone.
- `boundary_data_config`: any `AbstractBoundaryDataConfig`, or `nothing` for a setup whose lateral
  boundaries are closed. Defaults to `nothing`, so a setup opts in by naming one. Independent of
  `forcing_config`: the exterior state along the open edge is its own hourly file from its own
  pipeline, and either config can be named without the other.
- `atmosphere_config`: any `AbstractAtmosphereConfig`, or `nothing` for a setup that prepares no
  atmosphere. Defaults to `nothing`, so a setup opts in by naming one.
- `model_config`: any `AbstractModelConfig`, or `nothing` for a setup that is only prepared.
  Defaults to `nothing`; a setup naming a `run_config` must name one.
- `biogeochemistry_config`: any `AbstractBiogeochemistryConfig`, or `nothing` for a model with no
  biogeochemistry. Defaults to `nothing`.
- `boundary_condition_configs`: a tuple of `AbstractBoundaryConditionConfig`s, merged in order (last wins),
  or `nothing` for a setup that is only prepared. Defaults to `nothing`; a setup naming a
  `run_config` must name one, `()` for the model's defaults everywhere.
- `run_config`: any `AbstractRunConfig`, or `nothing` for a setup that is only
  prepared and not run. Defaults to `nothing`, so a setup opts in by naming one.

# Example

```julia
using FjordSim: Grids
using FjordSim.Configs: AbstractGridConfig

struct SingleColumn <: AbstractGridConfig
    depth::Float64
end

Grids.domain_grid(config::SingleColumn, architecture) = ...      # new behavior
Grids.simulation_grid(config::SingleColumn, file, architecture) = ...
```
"""
Base.@kwdef mutable struct FjordConfig{
    G<:AbstractGridConfig,
    B<:AbstractBathymetryConfig,
    F,
    I,
    R,
    D,
    A,
    M,
    BG,
    BC,
    S,
} <: AbstractFjordConfig
    grid_config::G
    bathymetry_config::B
    forcing_config::F = nothing
    initial_state_config::I = nothing
    river_config::R = nothing
    boundary_data_config::D = nothing
    atmosphere_config::A = nothing
    model_config::M = nothing
    biogeochemistry_config::BG = nothing
    boundary_condition_configs::BC = nothing
    run_config::S = nothing
end

grid_config(config::FjordConfig) = config.grid_config
bathymetry_config(config::FjordConfig) = config.bathymetry_config
forcing_config(config::FjordConfig) = config.forcing_config
initial_state_config(config::FjordConfig) = config.initial_state_config
river_config(config::FjordConfig) = config.river_config
boundary_data_config(config::FjordConfig) = config.boundary_data_config
atmosphere_config(config::FjordConfig) = config.atmosphere_config
model_config(config::FjordConfig) = config.model_config
biogeochemistry_config(config::FjordConfig) = config.biogeochemistry_config
boundary_condition_configs(config::FjordConfig) = config.boundary_condition_configs
run_config(config::FjordConfig) = config.run_config

end  # module Configs
