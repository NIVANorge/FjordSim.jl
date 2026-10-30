# Adding a new source

> Read this before subtyping any `Abstract*Config` supertype — it lists
> every hook a new source must or may implement, with defaults and rationale.

Every pipeline is a generic function on the config supertype plus a small set of hooks. Add a
source by subtyping and overloading the hooks — never by editing the generic function. The
adapter files (`src/Bathymetry/geonorge.jl`, `src/Forcing/norkyst.jl`,
`src/InitialState/norkyst_initial_state.jl`, `src/Rivers/nve_rivers.jl`,
`src/Boundaries/BoundaryData/norkyst_boundary_data.jl`) are the templates.

Neither the supertypes nor the hooks are exported — `using FjordSim` gives a setup only what it
names. Code adding a source imports the supertype from `FjordSim.Configs` and extends each hook
through the module that declares it — the pipeline's own (`Rivers`, `BoundaryData`, ...), or
`Grids` for `open_edges` and `Plotting` for `plot_path`:

```julia
using FjordSim: Rivers
using FjordSim.Configs: AbstractRiverConfig

struct MyRivers <: AbstractRiverConfig
    # ...
end

Rivers.river_locations(config::MyRivers) = ...
```

Where an adapter's code goes:

- **The module file holds only what every source of the pipeline shares.** `Forcing.jl`,
  `BoundaryData.jl` and the rest import `Configs`, `Grids`, `Regridding`, `Utils` and `Plotting`, never a
  dataset. Everything one dataset needs — its config type, its hook methods, its download, and the
  `using` lines for them — goes in the adapter file the module `include`s at its end.
- **A dataset more than one pipeline reads gets its own module in `src/DataSources/`**, holding what
  is the dataset's rather than a pipeline's: endpoints, variable names, source grid, subset download.
  Each adapter imports from it (`using ..NorKyst: …`). A dataset only one pipeline reads — NORA3,
  NVE, Geonorge — keeps all of that in its adapter.
- **A helper that names no dataset goes in `Regridding`** — a file layout, a mask, a rotation — not
  in the first pipeline that needed it, where a second pipeline would have to import a sibling.

Grid — `AbstractGridConfig`:

| Hook | Required |
|---|---|
| `domain_grid(config, architecture)` → the bare grid, no bathymetry | yes |
| `simulation_grid(config, bathymetry_file, architecture)` → the grid the simulation runs on | yes |

Both are declared in `Grids`, which is included before every pipeline, since `Bathymetry`,
`Atmospheres`, `Forcing`, `Rivers` and `BoundaryData` all call `domain_grid`. Neither names a result type, so a grid that is not a
`LatitudeLongitudeGrid` can implement them.

Bathymetry — `AbstractBathymetryConfig`, consumed by `prepare_bathymetry`:

| Hook | Required | Default |
|---|---|---|
| `bathymetry_dataset(target_grid, config)` → NumericalEarth dataset | yes | none |
| `regrid_options(config)` → NamedTuple for `regrid_bathymetry` | no | `(;)` |
| `smoothing_options(config)` → NamedTuple for `smooth_bathymetry_gaps!` — `open_boundary_land_cells`, `max_island_cells`, `close_narrow_passages`, `spike_ratio`, `minimum_cell_fraction`, `max_slope_factor`, `minimum_depth`, `preserve_sills`, `raise_boundary_pockets`. `preserve_sills = true` also has `prepare_bathymetry` read the source back on its own grid (`native_bathymetry`), which needs the dataset's `longitude_interfaces`, `latitude_interfaces`, `size` and a 2D variable in its file | no | `(;)`, leaving only the topological cleanup |

The `edges` keyword `process_bathymetry` takes is a pipeline *argument*, not a hook, exactly as
`regrid_forcing`'s is: the step reads it with `open_edges(boundary_data_config)` and every source
inherits `clear_open_boundary_land` unchanged.

Forcing — `AbstractForcingConfig`, consumed by `prepare_forcing`:

| Hook | Required |
|---|---|
| `forcing_time_steps(config)` → `Vector{SourceRecord}` | yes |
| `forcing_source_grid(config, filepath)` → source grid | yes |
| `forcing_variable_names(config)` → `Dict` source name => FjordSim name | yes |
| `fetch_forcing(target_grid, config)` | only if it downloads |
| `simulation_forcing(config, grid, filepath, tracers, reference_date)` → the forcing term object `coupled_model` consumes | no; defaults to `forcing_from_file`, the FjordSim NetCDF contract |
| `forcing_date_range(config, filepath)` → `(first, last)` `DateTime`s | no; defaults to reading `ds["time"]` from the prepared NetCDF. A source that overrides `simulation_forcing` overrides this too |
| `source_field_grid(source, architecture)`, `projected_target_nodes(longitude, latitude, source)` | only for a source grid that is not a regular projected grid; dispatch on the source-grid type, not the config |

Initial state — `AbstractInitialStateConfig`, consumed by `prepare_initial_state`:

| Hook | Required |
|---|---|
| `initial_state_record(config, date)` → the `SourceRecord` holding the state at `date` | yes |
| `initial_state_source_grid(config, filepath)` → source grid | yes |
| `initial_state_variable_names(config)` → `Dict` source name => FjordSim name | yes |
| `fetch_initial_state(target_grid, config, date)` | only if it downloads |

`date` is the run config's `start_date`, passed by the setup-level step: the dataset
decides which of its records stands for that instant, and fetches only that one.

Rivers — `AbstractRiverConfig`, consumed by `add_rivers`:

| Hook | Required | Default |
|---|---|---|
| `river_locations(config)` → `Vector{RiverLocation}` | yes | none |
| `river_discharge(config, times)` → `(river, time)` matrix in m³/s, finite and non-negative | yes | none |
| `river_series(config, times)` → `Dict` FjordSim tracer name => `(river, time)` matrix of what the water carries; finite for every river the grid carries (no "ambient" value — unstable), and a tracer it omits enters at zero | yes | none |
| `fetch_rivers(target_grid, config)` | only if it downloads | none |
| `river_search_radius(config)` → cells to search for a coastal cell | no | `config.search_radius` |
| `river_maximum_inflow_velocity(config)` → m/s, the fastest a river's water may cross the surface of its columns | no | `config.maximum_inflow_velocity` |

A river is a volume of water: `add_rivers` snaps each mouth to a coastal cell, gives it as many
connected surface columns as its peak discharge needs to enter no faster than
`maximum_inflow_velocity`, and writes the result — inflow speed and river values per column, on a
daily axis over the run window — to a river file of its own. The run reads it as a free-surface source
and a top-cell tracer source; see "Rivers as a volume source" in `docs/architecture.md` (Rivers
section) for why that pair, and why the columns are spread. The speed bound is a vertical Courant
number on the top cell, so it has no default: it belongs to the setup's grid and time step.

A river config is a `FjordConfig` field of its own, `river_config`, `nothing` for a setup with no
rivers, and independent of the forcing config. The file it writes stands beside the prepared forcing
rather than inside it, so `add_rivers` needs only the bathymetry and the run window, and a setup can
name rivers with no forcing at all. Every `river_series` variable
is written; the run uses those naming a model tracer and ignores the rest, so one river dataset can
serve models with different tracers.

Open-boundary data — `AbstractBoundaryDataConfig`, consumed by `prepare_boundary_data`:

| Hook | Required | Default |
|---|---|---|
| `boundary_time_steps(config)` → `Vector{SourceRecord}` | yes | none |
| `boundary_source_grid(config, filepath)` → source grid | yes | none |
| `boundary_variable_names(config)` → `Dict` source name => FjordSim name | yes | none |
| `fetch_boundary_data(target_grid, config)` | only if it downloads; reads its edges with `open_edges(config)`, and subsets the box `boundary_domain` derives from them | none |
| `boundary_date_range(config)` → `(first, last)` `DateTime`s | no | reads `ds["time"]` from the prepared NetCDF |
| `boundary_source_slab(config, reader, step, source_name)` → one source slab | no | `blended_slab`, a plain read. Override only for a variable *derived* from more than one of the source's own — a velocity component stated along the source's own grid axes is the case that needs it |

Like a river config, this is a `FjordConfig` field of its own, `boundary_data_config`, `nothing` for a
setup whose lateral boundary is not data-driven: the exterior state is its own hourly file from its own
collection read by its own pipeline, and either config is nameable without the forcing one.

The edges **are** a field of it, `open_edges`, read through the `open_edges` accessor and stated nowhere
else. Write one `Symbol` or a collection of them — a domain in the open ocean names all four — and the
built-in config's keyword constructor normalizes either to a `Vector{Symbol}` with `lateral_edges`, so
every consumer iterates one shape rather than branching on a scalar. The three pipeline entry points
take the config and read the edges off it rather than taking both, which is what keeps them from
disagreeing; `regrid_forcing` takes them as an `edges` keyword, since a forcing dataset has no
business naming them, and `nothing` there means a closed domain. Not to be
confused with `AbstractBoundaryConditionConfig`, which is the scheme rather than the data.

`prepare_boundary_data` reuses the `Regridding` core on a one-cell-thick target slab, so a boundary source on
a regular projected grid inherits `source_field_grid`, `projected_target_nodes`, `SourceFill` and the
interpolation kernel unchanged, exactly as a forcing source does. It loops the whole per-variable build
over the config's edges and writes them all into one file, so neither hook above becomes per-edge — the
cost of that is the download: `boundary_domain` returns the bounding box of the edges' bands, which for
two opposite edges is the whole domain plus the margin, so a multi-edge download is the full box at
hourly cadence and its interior is fetched and never read. Per-edge downloads are the optimization, at
the price of `boundary_time_steps` and `boundary_source_grid` becoming per-edge. What the prepared file's variables
are *named* and *shaped* is not a hook — `boundary_variable_name`, `boundary_dimension_names` and
`boundary_location` fix it, because the read side does.

Atmosphere — `AbstractAtmosphereConfig`, consumed by `prepare_atmosphere`:

| Hook | Required |
|---|---|
| `atmosphere_time_steps(config)` → `Vector{AtmosphereRecord}` | yes |
| `atmosphere_source_grid(config, filepath)` → source grid, e.g. `ProjectedAtmosphereGrid` | yes |
| `atmosphere_variable_names(config)` → `Dict` downloaded name => prepared name | yes |
| `fetch_atmosphere(target_grid, config)` | only if it downloads |
| `prescribed_atmosphere(config, architecture; reference_date)` → `PrescribedAtmosphere` | only if the setup is simulated |
| `prescribed_radiation(config, architecture; reference_date)` → `PrescribedRadiation` | only if the setup is simulated |
| `atmosphere_date_range(config)` → `(first, last)` `DateTime`s | no; defaults to `nothing`, skipping the coverage check |

The last three are the read side, consumed by `build_simulation` rather than `prepare_atmosphere`,
and are what keep `Simulations` from naming a dataset. The first two take no float type on
purpose: both NumericalEarth constructors default to `Float32`, and passing
`Oceananigans.defaults.FloatType` would silently promote the atmosphere to `Float64`. A `nothing`
atmosphere config yields `nothing` for all three.

`reference_date` is the instant the returned time axes are zeroed at, and is *not*
`NORA3PrescribedAtmosphere`'s `start_date`: that one selects which records to load, this one
selects where t = 0 sits. `build_simulation` passes the run config's `start_date` to both
this hook and `simulation_forcing`, which is the only thing keeping the two in phase — see
`docs/architecture.md` (Simulations section).

The prepared variable names and units are *not* a hook — they are fixed by the read side in
`ATMOSPHERE_VARIABLES`. A source whose download already normalizes names (as NORA3's does, since
five of its eight variables are derived rather than copied) returns the identity mapping.
`Regridding`'s `grid_rotation_angle` and `rotate_to_east_north` are available to any adapter whose
source gives wind relative to its own grid axes.

Model — `AbstractModelConfig`, the setup's `model_config`, consumed by `build_simulation`. Its
biogeochemistry and its boundary conditions are setup components of their own, beside it:

| Supertype | `FjordConfig` field | Hooks |
|---|---|---|
| `AbstractModelConfig` | `model_config` | `coupled_model(model, grid; forcing, boundary_conditions, initial_conditions, biogeochemistry, atmosphere, radiation, boundary_data_config, river_columns, stop_time, initial_time_step)` → the model, and `model_tracers(model)`: required. `model_boundary_conditions(model, pieces, grid, forcing, boundary_data_config, tracers, boundaries)`: optional, defaults to merging the pieces into Oceananigans `FieldBoundaryConditions` |
| `AbstractBiogeochemistryConfig` | `biogeochemistry_config` | `model_biogeochemistry(config, grid)` → the biogeochemistry `coupled_model` receives; required. `examples/oslofjorden_npzd.jl` has one |
| `AbstractBoundaryConditionConfig` | `boundary_condition_configs` (a tuple, merged in order) | `boundary_condition_sides(config, grid, forcing, boundary_data_config, tracers[, boundaries])` — implement the six-argument form only if the piece reads the prepared open-boundary state; the five-argument one is reached through a forwarding fallback |

The boundary-condition pieces are a plain tuple because how they combine, and what shape the result
takes, is the model's business: a model whose boundary objects are not `FieldBoundaryConditions`
overloads `model_boundary_conditions` on its own config.

Run — `AbstractRunConfig`, the setup's `run_config`: how a built model is run. The pipeline reads it
only through accessors (`FjordSim.Configs`), never its fields, so a subtype stores or derives each
value as it likes: `results_root`, `writers`, `callbacks`, `time_stepping`, `initial_conditions`,
`start_date`, `stop_time`, `loops`, `pickup`, and `Simulations.run_architecture`. A subtype
implementing only `results_root`, `start_date` and `stop_time` still gets `results_path`, `run_tag`
and `coverage_window`, so it can be prepared for; `build_simulation` needs the rest. A different way
of *running* — an ensemble, a nested run — is not a run config but a setup type of its own
overloading `run_simulation` (see below). What a run config returns is dispatched on one hook table
each:

| Supertype | Accessor | Required hooks |
|---|---|---|
| `AbstractCallbackConfig` | `callbacks` (a tuple) | `attach_callback!(simulation, callback, config)` |
| `AbstractWriterConfig` | `writers` (a tuple) | `attach_writer!(simulation, writer, config, loop)`; optionally `checkpoint_trait`, `output_path_trait` and `station_output`, all defaulting to the negative. `SnapshotWriter` (NetCDF), `FieldSnapshotWriter` (JLD2, for z-reduced fields NetCDF rejects), `SurfaceSnapshotWriter` (JLD2, the top level only, split into file parts) and `CheckpointWriter` are the built-ins; `StationWriter` and `FieldStationWriter`, whose `station_output` is `true`, are what `validate_simulation` scores |
| `AbstractTimeSteppingConfig` | `time_stepping` | `attach_time_stepping!(simulation, config)`, `initial_time_step(config)` |
| `AbstractInitialConditionsConfig` | `initial_conditions` (or a plain `NamedTuple`) | `resolve_initial_conditions(initial_conditions, grid, config, tracers)` → the `NamedTuple` `set!` takes. `config` is the whole setup, so a source reads whichever component holds its data. `FromForcing`, `FromInitialState` and `FromResults` are the built-ins |

A writer needs an `output_file` field only when its `output_path_trait` is `NamesOutputFile`; on any
other, `results_path` raises a stated `ArgumentError` rather than failing on a missing field.

`CoupledHydrostaticModel`, the built-in `AbstractModelConfig`, nests two further supertypes of its
own, so a different free surface or closure is a new subtype and a new method, not an edit to
`coupled_model`:

| Supertype | Field | Required hook | Built-in |
|---|---|---|---|
| `AbstractFreeSurfaceConfig` | `free_surface` | `free_surface(config, grid)` → the free-surface object `coupled_model` passes to `HydrostaticFreeSurfaceModel` | `SplitExplicitFreeSurfaceConfig` |
| `AbstractClosureConfig` | `closure` | `model_closure(config, grid, boundary_data_config, river_columns)` → the closure, or tuple of closures, it passes | `BoundarySponge`, `RiverMouthMixing` |

They differ in one way that matters. `free_surface` is *always* a config — there is nothing else
it could hold. `closure` holds a config **or** a pre-built Oceananigans closure, and usually the
latter; `model_closure`'s fallback is the identity, so the two shapes are resolved by dispatch and a
setup naming a plain closure needs no change and never sees the hook. `model_biogeochemistry` has the
same identity fallback, so a `nothing` biogeochemistry passes through. Same pattern as
`initial_conditions` and `resolve_initial_conditions`.

`model_closure` takes `boundary_data_config` as well as the grid because a closure may act on the open
edges, which are stated only on the `AbstractBoundaryDataConfig`, and `river_columns` — what
`read_river_columns` reads from the river file, or `nothing` — because one may act on the river
mouths. That is why `coupled_model` carries both too.

`CoupledHydrostaticModel` also carries `extra_kwargs`, which is not a hook but the reason a
setup rarely needs a *new* `AbstractModelConfig` subtype: it splats one `NamedTuple` per
slot into each of the three constructors `coupled_model` calls, so a model differing from the
built-in one by a keyword is a config value rather than a new type and a new method. See
`docs/architecture.md` (Simulations section) for the slot names and for why the constructor rejects
keys that would shadow the fields.

## A new kind of input, or a setup type of its own

Every pipeline step and `build_simulation` has exactly one method, on `AbstractFjordConfig`, and
holds no logic: it picks out the components it reads through accessor functions —
`grid_config(config)`, `forcing_config(config)`, … — and hands them to a work function of another
name, which dispatches on those components alone: `prepare_bathymetry(config)` is
`process_bathymetry(grid_config, bathymetry_config, boundary_data_config)`. Each optional component
has a `::Nothing` method of the work function, which is how a step becomes a no-op. The table of
steps and work functions is in `docs/architecture.md` (Configs section). So a setup type of its own
reuses every step by implementing the eleven accessors (`FjordSim.Configs`), and the test suite's
`DelegatingFjord` is the minimal example; and it replaces one step by overloading that step on
itself, `FjordSim.run_simulation(config::EnsembleFjord) = …`, which nothing else answers to.

An input of a kind `FjordConfig` has no field for — nutrient loads, wave forcing — is such a type
with one more field, and three overloads, none of them an edit to FjordSim:

```julia
using FjordSim: Configs, Simulations, CLI

struct NutrientFjord{C<:FjordConfig,N} <: Configs.AbstractFjordConfig
    base::C
    nutrient_config::N
end

Configs.grid_config(config::NutrientFjord) = Configs.grid_config(config.base)   # … every accessor

# The keyword the model receives: the default inputs plus the new one.
Simulations.model_inputs(config::NutrientFjord, grid, run) = merge(
    Simulations.model_inputs(config.base, grid, run),
    (; nutrients = nutrient_forcing(config.nutrient_config, grid)),
)

# The step that prepares it, offered on the command line after the built-in ones.
CLI.pipeline_steps(config::NutrientFjord) =
    [CLI.pipeline_steps(config.base); "prepare_nutrients" => prepare_nutrients]
```

The model config it pairs with takes the new keyword in its own `coupled_model` method. A new step
follows the same shape as the built-in ones: a one-line method on the setup type that picks out the
components, and a work function of another name on those components.

The `coverage` keyword `regrid_forcing`, `place_rivers`, `regrid_boundary_data` and
`regrid_atmosphere` take is a pipeline *argument*, not a hook: the steps derive it from the run config with `coverage_window`
and every source inherits the padding unchanged. No source gains a hook for it.

Required fields are listed in each supertype's docstring in `src/Configs.jl`, and its hooks under
"# Methods a subtype provides", one bullet per hook, saying `Required.`, `Optional` or `Only if …`.
`docs/make_dispatch_map.jl` reads that section to write `docs/dispatch-map.md`, so a new or renamed
hook goes there, and the map is regenerated with `julia --project docs/make_dispatch_map.jl`; the
"Dispatch map" test fails until it is. Path resolution
(`bathymetry_path`, `forcing_path`, `forcing_directory`, `initial_state_path`,
`initial_state_directory`, `river_forcing_path`, `boundary_data_path`,
`boundary_data_directory`, `atmosphere_path`,
`atmosphere_directory`, `results_path`, `plot_path`) and the plots come for free. A missing
required hook surfaces as a `MethodError` naming it, which the "Config extensibility" testset
asserts.

## Observations — `AbstractObservationConfig`

The one supertype in `src/Validation/`, and the only one whose subtypes run *after* a simulation
rather than before it. Same shape as every source above: subtype it, overload the hooks, and nothing
in the scoring or plotting code changes.

| Hook | Required | Default |
|---|---|---|
| `observation_stations(config)` | yes | — |
| `observation_series(config, station, variable; window)` | yes | — |
| `download_observations(config, window)` | no | no-op |
| `observations_directory(config)` | no | `joinpath(data_root, output_directory)` |

Required fields: `data_root` and `output_directory`.

Two conventions the hooks rely on, both because a validation sweep asks every source for every
station and variable it might hold:

- **`observation_series` returns `nothing`, not an error, for a pair it has nothing for.** Most
  combinations are empty — an instrument that was in the water for ten weeks of a two-year run, a
  station outside a programme's remit — and that is a blank row rather than a failure.
- **`download_observations` defaults to a no-op.** A source read from files someone supplies has
  nothing to fetch and should not have to say so. `KartverketSeaLevel` is the only built-in source
  that fetches, because public tide-gauge data is the exception: most programmes are held by their
  owners and want a reader, not a client.

Every source returns the same `ObservationSeries` — `(station, variable, times, depths, values,
units, source)`, with `values` shaped `(time, depth)` — which is what lets `metrics.jl` stay ignorant
of where a number came from. Profiles are interpolated onto one depth axis at read time, since no
instrument samples the same depths twice and a Hovmøller diagram needs a rectangular array;
`interpolate_to_depths` does not extrapolate, so a shallow cast cannot invent deep water.

`variable` is the FjordSim field name the observation is scored against — `"T"`, `"S"`, `"u"`,
`"v"`, `"eta"` — not the source's own name for it, so unit and naming conversions belong inside the
adapter where a mistake is visible.

## Station output — `StationWriter` and `FieldStationWriter`

Not a new supertype: both are `AbstractWriterConfig`s, listed in the writer table above. They are
mentioned here because a validation needs them and they are the only writers that do not write the
whole domain.

Each takes a `stations::Vector{Station}` of lon/lat positions and writes `indices = (i, j, :)` — one
whole water column per station, one file each. Horizontally staggered fields such as `u` and `v`
are interpolated to the cell centre first (`station_outputs`), since the column's own west or south
face may be a closed coastal one that holds zero. Snapping to the grid happens at attach time, where
the grid exists, and `report_station_cells` logs the grid indices, the model's water depth and the
offset in cells and metres for every station. That offset is the number a reader of the validation needs and
nothing else records: at a couple of hundred metres per cell a station in a narrow sound can be
placed a kilometre from its instrument, and a disagreement that large is a property of the
comparison rather than of the model.

`writer_keys` is overloaded to report one `output_writers` key per station, so `validate_writers`
still catches two writers that would replace each other.
