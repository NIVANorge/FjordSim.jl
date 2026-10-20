# Architecture

> The *why* behind each module: design rationale, ordering constraints, measured numbers, and the
> traps that fail silently. Not a map of the code — signatures, call order and file layout are in the
> code and its docstrings. Read the section for the module you are changing. `CLAUDE.md` has the
> operational essentials and a one-line summary of each module.

## What a run is assembled from

`build_simulation(config::FjordConfig)` reads only prepared files, each written by one pipeline step
onto (or along) the model grid. Only the bathymetry is required; every other input belongs to a
`FjordConfig` field defaulting to `nothing`, and a setup naming none skips its steps and its read:

| input | written by | used for |
|---|---|---|
| `bathymetry.nc` | `prepare_bathymetry` | the `ImmersedBoundaryGrid` itself, via `simulation_grid` |
| atmosphere file | `prepare_atmosphere` | `PrescribedAtmosphere`/`PrescribedRadiation` (NumericalEarth) |
| boundary file | `prepare_boundaries` | the exterior state on every open edge, hourly, tide included |
| river file | `add_rivers` | volume sources `η` and top-cell tracer sources |
| forcing file | `prepare_forcing` | the `FromForcing` initial state — and little else (below) |
| initial-state file | `prepare_initial_state` | the `FromInitialState` initial state: one record at `start_date` |
| a snapshot file | a previous run | the `FromResults` initial state |

**When a setup names a forcing config, its file is read and its terms are built** — a prerequisite
with its coverage validated, even for a run starting `FromResults`. `prepare_forcing` writes `λ = 0`
everywhere (see Forcing), so today those `ForcingFromFile` terms are zero and what the file supplies
is the `FromForcing` state. A setup that needs neither names no forcing config: nothing else depends
on it — the initial state, rivers, open boundaries and atmosphere each have their own config and file.

Every input is dispatched through a hook on an abstract supertype, so a new source is a subtype plus
methods (`docs/adding-a-source.md`).

## 1. Configs (`src/Configs.jl`)

Everything that more than one later module needs lives here, because the modules are included in a
line and none can see a later one:

- **The grid hooks** `domain_grid`/`simulation_grid` are declared here, not in `Grids`, because
  `Bathymetry`, `Atmospheres`, `Forcing`, `Rivers` and `Boundaries` call `domain_grid` and are all
  included before `Grids`. `model_closure` is declared here too: its fallback is the identity, which
  is only useful if it is the first method anyone sees.
- **The edge vocabulary** (`LATERAL_EDGES`, `validate_open_edge`, `lateral_edges`, `open_edges`):
  `Bathymetry` clears land along open edges, `Regridding` restores their velocity face rows,
  `Boundaries` regrids along them, `OpenLateralBoundaryFromData` and `BoundarySponge` act on them.
  `Configs` is the one module all of them can see.
- **The path helpers** (`bathymetry_path`, `forcing_path`, …) are defined on the supertypes so a new
  source inherits them without loading the built-in source's module.

Rivers and open-boundary data are `FjordConfig` fields of their own (`river_config`,
`boundary_config`), not parts of the forcing config: each is a separate file written by a separate
pipeline, and a setup can have either with no interior forcing. `AbstractBoundaryDataConfig` (the exterior *values*) and
`AbstractBoundaryConditionConfig` (the *scheme* acting on them) are configured separately for the
same reason.

**The run tag is the wall-clock launch instant** (`yyyymmddTHHMMSS`, `Configs.LAUNCH_TAG`, filled by
`__init__` — a `const` initialized with `now()` would freeze the precompilation instant). So it
distinguishes invocations, not configurations: re-running a config never overwrites an earlier run's
snapshots or log. Seconds are in it because a crash and an immediate relaunch otherwise collide.
Consequences: which window a file covers is not in its name but in the snapshot's `start_date`
attribute, and checkpoints carry no tag at all (see Checkpointing).

**No setup is part of the package.** A setup is a script ending in its `FjordConfig`, and
`fjord_config(path)` `Base.include`s it into `Main` on every call. So each step gets a fresh config —
the structs are mutable and `native_region!` edits the bathymetry config, so a shared instance would
leak state between steps — built after `__init__` has filled `DybdedataConfig`'s scratch
`raw_directory`, which is still `""` during precompilation. How to write one is `docs/setups.md`;
the measurements behind a setup's values are in its own comments.

`coverage_window(simulation_config)` is the `(first, last)` interval the run needs. Every prepare
step that has a time axis — forcing, rivers, boundaries, atmosphere — and `validate_simulation` take
it as their `coverage`.

## 2. Utils (`src/Utils.jl`)

`progress` takes every reduction on the `Field`, never on `interior(field)`. A reduction over a
`Field` on an `ImmersedBoundaryGrid` excludes the immersed periphery; `interior` returns a bare array
that includes it. Oceananigans writes `zero(eltype)` into every immersed tracer cell each
`update_state!`, and most of a fjord grid is land — so the array form reported the land mask as the
minimum temperature for the whole run. It also keeps the reduction on the GPU.

## 3. Plotting (`src/Plotting.jl`)

Included before the pipelines because each setup-level driver plots as its last step; `Plotting`
needs only `Configs`, so there is no cycle. Every `plot_*` dispatches on a config supertype so a new
source inherits it. Not tested, by rule.

## 4. Bathymetry (`src/Bathymetry/`)

`prepare_bathymetry` is regrid (NumericalEarth) → `smooth_bathymetry_gaps!` → write. The one
source-specific step is `bathymetry_dataset`; `geonorge.jl` implements it for Geonorge Sjøkart
Dybdedata.

### Gridding the soundings

Soundings are gridded by linear interpolation on their Delaunay triangulation
(`GEONORGE_GRIDDING_ALGORITHM`), so a basin keeps its deepest sounding. The `invdistnn` it replaced
had `smoothing=0.2` in *degrees* — kilometres, against soundings metres apart — so it was a plain
mean of 16 soundings that took 40 m off Bunnefjorden and 35 m off the Drøbak Sound before any
smoothing ran. Linear takes about a minute at `raw_resolution_factor = 2` against three hours. The
raw file's name carries `_linear`, so a cached inverse-distance file is never read as the new one.

### The smoothing stages, and why their order is forced

`smooth_bathymetry_gaps!` runs the topological cleanup every source gets (diagonal pairs, isolated
sea/land cells), then seven opt-in stages, each skipped when its `smoothing_options` parameter is
`false` or zero:

1. `clear_open_boundary_land`
2. `fill_small_islands`
3. `remove_narrow_passages`
4. `fill_shallow_spikes`
5. `snap_partial_bottom_cells` (first pass)
6. `limit_bottom_slope`
7. `preserve_sill_depths`
8. `snap_partial_bottom_cells` (second pass)

- Topological cleanup first: checkerboard noise would skew the neighbour medians.
- **Land-mask stages before depth stages**, or the cells about to change would already have
  contributed to a neighbour median and a slope pair. Applied after the depth stages, island
  flooding leaves `max r = 0.70` and two new spikes.
- `clear_open_boundary_land` and `fill_small_islands` only turn land into sea, so they cannot create
  a passage or strand a cell; running them first lets the passage stage judge the field as it will
  be. Closing a passage, by contrast, adds land exactly where islands sit and can bridge one to the
  mainland where the size test never sees it again.
- Despiking before slope limiting: a spike is the one-cell feature slope limiting would smear into
  its neighbours instead of removing.
- `preserve_sill_depths` after slope limiting, which would raise the channel it carves.
- **The snap runs twice.** Slope limiting ignores the vertical faces, so it puts back slivers the
  first snap removed: 2628 of them — 6 % of bottom cells — at `max_slope_factor = 0.25` on the
  24-level grid. The second snap leaves none, at the cost of raising the steepest slope from 0.250 to
  0.307 on 16 pairs of 84 000. The first snap is still needed: it gives the limiter less to move, and
  is the only snap when the limiter is off.

### Why the depth stages exist

`PartialCellBottom` bounds how *thin* a cell may be, not how much depth changes between adjacent
columns. A shallow cell beside a deep one carries the same transport in a fraction of the column, so
velocity grows until the time-step wizard goes NaN and `calculate_substeps` throws
`InexactError: Int64(NaN)` — the last symptom, not the cause. A `minimum_depth` floor alone makes
this worse by leaving a spike, so the two are configured together. `limit_bottom_slope` moves each
pair symmetrically, conserving volume, which in a narrow trench means the deep cell is raised once
per shallow neighbour: on `oslofjorden_validation`, `max_slope_factor = 0.25` took 34 m off
Vestfjorden at Fl1; 0.5 takes 8 m.

### `clear_open_boundary_land`

Floods land within `open_boundary_land_cells` rows of an open edge, so the boundary scheme sees the
clean channel it is derived for rather than a headland splitting the prescribed inflow.

- **The edges come from `open_edges(boundary_config)`**, passed as `prepare_bathymetry`'s `edges`
  keyword — not a bathymetry field, which could disagree.
- **Depths come by repeated relaxation**, each round computed in full before it is written, so a
  cleared headland ramps from the coast rather than stepping and the result does not depend on walk
  order.
- **The price is exterior data**: a newly wet cell has no profile in the boundary file and
  `fill_boundary_gaps!` fills it from the nearest wet cell. Keep the width to a few cells.

### `remove_narrow_passages`

Bounds channel **width**, which the depth stages cannot see. A one-cell-wide passage usually closes a
loop between already-connected basins, and the head difference around the loop is forced through one
shallow cell. On `oslofjorden` a 2.35 m canal at 10.43°E, 59.09°N carried a 47 m/s jet and collapsed
the time step from 180 s to 0.3 s.

- **A passage that is a basin's sole link is kept** (12 of 78 on `oslofjorden`).
- **Candidates are tested one at a time against the partly closed field**, so two passages that are
  each redundant while the other is open cannot together sever a basin. Deterministic, row-major.
- **One pass, not iterated**: iterating would erode a genuine narrow arm cell by cell.
- `wet_component` is hand-rolled rather than `ImageMorphology.label_components`, which reaches
  FjordSim only transitively.

### `fill_small_islands`

The dual: floods interior land patches of at most `max_island_cells`. On `oslofjorden` a three-cell
island near Horten carried the domain's velocity maximum, 2.67 m/s.

- `fill_isolated_land_cells` is its one-cell case and cannot be stretched to it, and the depth stages
  never see an island, which agrees with its surroundings in depth.
- **A component touching the domain edge is kept**: land continuing outside the domain may have only
  a few cells inside it.
- **No cleanup pass follows**: every neighbour of a flooded component is sea by maximality.

### `snap_partial_bottom_cells`

`PartialCellBottom` never makes a cell thinner than `ϵ Δz`; it pushes the bottom *down* until the
cell is that thick, so a sounding just below a face becomes a sliver whose floor is nowhere the
sounding was. Such a cell holds a fraction of its neighbours' water and is often cut off laterally,
so it is a **tracer** hazard: on `oslofjorden`'s open southern boundary one went 33 → 65 psu and
5 → 28 °C in 4.5 days while the boundary file asked for 33.2 psu. This stage snaps the bottom *up* to
the face instead. `minimum_cell_fraction` must equal the grid's `PartialCellBottom` fraction (0.2).

- **Refining `z_faces` does not fix this**: more faces means more laterally isolated bottom cells
  (3.9 % of columns on 18 levels, 5.5 % on 24).
- **Up, not down**: up rejoins the level the neighbours are on. Cost: 0.95 % of domain volume.
- **Guards**: a column is left alone when the face above is the surface or snapping would breach
  `minimum_depth`.
- **It snaps in the file's precision.** Smoothing runs in `Float64`; the file stores `Float32`
  (`BATHYMETRY_ELTYPE`). Where a face is not representable in `Float32`, narrowing puts the column a
  hair below it and recreates a 5e-8-thick sliver — 1265 columns on `oslofjorden`. `snap_to_face`
  rounds away from the water so the value survives.

### `preserve_sill_depths`

A cell mean is the wrong statistic for a channel narrower than the cell: Svelvik (180 m wide, 11 m
deep) came out 5.6 m on a 194 m grid. Taking each cell's deepest source water instead ruins slopes
(Drøbak 43 m against ~20 m). So the stage works on connections: the source's merge tree lists every
sill and the two basins it separates; where the model grid resolves both basins but does not connect
them at sill depth, a 0-1 BFS finds the path crossing the fewest too-shallow cells — using only cells
whose deepest source water reaches the sill — and deepens it. No land changes and nothing goes
deeper than the source. On `oslofjorden_validation`: 1467 cells, Svelvik 10.8 m, mean depth changed
by under 0.1 m. It only deepens, so a connection the regrid made too *deep* is out of reach.

## 5. Regridding (`src/Regridding.jl`)

The core `Forcing`, `InitialState` and `Boundaries` share; `Rivers` uses only the masks and the file
layout, and `Atmospheres` only the source-grid rotation.

- **The kernel runs where the config's `architecture` says** (`:auto`/`:cpu`/`:gpu`); `:auto` picks
  the GPU, ~12x faster than one CPU thread, and `:gpu` errors rather than falling back. A `Symbol`
  keeps config types concrete and configs loadable without a GPU. `target_grid` stays on the CPU
  regardless: building the masks walks `peripheral_node` cell by cell.
- **The land mask comes from `peripheral_node`** on the model grid, so it matches the model exactly,
  including `PartialCellBottom` and the velocity-face wall convention; open edges get their velocity
  face row restored.
- **Two library functions are deliberately not used.** NumericalEarth's dataset path always builds a
  `LatitudeLongitudeGrid`, but NorKyst's grid is rotated ~59° here, so the source is expressed as a
  `RectilinearGrid` in projected metres (`ProjectedSourceGrid`). And `inpaint_mask!` cannot fill a
  fully masked level — on this subset it never terminates or writes zeros — hence `SourceFill`.

### Time-axis padding

`pad_time_steps` extends an axis to both ends of the run window by replicating the nearest step.
Needed because NorKyst's daily records sit at 12:00, so a run starting at 00:00 has no record for its
first half-day, and every reader uses `Cyclical()` indexing, which wraps instead of failing.

- It runs **after** `daily_time_steps`: a pad at the record level would break that function's
  whole-day cadence test and silently disable gap-filling.
- A pad reaches **at most one record spacing** past the data (`forcing_record_spacing`). Unbounded,
  it would manufacture the coverage `validate_time_coverage` exists to verify. Overshooting is
  reported as a missing download.
- Padding changes `Cyclical`'s inferred period. That is harmless only because the loop resets the
  clock and never leaves `[t¹, tᴺ]`, so the cycling branch is dead — **do not "fix" the period.**

## 6. Atmosphere (`src/Atmospheres/`)

The prepared file is a **contract fixed by the read side** (`MultiYearNORA3`): `Float32`
`(lon, lat, time)`, uniform 1D axes (`compute_faces` infers spacing from the first difference),
CF `time`, the eight names and units in `ATMOSPHERE_VARIABLES`. Air temperature is **Kelvin** and
both radiative fluxes are **downwelling**: `PrescribedRadiation` applies its own albedo, so a net
flux would count it twice.

Interpolation is deliberately not shared with `Regridding`: that machinery is 3D, masked and
depth-oriented, atmospheric fields are 2D and defined over land, and the prepared grid is ~50x60.
The cost is ~15 duplicated lines of ArchGDAL transform.

`pad_atmosphere_records` runs after `validate_atmosphere_records`, whose hourly-gap warning is about
download integrity and would fire on the deliberately off-cadence pad.

### The NORA3 download (`nora3_source.jl`)

The download owns everything about the forecast structure and writes gap-free hourly months already
carrying the prepared names, so `prepare_atmosphere` is a pure regrid and re-running it at another
`resolution` costs no download.

- **Which lead supplies which hour**: leads 4–9 of each 6-hourly run tile the day with no overlap —
  04–09 from 00Z, …, and 00–03 from the *previous* day's 18Z.
- **De-accumulation**: flux accumulators restart every run, so `process_run` walks leads 3–9 and
  differences stay inside one run. Increments are labelled at the *end* of their interval, which puts
  fluxes on the instantaneous fields' axis (a midpoint label needs a second time axis
  `MultiYearNORA3` cannot read).
- **One DAP2 request per lead** (`read_lead`, `parse_dods`), paced by
  `NORA3_REQUEST_INTERVAL_SECONDS`: NetCDF's own OPeNDAP reader makes ~12 requests per lead and
  thredds.met.no returns 429 after ~6 rapid ones.
- **Winds are rotated on the native grid** before interpolation; `grid_rotation_angle` derives the
  angle from the 2D lon/lat fields (about −48° over Oslofjord).
- **The projection is the constant `NORA3_PROJ4`**: files from 2024 on describe it only by CF
  attributes, without the `proj4` string.
- met.no's `nora3_subset_atmos` collections are unusable: no `specific_humidity_2m`, and only *net*
  radiation, which cannot be inverted to downwelling.

## 7. NorKyst (`src/DataSources/NorKyst.jl`)

**One module for the dataset, one adapter per pipeline.** The forcing, the initial state and the
open boundaries all read NorKyst-800m, from four collections — the operational daily means and
hourly instants, and the Norkyst-v3 hindcast's `zdepth` and `sdepth` — that share one
polar-stereographic grid, projection variable and depth axis. What is the same for all of them
lives here: the endpoints, the variable names, `norkyst_source_grid`, and the download of a
longitude/latitude box (`subset_ranges`, `define_output_file`, `write_parameter_chunk!`). A config
type and its hooks belong to the adapter in the directory of the pipeline that consumes it.

`src/DataSources/` holds only datasets more than one pipeline reads. NORA3 (atmosphere), NVE
(rivers) and Geonorge (bathymetry) each feed one pipeline and stay in its directory.

## 8. Forcing (`src/Forcing/`)

**There is no interior relaxation.** `prepare_forcing` writes `λ = 0` everywhere: the open boundary
nudges towards hourly exterior data *at* the boundary, and a band relaxing the same variables a few
cells inside would fight it. So nothing damps the interior near the boundary except `BoundarySponge`
(Simulations); if more is needed, the candidate is a `u`/`v`-only band, which touches no tracer the
boundary nudges.

Forcing is optional: a setup naming no `forcing_config` skips both steps, builds no forcing terms,
and needs no forcing file — only `FromForcing` initial conditions do. One that names it gets its file
read: `ForcingFromFile` dispatches on `λ` (`> 1` x-flux, `< -1` y-flux, `|λ| < 1` relaxation) for
every variable, which with `λ = 0` costs a read of two series per variable each step for no effect.

`norkyst.jl` reads the operational daily-mean collection (from 2017-02-20); `norkyst_hindcast.jl` the
Norkyst-v3 hindcast (from 2012), keeping hour 12 of its 24 hourly records to match the operational
12:00 stamp. Both share everything through `AbstractNorKystConfig`.

## 9. Initial state (`src/InitialState/`)

**Why it is not the forcing.** The forcing file has to span the whole run window, because the run
reads it as forcing terms at every step and `build_simulation` checks its coverage, and the NorKyst
forcing downloads whole years. A realistic start needs one record. So the initial state is its own
config, two steps and a one-record file, and a setup starting from it needs no forcing config.

- **The date is the simulation config's `start_date`**, not a field of its own, so the two cannot
  disagree. The record is picked when the file is prepared, so the file stores the `start_date` it was
  prepared for and `FromInitialState` refuses a mismatch. A changed `start_date` therefore fails loudly
  instead of starting from the old day.
- The record is whatever the source holds for that instant: NorKyst's daily mean for the day
  `start_date` falls in, stamped 12:00, so a run starting at midnight starts from a mean centred half
  a day later.
- The file has the forcing file's layout (`define_forcing_dimensions!`, `forcing_dimension_names`,
  `NaN` on land) with one record and no `_lambda` twins, so it is read by the same `state_variables`
  as `FromForcing`. It regrids on the CPU: one record is cheap.
- The NorKyst adapter downloads with the `NorKyst` module's subsetting, as the forcing's does. Only
  the operational collection (from 2017-02-20) has an adapter.

## 10. Rivers (`src/Rivers/`)

### Rivers as a volume source

A river is water entering the fjord, not a property of the cells it enters. Relaxing T and S towards
river values diluted the surface but added no water: no sea-level rise, no estuarine outflow, no
compensating deep inflow over the sills, so basins were diluted from above and never renewed (D-2 at
50 m fell from 32 to 27.7 g/kg against an observed 31.0–31.4).

- **The only volume hook is the surface.** Oceananigans takes `forcing = (; η = F)` (m/s); `w` is
  diagnosed from the floor with no per-cell source, and an immersed face carries no normal velocity.
  A ROMS-style level or face source would need type piracy. NumericalEarth uses the same hook.
- **The tracer correction.** The inflow enters through the surface face carrying the top cell's own
  value, so `RiverTracerSource` adds `w (c_r − c) Az / V` in the top cell and the column gains
  exactly `Q c_r`. Every model tracer gets it, and a river must carry its own value: one carrying the
  water it enters runs away under upwind-biased WENO (59 °C beside Numedalslågen before a crash). So
  `add_rivers` refuses non-finite values, an unnamed tracer enters at zero, and an ungauged NVE river
  gets the gauged rivers' mean temperature. The "river budget" test closes volume and salt on a
  closed basin.
- **Columns are spread** so no column's inflow exceeds `maximum_inflow_velocity`: one 193 m column
  taking Glomma's eastern arm would be `w ≈ 1.2 cm/s`, a vertical Courant number of 2 in a 1 m top
  cell at 180 s. Grown through edge-connected water (not a radius, which crosses peninsulas), largest
  river first, one river per column.
- **A mouth moves up its river bed.** The land mask often carries the lower river as water (both
  Glomma arms run 16–22 cells up), while NVE puts the outlet at the seaward end; entering there leaves
  the bed a stagnant dead end. `narrow_water` finds channels by morphological opening and
  `river_bed_head` moves the mouth to the channel's head.
- **A file of its own**, `(Nx, Ny, time)` and daily (the records are daily, matched by calendar
  date). Everything about *where* a river enters is decided once, against the model's own land mask.
- `RiverInflow` reads nothing but its own series: the split-explicit substep calls `forcing.η` with
  only `(; η, U, V)` for fields.
- `RiverMouthMixing` (opt-in) adds NumericalEarth's river-mouth diffusivity over the river columns,
  since CATKE does not mix down a 1 m fresh layer over salt water.

### NVE: HydAPI says what, the map services say where

- **HydAPI status codes each mean one thing** (404 no such series, 401 bad key with an empty body,
  400 validation errors, 429 rate limit), so `nve_request` translates them.
- **"No data" can be HTTP 200** with `observationCount: 0`; `download_rivers` checks
  `nve_observation_count`, or an out-of-record year would be cached as a row of `NaN`. A river with
  genuinely no temperature says `temperature_station = ""`.
- **Daily records are stamped 11:00 UTC**, so series are matched by calendar date.
- **A gauge is never the mouth** — it sits tens of km upstream — and HydAPI's outlet fields are null
  or placeholders. Discharge and temperature often come from different stations (Glomma: Solbergfoss
  and Sarpfossen).
- **Discharge precedence**: the outlet's `mean_discharge`, the daily series, the REGINE catchment
  normal, `annualRunoff`, else an error — a river is a volume and has no neutral size. The catchment
  normal is what lets most discovered rivers need no gauge.
- No salinity parameter exists, hence `S = 0`. Water speed (9 stations) and conductivity
  (16 stations) are not used.
- `/Stations?Polygon=…` takes **latitude first**; `Active` is an inverted enum; a decommissioned
  station is absent from `/Stations` but present in `/Series`.
- Sanity-check small-stream temperature: Akerselva's 2020 series averages 21.2 °C and peaks at
  31.1 °C, a dry or sun-exposed sensor.

### Outlet discovery

With a finite `minimum_discharge`, `nve_outlets` derives every NVE river mouth in the domain from the
ELVIS network, and `outlets` become overrides keyed by `vassdragsnr`. Hand-copied coordinates were
wrong in ways nothing could detect: Mosseleva 7 km off, Drammenselva 602 m from water, Glomma's
western arm missing.

- **ELVIS segments are directed downstream and share endpoints exactly**, so a mouth is an end-vertex
  that starts no segment.
- **`margin` matters**: a spatial query returns whole features, so a segment crossing the box edge
  without its downstream neighbour looks like a mouth. Query wider, filter back to the domain.
- **One mouth per catchment** (highest Strahler order). Glomma's western arm survives because NVE
  files it under its own catchment.
- **The land mask is not tested here**; `river_mouths` already drops an outlet with no coastal water
  in reach.
- **An unmatched override is an error**, as is a duplicate: a typo would otherwise silently put a
  river on its catchment normal.
- The caches record the domain they were fetched for, so a grid change re-downloads.

## 11. Boundary conditions (`src/BoundaryConditions.jl`)

A piece (`AbstractBoundaryConditionConfig`, hook `boundary_condition_sides`) says what one part of
the boundary is; a set (`AbstractBoundaryConditionSetConfig`, hook `field_boundary_conditions`) says
how pieces combine. `MergedBoundaryConditions` merges with `recursive_merge`, last wins. The two
hooks have different names because they return different shapes, and the per-piece hook is not
`boundary_conditions` because that collides with `Oceananigans.Fields.boundary_conditions`.

The six-argument `boundary_condition_sides` (with the prepared `boundaries`) forwards to the
five-argument form, so a piece that carries no exterior state implements the shorter signature.

### `QuadraticBottomDrag` needs both halves

A `bottom` condition acts at the underlying grid's `k = 1`, which is immersed everywhere because
`z_faces` reaches below the deepest sounding — so the flux is discarded. The fjord's real seabed is
the *immersed* boundary, free-slip by default. With no immersed drag there was no momentum sink below
the surface: a circulation near Horten ratcheted to a 2.84 m/s jet over four days and kept
accelerating after the wind dropped. `coupled_simulation` builds `HydrostaticFreeSurfaceModel`
directly, so it inherits neither half from NumericalEarth's `ocean_simulation`. Immersed side walls
stay free-slip.

### Air-sea fluxes

Only `T` and `S` get a top condition: NumericalEarth assembles heat and salt fluxes only, so a
biogeochemical tracer's surface flux belongs in its own config. The two go through NumericalEarth's
`build_tracer_top_bc` with a `FreshwaterExchange`, because `net_fluxes` reads the exchange back out of
the condition — a bare flux condition is a `MethodError`. Its volume flux stays zero on a static-z
grid.

## 12. Boundaries (`src/Boundaries/`)

**A separate file because it is hourly.** A Flather boundary compares the model's `η` with an
exterior one, and the daily-mean collection has no tide. Only the boundary rows need that cadence:
a few hundred MB where an hourly interior file would be hundreds of GB.

- **One file for every open edge**, variables named for their side (`south_T`), on one time axis,
  with an `open_edge` attribute listing the sides; `validate_boundary_file_edges` tells a setup that
  opened another edge to re-run the step.
- **It reuses the forcing machinery on a one-cell-thick slab**: a south edge's mask is `(Nx, 1, Nz)`.
  The three surface variables get a single-level source grid centred on `SURFACE_TARGET_DEPTH`.
- **`fill_boundary_gaps!` fills dry boundary cells on read**: the schemes do arithmetic on the
  values, and on a south edge the `u` faces at `i = 1` and `Nx + 1` are peripheral (holes in the file)
  yet not immersed, so a `NaN` would be nudged into the halo.
- **The read side returns reduced `FieldTimeSeries`** keyed by edge, then by bare name. Their
  locations are exactly what `getbc` accepts, so `OpenLateralBoundaryFromData` passes them straight
  through. `BoundaryBackend` is separate from the forcing backend to keep a rank branch out of that
  hot path.
- **A large chunk cache.** The downloads' chunks span 180 hours, so under the default cache every
  hourly slab decompressed ~106 MB again: 2.4 h for 640 days. `open_source_file` gives each variable
  256 MB, bringing it to 0.7 h.
- **The operational collection's `ubar`/`vbar` are grid-relative**, unlike its derotated
  `u_eastward`/`v_northward`. With NorKyst rotated ~59°, `vbar` here is nearly *eastward*: unrotated,
  `cor(vbar, depth-mean v)` was −0.034, and a run grew a ~0.5 m/s spurious band inside the southern
  boundary. `boundary_source_slab` rotates the pair with NorKyst's `angle` **on the native grid
  before interpolation** (afterwards the components sit on different rows under different masks).
  It is a hook because the pipeline prepares one variable at a time and a rotation needs both. The
  hindcast's `ubar_eastward`/`vbar_northward` are already rotated, so
  `norkyst_hindcast_boundaries.jl` needs no such method; it reads the barotropic pair from the
  `sdepth` files and everything else from `zdepth`.

### `OpenLateralBoundaryFromData`

Per open edge, four groups, each a generic function with one method per `Val{edge}`:

| field (south edge) | condition | scheme |
|---|---|---|
| `v` | the `v` series | `NormalFlowBoundaryCondition`, `NormalRadiation` |
| `V` | `south_boundary_transport` | `GravityWaveRadiationBoundaryCondition` (Flather) |
| every model tracer | its series | `ValueBoundaryCondition`, `NormalRadiation` |
| `u` | the `u` series | `ValueBoundaryCondition`, `NormalRadiation` |

`NormalRadiation` is Orlanski with Marchesiello nudging; `GravityWaveRadiation` is what actually
opens the boundary. NorKyst's own attributes name the same set (`Che`, `Shc`, `RadNud`).

- **Do not resolve Oceananigans below 0.110.15** (PR #5351): before it, data-driven open velocity
  boundaries were broken under `SplitExplicitFreeSurface`. There is no `[compat]` entry for it.
- **`η` is not stated here**: Oceananigans pairs a Chapman condition onto every side with a
  gravity-wave `U`/`V` condition, and naming `U`/`V` switches substepping to `LocalHaloFilling`,
  which is what makes Flather act at all.
- **Edges merge as different leaves of one nested tuple.** A corner cell is written by two edges'
  conditions in Oceananigans' halo-fill order — the part of a multi-edge run to watch.
- **The transport uses the model's own column depth** (`column_depthᶜᶠᵃ`, the scheme's operator), so
  it is in the model's bathymetry, not NorKyst's. No land guard: `fill_boundary_gaps!` made the series
  finite.
- `inflow_timescale` and `outflow_timescale` are separate fields with no defaults.

## 13. Grids (`src/Grids.jl`)

Two hooks because they are needed at different moments: `domain_grid` (bare geometry, for the prepare
steps, which have no bathymetry file yet) and `simulation_grid` (with the bathymetry immersed). They
are FjordSim generics rather than methods on `LatitudeLongitudeGrid`/`ImmersedBoundaryGrid`, whose
names state a result type a rectilinear or curvilinear config could not honour. Not named
`target_grid`, which is a parameter name throughout the pipelines and would be shadowed. The loader
still accepts legacy bathymetry files with positive depths or swapped axes.

## 14. Simulations (`src/Simulations.jl`)

`coupled_simulation` lives here, not in `src/FjordSim.jl`, because a function defined after the
`include` block cannot be imported by a module inside it.

**No field anywhere in this stack has a default**: a defaulted closure would be one fjord's physics
quietly applied to another's.

**Not fields, because stated elsewhere**: the forcing, river and boundary files and the atmosphere
all come from their own configs. **`start_date` is a field because it cannot be derived**: each file
has its own first record, and when each reader zeroed its own axis the forcing (12:00) and the
atmosphere (00:00) ran twelve hours out of phase. `build_simulation` passes it to every reader as
`reference_date`; a `start_date` before a file's first record is an error, not a shift.

**`validate_time_coverage`** checks every prepared file spans `[start_date, start_date + stop_time]`,
which is what makes `Cyclical()`'s wrap unreachable. One `stop_time`, not `loops` of them — every
loop replays the same window. Each check goes through a hook on its own config
(`forcing_date_range`, `atmosphere_date_range`, …), so the module names no file layout.

`model_tracers(model)` is a hook because the tracer list is needed before the model exists: the
forcing terms, the open tracer boundaries and the initial conditions each build one thing per tracer.

### `extra_kwargs`

`CoupledHydrostaticSimulation` stores nine grid-independent model components plus `extra_kwargs`:
one `NamedTuple` per constructor `coupled_simulation` calls (`ocean_model`, `ocean_simulation`,
`coupled_simulation`, `coupled_model`), splatted after the explicit keywords. Ten fields is the
ceiling; anything further goes in `extra_kwargs`. Because a duplicate key would silently win,
`validate_extra_kwargs` rejects a key `coupled_simulation` already passes (`EXTRA_KWARG_SLOTS`), at
construction, where the error still points at the setup file.

`free_surface` is a config with its own hook because `SplitExplicitFreeSurface` needs the grid.

### `closure` and `BoundarySponge`

`closure` holds a built Oceananigans closure or an `AbstractClosureConfig`, resolved by
`model_closure(closure, grid, boundary_config, …)` whose fallback is the identity.

`BoundarySponge` appends a `HorizontalScalarDiffusivity` ramping to zero over `width_cells` from every
open edge. It exists because an open boundary radiates and nothing absorbed it: on the 2020
`oslofjorden` run the boundary row's velocities were 3–7x the interior's, and after fifty days
near-boundary bottom salinity ran away (54 → 192 psu) while the domain mean moved 1 psu and the
boundary file asked for 24–34.

- **Viscous, not a tracer band**: viscosity names no target and cannot fight the boundary nudging.
- **The edges come from `open_edges(boundary_config)`**; a setup naming none gets `base` back.
- **Stability bounds the coefficient**: `Δt ≤ Δx²/4ν` is 310 s at `ν = 30` and 155 s at `ν = 60` on
  a 193 m cell. Raise `ν` and check `max_time_step`.
- **The ramp is `max(0, 1 − d/w)²`**, smooth at its inner end and compactly supported.
- **It must be `discrete_form = true`.** The continuous form drops `parameters` and calls
  `ν(λ, φ, z, t)`; on GPU the resulting `MethodError` surfaces as an `InvalidIRError` naming an
  unrelated momentum-kernel frame. The "boundary sponge" testset calls the coefficient through
  Oceananigans' `νᶜᶜᶜ`/`κᶜᶜᶜ` so that is a test failure instead.

### Initial conditions

`initial_conditions` is a `NamedTuple` (passed to `set!` as is), `FromForcing(date)`,
`FromInitialState()` or `FromResults(path, date)`, resolved by `resolve_initial_conditions`. Every
file is already on the model grid, so this is a read and a type conversion, not an interpolation.

- Forcing-file and initial-state land is `NaN` and is zeroed — exact, because the prepare step wrote
  `NaN` at exactly the cells `peripheral_node` calls dry.
- Arrays are converted to `eltype(grid)`: a `Union{Missing,Float32}` array cannot become a `CuArray`.
- `FromForcing` needs the date exactly on the axis. `FromResults` takes the nearest record within
  half the writer's `interval`, because an adaptive step writes each record late (up to 35 s on a
  640-day run). A snapshot's `time` has no calendar, so the writer records `start_date` as a global
  attribute.
- What is set is `(tracers ∪ {u, v}) ∩ variables in the file` (`state_variables`). `η`, CATKE's `e`
  and the tendencies keep their defaults, so all three are a **lossy warm start**; `pickup` is the exact
  restart. Do not conflate them.
- `scripts/extract_initial_state.jl` copies one snapshot record into a small file, because a cloud
  run stages prepared inputs, not an earlier run's tens-of-GB snapshots.

### Looping

`loops` replays the window with the state carried over (a spin-up), one output file per loop.
`restart_loop!` resets every clock and nothing else. Four traps, all silent:

- **NumericalEarth's `reset_clock!` cannot be used**: it reads `sea_ice.clock`, which
  `FreezingLimitedOceanTemperature` lacks. `rewind_clock!` dispatches on `hasfield(…, :clock)`; the
  "loop restart clocks" testset asserts the upstream method still throws.
- **The atmosphere and radiation need `update_state!`**, or the first step of each loop is forced by
  the end of the window.
- **`ocean_sim.initialized` must be cleared by hand**, or the new writer never writes its `t = 0`
  record.
- **The writer is replaced and the old one closed**: a `TimeInterval` keeps its actuation count, and a
  `NetCDFWriter` holds an open file.

A monotonic clock over `loops * stop_time` was rejected: `Cyclical`'s inferred period stops matching
the loop once a file is padded, and fixing it would change the atmosphere hook signatures.

### Writers and callbacks

`validate_writers` rejects two checkpointing writers and two writers under one key;
`validate_callbacks` rejects two callbacks under one name — in each case the second silently replaces
the first.

- **`η` cannot go through Oceananigans' NetCDF writer**: a `(Center, Center, Nothing)` user output
  asks for a `z_aaf` the grid already owns, and `create_spatial_dimensions!` raises.
  `FieldSnapshotWriter` and `FieldStationWriter` write the same thing to JLD2. The split is by what a
  format can hold.
- **`StationWriter`** writes whole-column profiles at named positions, one file per station
  (Oceananigans takes one `indices` per file), at the cadence validation needs: 150 MB of station
  series beside 59 GB of daily 3D fields on the Oslofjord validation run. A dry station cell searches
  `search_radius` cells for water and logs how far it moved; a station with none is dropped with a
  warning.
- A name the model lacks is an **error** in a writer but intersected away in `state_variables`: the
  setup wrote the writer's names down, while the reader serves two kinds of file.
- `ProgressCallback`'s `report` is a field because `Utils.progress` reads `tracers.T`, and a model
  without `T` would crash at the first report, after compiling.
- `mkpath`: `build_simulation` creates `results_root`, the snapshot writer its own file's directory —
  `NetCDFWriter` creates only a `dir` keyword, and FjordSim passes the whole path as `filename`.

### Checkpointing

A `CheckpointWriter` goes on the **coupled** simulation: its `prognostic_state` is the resumable
state, and `run!(…; pickup)` looks there for exactly one checkpointer. Only prognostic state is
serialized, so no `ForcingFromFile` or `FieldTimeSeries` backend is ever written; expect a few hundred
MB per checkpoint on the Oslofjord grid, hence `cleanup = true`.

- The prefix carries the **loop index** (`checkpoint_prefix`), which `resume_loop` reads back, or a
  pickup could not tell loop 3 from loop 1.
- It carries **no run tag**, or a later launch could not find the previous one's checkpoints. So
  there is **one resumable run per results directory**: `Checkpointer` writes with `"w"` whatever
  `overwrite_files` says, and `pickup` resumes whatever is there, even from another `start_date`.
- Snapshots *are* per launch, so a resumed run's file starts at the resume point.
- `pickup` supersedes `initial_conditions`.

## 15. Validation (`src/Validation/`)

The only module that runs after a simulation, built the same way as the others: an
`AbstractObservationConfig` per source. Every stage is permissive about missing data — most
station-variable combinations in a sweep are empty, and each is a blank row, not a failure.
**Stations are matched by name** between a source and the run's station writers;
`default_observations` builds its station list from the writers, and a mismatch shows up as "no model
output", not an error.

## 16. CLI (`src/CLI.jl`)

`parse_arguments` never `exit`s, so the help path is testable; exit codes are 0 success, 1 step
failed, 2 bad arguments.

**Only `run_simulation` is logged**, by `tee_output` to `fjordsim_<run_tag>.log` under
`results_root`: its failure is a stacktrace through the whole coupled model, its output is not
reproducible by re-running, and it is the only step guaranteed a `results_root`.

- `redirect_stdio` takes only fd-backed streams, so the tee is a `Pipe` (needing
  `Base.link_pipe!`) plus a copying task; output may arrive in bursts but is never lost.
- **The failure is caught inside the redirect** by `run_step`, or `Base._start` would print it after
  the redirect was torn down and the error would be missing from the log. Parsing stays outside, so a
  usage error leaves no log.
- `show_compact_error` sets `:stacktrace_types_limited`, which is what abbreviates type parameters
  to `{…}`. `STACKTRACE_WIDTH` below 120 changes nothing (`Base.type_depth_limit` floors it).
- During a run `stdout` is a `Pipe`, so `displaysize` is 24x80 and wide `show` output may wrap.

## 17. Top-level (`src/FjordSim.jl`)

**It overrides one NumericalEarth method**, `update_net_ocean_fluxes!` on `NoSeaIceInterfaceModel`,
with a copy of upstream's body. Without sea ice, NumericalEarth's zero fluxes are an untyped
`ifelse` branch that LLVM's NVPTX back-end segfaults compiling, so a GPU run dies before its first
step. The copy warns at precompile when NumericalEarth is not the version it was written against;
delete it once upstream types those zeros.

`main` is **not exported**: Julia runs `Main.main` after a script whenever it resolves to an entry
point, so exporting it would run the CLI at the end of every script that does `using FjordSim`. The
`@main` is bare and comes *after* the definition: `@main function main(args) … end` expands to a call
and would run the CLI during precompilation.
