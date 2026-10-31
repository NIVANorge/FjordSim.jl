# CLAUDE.md

> This file holds what's needed on every task: commands, hard rules, and gotchas. Everything else
> lives in the docs below and is read on demand.

## Docs — read before the matching task

| Doc | Read it before |
|---|---|
| [`docs/running.md`](docs/running.md) | running any pipeline step; changing `start_date`, `stop_time`, `z_faces`, the grid `size` or a bathymetry knob; a fresh rerun into a used results directory; resuming |
| [`gcp/README.md`](gcp/README.md) | anything with `gcp/fjordsim-gcp` or the cloud VM |
| [`docs/architecture.md`](docs/architecture.md) | changing a module in `src/` — its section only (design rationale, ordering constraints, measured numbers) |
| [`docs/dispatch-map.md`](docs/dispatch-map.md) | looking up which hooks a supertype has and who implements them; it is generated — after changing a hook, a supertype docstring or a step, rerun `julia --project docs/make_dispatch_map.jl` |
| [`docs/adding-a-source.md`](docs/adding-a-source.md) | adding a grid, bathymetry, forcing, river, boundary, atmosphere or observation source |
| [`docs/setups.md`](docs/setups.md) | writing a new setup or changing a setup's parameters |
| [`README.md`](README.md) | user-facing overview, quick start, `run_simulation` fields |
| [`TODO.md`](TODO.md) | open work and the measurements behind it |

The `run-simulation` skill (`.claude/skills/run-simulation/`) is the checklist for running steps,
and the `new-setup` skill (`.claude/skills/new-setup/`) the one for writing a new setup.

## Commands

```bash
julia --project test/runtests.jl                                   # all tests
julia --project -m FjordSim --help                                 # subcommands
julia --project -m FjordSim <step> --config <setup.jl>             # one pipeline step
```

The steps, in order: `prepare_bathymetry`, `download_forcing`, `prepare_forcing`,
`download_initial_state`, `prepare_initial_state`, `add_rivers` (needs `NVE_API_KEY`),
`download_boundary_data`, `prepare_boundary_data`, `download_atmosphere`, `prepare_atmosphere`,
`run_simulation` (needs a GPU), `validate_simulation`. `--config` is the only
option and is required. From the REPL each step is a function of the same name on a setup (any
`AbstractFjordConfig`); `build_simulation(config)` assembles a `Simulation` without running it. Details: `docs/running.md`.

`run_simulation` on a cloud GPU goes through `gcp/fjordsim-gcp`; read `gcp/README.md` first.

`-m` is Julia 1.12's package entry point, which is why `Project.toml` has `[compat] julia = "1.12"`.

## Hard rules for runs

- **One run per results directory.** Checkpoints carry no run tag. A finished run is archived as
  `<name>_<YYYY-MM-DD>`, in the bucket and locally, and the next run — rerun, test or changed
  parameter — uses the same config file, edited, not a new one. Procedure in `docs/running.md`.
- **Prepared files go stale silently.** Changing `start_date`/`stop_time` or `z_faces`/`size`/a
  bathymetry knob does nothing until the listed prepare steps are re-run — `simulation_grid` reads
  `bathymetry.nc`, not the config. Lists in `docs/running.md`.

## Branches

- **Never change code on `main`.** Before the first edit for a task, check the current branch
  (`git branch --show-current`). If it is `main`, create a branch named for the task
  (e.g. `river-mouth-salinity`) and switch to it with `git switch -c <name>`; every change for the
  task goes there.
- When the work is ready, ask before merging the branch back into `main` — never merge on your own.
- Once the merge has succeeded, delete the merged branch (`git branch -d <name>`).

## Architecture

FjordSim wraps [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) and
[NumericalEarth.jl](https://github.com/NumericalEarth/NumericalEarth.jl) to set up regional ocean
simulations of Norwegian fjords. A simulation reads prepared files only. The bathymetry (which
defines the grid) is the one every setup needs; the forcing, initial state, rivers, open-boundary
data and atmosphere are each optional — a `FjordConfig` field defaulting to `nothing`, whose steps are then
no-ops (`docs/architecture.md`, "What a run is assembled from").

**A new grid, bathymetry source, forcing dataset, initial-state dataset, river dataset,
open-boundary dataset or atmosphere dataset is added by subtyping the matching `Abstract*Config` supertype and overloading
methods on it — never by editing the generic pipeline.** Read `docs/adding-a-source.md` first.

**A step has one method, on the setup type, and holds no logic.** Every pipeline step and
`build_simulation` is defined only as `step(config::AbstractFjordConfig)`, which picks out the
components it reads through the accessors (`grid_config(config)`, …) and calls a work function of
another name (`regrid_forcing`, `place_rivers`, …) that dispatches on those components, typing only
its own; a missing component is a `::Nothing` method, not an `isnothing` branch. Never add a method
of a step's own name on anything but a setup type. The model (`model_config`), its biogeochemistry,
its boundary conditions and the run (`run_config`, read only through its accessors) are separate
setup components (`docs/architecture.md`, Configs).

The modules, in `include` order from `src/FjordSim.jl`:

1. **Configs** (`src/Configs.jl`) — the abstract supertypes, `AbstractFjordConfig` and its
   accessors, and `FjordConfig`.
2. **Utils** (`src/Utils.jl`) — `progress` callback, `recursive_merge`, time-step helpers,
   `coverage_window`, `fjord_data_root`/`fjord_results_root`.
3. **Grids** (`src/Grids.jl`) — `EvenGrid`, the `domain_grid`/`simulation_grid` hooks, and the edge
   vocabulary (`LATERAL_EDGES`, `open_edges`).
4. **Plotting** (`src/Plotting.jl`) — one `plot_*` function per pipeline, dispatching on config
   supertypes, and the `plot_path` generic each pipeline adds a method to.
5. **Bathymetry** (`src/Bathymetry/`) — `prepare_bathymetry` and the eight-stage
   `smooth_bathymetry_gaps!` pipeline (order is forced; see the doc before reordering).
6. **Regridding** (`src/Regridding.jl`) — the core `Forcing`, `InitialState` and `BoundaryData`
   share: source records and time axes, water masks, the projected source grid and its rotation,
   the interpolation kernel, and the prepared files' layout.
7. **Atmosphere** (`src/Atmospheres/`) — `prepare_atmosphere`, NORA3 read/download.
8. **NorKyst** (`src/DataSources/NorKyst.jl`) — the NorKyst-800m dataset the forcing, initial-state
   and boundary adapters share: endpoints, variable names, source grid, subset download.
9. **Forcing** (`src/Forcing/`) — the interior forcing pipeline (`prepare_forcing`,
   `download_forcing`) and its NorKyst adapters.
10. **Initial state** (`src/InitialState/`) — `download_initial_state`, `prepare_initial_state`:
    one source record at `start_date` for `FromInitialState`, and its NorKyst adapter.
11. **Rivers** (`src/Rivers/`) — the river pipeline (`add_rivers`) and its NVE adapter.
12. **Boundary conditions** (`src/Boundaries/BoundaryConditions/`) — the `boundary_condition_sides`
    hook, `model_boundary_conditions` (the default merge on any model config), `AirSeaFluxes`,
    `QuadraticBottomDrag`, and `OpenLateralBoundaryFromData`, the condition that reads the prepared
    open-boundary data.
13. **Boundary data** (`src/Boundaries/BoundaryData/`) — the open-boundary data pipeline
    (`download_boundary_data`, `prepare_boundary_data`, `boundary_series`) and its NorKyst adapters.
    Independent of `BoundaryConditions`; `build_simulation` hands one's output to the other.
14. **Advection schemes** (`src/AdvectionSchemes.jl`) — `RiverMouthAdvection` (first-order upwind
    horizontally around the river mouths) and `model_tracer_advection`, the hook `coupled_model` builds
    the tracer advection through.
15. **Simulations** (`src/Simulations.jl`) — `CoupledHydrostaticModel` and `coupled_model`,
    `RunConfig`, `build_simulation` (staged in `assemble_simulation`: `validate_model_inputs`,
    `model_inputs`, `coupled_model`)/`run_simulation` (`run_loops`), writers, callbacks, looping,
    checkpointing, initial conditions, the run tag.
16. **Validation** (`src/Validation/`) — the only module that runs *after* a simulation:
    observation sources, skill and tidal statistics, and the comparison figures. `StationWriter`
    (in `Simulations`) is what a run writes for it.
17. **CLI** (`src/CLI.jl`) — `fjord_config` (loads a config file), subcommand dispatch, logging.
18. **Top-level** (`src/FjordSim.jl`) — re-exports and `main`.

## Setups

No setup is part of the package. A setup is a script whose last expression is a `FjordConfig`,
passed by path to `--config` and loaded by `fjord_config`; `examples/` holds the worked ones
(`oslofjorden.jl`, `drammensfjorden.jl`, `lista.jl`, `moro.jl`, `mjosa.jl` — a closed lake — and
variants that `include("oslofjorden.jl")`). Adding a fjord is one new file. `CoupledHydrostaticModel`,
`RunConfig` and everything they nest have **no defaults** — every field is a scientific
choice a setup must state. How to write
one is `docs/setups.md`; why a setup has the values it has belongs in its own comments.

## GPU & Kernel Compatibility

Forcing callables (e.g. `ForcingFromFile`) and any function called inside an Oceananigans
kernel must follow these rules:

- Mark with `@inline`
- Use `ifelse` — never short-circuiting `&&`/`||` or ternary `?:` inside kernels
- Must be type-stable and allocation-free
- Never loop over grid points outside kernels — use `launch!` via KernelAbstractions
- Use literal zeros: `max(0, a)` not `max(zero(FT), a)`
- Never hardcode `Float64` literals (`0.0`, `1.0`) in kernels or constructors — use `zero(grid)`,
  `one(grid)`, `convert(FT, val)`, or rational literals like `1//2`
- Use `on_architecture` for CPU↔GPU data transfers — never call `Array()` or `CuArray()` directly
- Always index fields with three indices: `field[i, j, k]` — 2D indexing appears to work on some
  fields by coincidence but is unsupported and will break

`src/Forcing/Forcing.jl` already demonstrates the correct pattern; match it when adding new
forcings.

## Type Stability

- All structs must have concrete field types — no `Any`, no abstract field types
- User-facing constructors may accept flexible arguments, but the returned struct must be
  fully typed (follow the materialization pattern from NumericalEarth if needed)
- Type annotations are for dispatch only, not documentation
- Profile with `@code_warntype` when touching GPU-executed paths

## Import Conventions

- Extend functions via `function Mod.foo(...) ... end` — never `import Mod: foo` to extend
- Exports go at the top of each module file, before other code
- Export only what a setup names: built-in configs, writers, boundary conditions, pipeline steps.
  Abstract supertypes, hooks and internals stay unexported and are reached through their module
  (`using FjordSim.Configs: AbstractRiverConfig`, `Rivers.river_series(config::MyRivers, times)`).
  A submodule exports nothing the top level does not re-export
- Keep imports explicit so the dependency surface stays auditable

## Tests

- **Tests are general, never setup-specific.** Never assert a value a setup states — a grid size, a
  box, a year, a `start_date`, a station list, which dataset it picks. Those are the setup author's
  to change, and a test that pins them turns every scientific choice into a test failure. The suite
  never loads a file from `examples/`.
- Build fixtures from `test/utilities.jl` (`test_run_config`, `immersed_test_grid`, ...).
  A config file a test writes is evaluated in `Main`, so it can call those fixtures too.
- **Do not test plotting.** No test asserts that a figure was drawn or how it looks.
- The suite is offline: write a fixture file rather than fetching one. Anything needing the network
  is a manual integration check, not a test.

## Code Style

- Avoid abbreviations: `latitude` not `lat`, `temperature` not `temp` (exception: physics
  symbols `T`, `S`, `u`, `v`, `λ`, `φ` are conventional)
- Keyword args: no-space inline `f(x=1)`, single-space multiline `f(\n    a = 1,\n    b = 2\n)`
- No trailing whitespace, no trailing blank lines; files end with exactly one newline
- Always use explicit `return` in functions longer than one expression
- Delete commented-out code — git is the history; no debugging artifacts or stale copy-paste
- Comments are for any reader, on any machine, at any later time: state what holds for everyone —
  what the code does and why, and measured numbers with what they were measured on. Never the
  author's situation or the file's history: not "the FileGDB is already on this machine", "Lista
  has no run yet", "the 42 m layers this file had", "from your environment". A machine's state
  belongs nowhere in the repo; a task's progress goes in `TODO.md`; history goes in the commit
  message.

## Writing docs, CLAUDE.md and TODO.md

- **No user- or account-specific information in any checked-in file** — no names, email addresses,
  personal account permissions ("this account cannot …"), local paths of one machine, or one
  person's preferences. Describe what the code or workflow does ("the workflow starts the VM through
  two Cloud Run jobs"), not who can do what. Personal facts and preferences go in Claude's memory.
- **CLAUDE.md holds only what every task needs.** Procedures and rationale go in `docs/` or
  `gcp/README.md`, with a row in the docs table above; don't copy their content here.
- **Remove completed items from `TODO.md`** when finishing them — don't tick them off and leave them.
  If an open item needs a result the finished one produced, fold that result into the open item in a
  line or two; the full record of what was done goes in the commit message.

## Common Pitfalls

- Never set `Oceananigans.defaults.FloatType` in a setup: it is process-wide, so every config built
  after it — other setups, test fixtures — silently changes precision. The simulation's float type
  is `EvenGrid`'s `float_type` (default `Float32`, for GPU speed), and a setup passes that same `FT`
  to every model component that takes one (`WENO(FT)`, `CATKEVerticalDiffusivity(FT; ...)`,
  `FreezingLimitedOceanTemperature(FT)`, ...), since those read the global when not given one. The
  model clock stays `Float64` whatever the grid is (`coupled_model`): a `Float32` kernel time
  rounds to a different `FieldTimeSeries` slice than the host near an exact hour and faults on GPU.
- Never extend `getproperty` to fix undefined-property bugs — fix the caller instead
- A variable named the same as a function produces "type is not callable" — rename the variable
- Never add/remove/change `[deps]` in `Project.toml` on your own initiative; only touch
  `[compat]` when explicitly asked. But if a dependency would meaningfully simplify the
  implementation, *ask* — do not silently contort the design to avoid it. This applies to packages
  already present transitively (e.g. `KernelAbstractions` via Oceananigans), where a direct
  `[deps]` entry costs nothing but is still required to `using` them.
- `Pkg.resolve()` only *checks* the manifest against `[compat]`; it will not move a version to
  satisfy a widened bound, and reports `empty intersection between X@old and project compatibility`
  instead. Use `Pkg.update()` after raising a compat bound. (`Pkg.up` is not a name in current Pkg.)
- `Base.@kwdef` on a *parametric* struct does not convert its arguments, unlike on a
  non-parametric one: the generated constructor keeps the declared field types in its signature,
  so `RunConfig(stop_time = 3600, ...)` is a `MethodError` where `stop_time = 1hour` is
  fine. Every `Oceananigans.Units` constant is already a `Float64`, so writing durations with them —
  `365days`, `1hour`, `3minutes` — sidesteps this. `loops` is the one field this does *not* apply to,
  being an `Int` already. `RunConfig` is now the **only** config the rule reaches: every other
  one is either non-parametric, or — like `SnapshotWriter`, `ProgressCallback` and
  `CoupledHydrostaticModel`, all of which *are* parametric — has a hand-written keyword
  constructor that converts. `Regridding.PreparedVariable` is the one internal type the rule reaches
  too, for the same reason.

## Key conventions

- Bathymetry convention: `h < 0` = below sea level (bottom height), `h >= 0` = land.
- Data files default to `~/FjordSim_data/<fjord>/` and results to `~/FjordSim_results/<fjord>/`,
  the latter from the run config's `results_root`. A setup must build both with
  `fjord_data_root("<fjord>")` / `fjord_results_root("<fjord>")` (`src/Utils.jl`) rather than
  `joinpath(homedir(), ...)`: those read `FJORDSIM_DATA_ROOT` / `FJORDSIM_RESULTS_ROOT`, which is
  what lets one setup run unchanged on a laptop and on a cloud VM's staging disk. A setup sharing
  another fjord's downloads passes that other name, so a relocation moves both.
- Open-boundary convention: the domain is open on any subset of its four lateral edges — none, one,
  or all four for a region in the open ocean — named once by the boundary data config's `open_edges`
  and read through the `open_edges` accessor as a `Vector{Symbol}`, empty for a setup naming no
  boundary data config, whose lateral boundaries are then all closed walls. Everything that acts on an edge
  dispatches on `Val(edge)` and every consumer *iterates*, so the four edges are independent; the
  prepared boundary file names every variable for its side (`south_T`), which is how they all fit in
  one file on one time axis.
