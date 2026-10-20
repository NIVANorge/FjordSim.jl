# CLAUDE.md

> This file holds what's needed on every task: commands, hard rules, and gotchas. Everything else
> lives in the docs below and is read on demand.

## Docs — read before the matching task

| Doc | Read it before |
|---|---|
| [`docs/running.md`](docs/running.md) | running any pipeline step; changing `start_date`, `stop_time`, `z_faces`, the grid `size` or a bathymetry knob; a fresh rerun into a used results directory; resuming |
| [`gcp/README.md`](gcp/README.md) | anything with `gcp/fjordsim-gcp` or the cloud VM |
| [`docs/architecture.md`](docs/architecture.md) | changing a module in `src/` — its section only (design rationale, ordering constraints, measured numbers) |
| [`docs/adding-a-source.md`](docs/adding-a-source.md) | adding a grid, bathymetry, forcing, river, boundary, atmosphere or observation source |
| [`docs/setups.md`](docs/setups.md) | adding a fjord or tuning a setup's parameters |
| [`README.md`](README.md) | user-facing overview, quick start, `run_simulation` fields |
| [`TODO.md`](TODO.md) | open work and the measurements behind it |

The `run-simulation` skill (`.claude/skills/run-simulation/`) is the checklist for running steps.

## Commands

```bash
julia --project test/runtests.jl                                   # all tests
julia --project -m FjordSim --help                                 # subcommands and setups
julia --project -m FjordSim <step> --config <setup-or-path.jl>     # one pipeline step
```

The steps, in order: `prepare_bathymetry`, `download_forcing`, `prepare_forcing`, `add_rivers`
(needs `NVE_API_KEY`), `download_boundaries`, `prepare_boundaries`, `download_atmosphere`,
`prepare_atmosphere`, `run_simulation` (needs a GPU), `validate_simulation`. `--config` is the only
option and is required. From the REPL each step is a function of the same name on a `FjordConfig`;
`build_simulation(config)` assembles a `Simulation` without running it. Details: `docs/running.md`.

`run_simulation` on a cloud GPU goes through `gcp/fjordsim-gcp`; read `gcp/README.md` first.

`-m` is Julia 1.12's package entry point, which is why `Project.toml` has `[compat] julia = "1.12"`.

## Hard rules for runs

- **One run per results directory.** Checkpoints carry no run tag. A different experiment gets its
  own `results_root`; a from-scratch rerun of the same config first archives the old run as
  `<name>_<YYYY-MM-DD>`. Procedure in `docs/running.md`.
- **Prepared files go stale silently.** Changing `start_date`/`stop_time` or `z_faces`/`size`/a
  bathymetry knob does nothing until the listed prepare steps are re-run — `simulation_grid` reads
  `bathymetry.nc`, not the config. Lists in `docs/running.md`.

## Architecture

FjordSim wraps [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) and
[NumericalEarth.jl](https://github.com/NumericalEarth/NumericalEarth.jl) to set up regional ocean
simulations of Norwegian fjords. A simulation reads prepared files only. The bathymetry (which
defines the grid) is the one every setup needs; the forcing, rivers, open-boundary data and
atmosphere are each optional — a `FjordConfig` field defaulting to `nothing`, whose steps are then
no-ops (`docs/architecture.md`, "What a run is assembled from").

**A new grid, bathymetry source, forcing dataset, river dataset, open-boundary dataset or
atmosphere dataset is added by subtyping the matching `Abstract*Config` supertype and overloading
methods on it — never by editing the generic pipeline.** Read `docs/adding-a-source.md` first.

The modules, in `include` order from `src/FjordSim.jl`:

1. **Configs** (`src/Configs.jl`) — the abstract supertypes, `FjordConfig`, the edge vocabulary
   (`LATERAL_EDGES`, `open_edges`), and the path helpers every source inherits.
2. **Utils** (`src/Utils.jl`) — `progress` callback, `recursive_merge`, time-step helpers.
3. **Plotting** (`src/Plotting.jl`) — one `plot_*` function per pipeline, dispatching on config
   supertypes.
4. **Bathymetry** (`src/Bathymetry/`) — `prepare_bathymetry` and the eight-stage
   `smooth_bathymetry_gaps!` pipeline (order is forced; see the doc before reordering).
5. **Atmosphere** (`src/Atmospheres/`) — `prepare_atmosphere`, NORA3 read/download.
6. **Regridding** (`src/Regridding.jl`) — the core `Forcing` and `Boundaries` share: source
   records and time axes, water masks, the projected source grid and the interpolation kernel.
7. **Forcing** (`src/Forcing/`) — the interior forcing pipeline (`prepare_forcing`,
   `download_forcing`) and its NorKyst adapters.
8. **Rivers** (`src/Rivers/`) — the river pipeline (`add_rivers`) and its NVE adapter.
9. **Boundaries** (`src/Boundaries/`) — the open-boundary data pipeline (`download_boundaries`,
   `prepare_boundaries`, `boundary_series`) and its NorKyst adapters.
10. **Boundary conditions** (`src/BoundaryConditions.jl`) — `AirSeaFluxes`, `QuadraticBottomDrag`,
    `OpenLateralBoundaryFromData`.
11. **Grids** (`src/Grids.jl`) — `EvenGrid`, the `domain_grid`/`simulation_grid` hooks.
12. **Simulations** (`src/Simulations.jl`) — `SimulationConfig`, `build_simulation`/`run_simulation`,
    writers, callbacks, looping, checkpointing, initial conditions.
13. **Setups** (`src/Setups/`) — the built-in fjords; see `docs/setups.md`.
14. **Validation** (`src/Validation/`) — the only module that runs *after* a simulation:
    observation sources, skill and tidal statistics, and the comparison figures. `StationWriter`
    (in `Simulations`) is what a run writes for it.
15. **CLI** (`src/CLI.jl`) — subcommand dispatch, logging.
16. **Top-level** (`src/FjordSim.jl`) — re-exports and `main`.

## Setups

`src/Setups/` holds one lowercase file per fjord, each a zero-arg function returning a fresh
`FjordConfig`, registered in `SETUPS`. Adding a fjord is a two-place edit: the new file, and its
entry in `SETUPS`. `SimulationConfig` and everything it nests have **no defaults** — every field is
a scientific choice a setup must state.

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
- Keep imports explicit so the dependency surface stays auditable

## Tests

- **Tests are general, never setup-specific.** Never assert a value a setup states — a grid size, a
  box, a year, a `start_date`, a station list, which dataset it picks. Those are the setup author's
  to change, and a test that pins them turns every scientific choice into a test failure. Assert the
  *conventions* a setup getting it wrong would break, over `setup_names()`, the way the "Setups"
  testset already does.
- Build fixtures from `test/utilities.jl` (`test_simulation_config`, `immersed_test_grid`, ...)
  rather than reading a registered setup with `fjord_config`.
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
  model clock stays `Float64` whatever the grid is (`coupled_simulation`): a `Float32` kernel time
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
  so `SimulationConfig(stop_time = 3600, ...)` is a `MethodError` where `stop_time = 1hour` is
  fine. Every `Oceananigans.Units` constant is already a `Float64`, so writing durations with them —
  `365days`, `1hour`, `3minutes` — sidesteps this. `loops` is the one field this does *not* apply to,
  being an `Int` already. `SimulationConfig` is now the **only** config the rule reaches: every other
  one is either non-parametric, or — like `SnapshotWriter`, `ProgressCallback` and
  `CoupledHydrostaticSimulation`, all of which *are* parametric — has a hand-written keyword
  constructor that converts. `Regridding.PreparedVariable` is the one internal type the rule reaches
  too, for the same reason.

## Key conventions

- Bathymetry convention: `h < 0` = below sea level (bottom height), `h >= 0` = land.
- Data files default to `~/FjordSim_data/<fjord>/` and results to `~/FjordSim_results/<fjord>/`,
  the latter from the simulation config's `results_root`. A setup must build both with
  `fjord_data_root("<fjord>")` / `fjord_results_root("<fjord>")` (`src/Configs.jl`) rather than
  `joinpath(homedir(), ...)`: those read `FJORDSIM_DATA_ROOT` / `FJORDSIM_RESULTS_ROOT`, which is
  what lets one setup run unchanged on a laptop and on a cloud VM's staging disk. A setup sharing
  another fjord's downloads passes that other name, so a relocation moves both.
- Open-boundary convention: the domain is open on any subset of its four lateral edges — none, one,
  or all four for a region in the open ocean — named once by the boundary data config's `open_edges`
  and read through the `open_edges` accessor as a `Vector{Symbol}`, empty for a setup naming no
  boundary config, whose lateral boundaries are then all closed walls. Everything that acts on an edge
  dispatches on `Val(edge)` and every consumer *iterates*, so the four edges are independent; the
  prepared boundary file names every variable for its side (`south_T`), which is how they all fit in
  one file on one time axis.
