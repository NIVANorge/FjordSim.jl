# FjordSim.jl

A framework for regional ocean simulations built on top of
[Oceananigans](https://github.com/CliMA/Oceananigans.jl) and
[NumericalEarth](https://github.com/NumericalEarth/NumericalEarth.jl).

A simulation is completely defined in a config file, for example
[`examples/oslofjorden.jl`](examples/oslofjorden.jl).
It contains the simulation parameters and the programming objects that define the logic of a simulation.
It is done using Julia's main feature/paradigm — multiple dispatch.
Providing new structs and overloading existing/new methods on them, a user can change anything.
At the same time, since the entire simulation is a config type file, all the parts are in one place,
compact, and readable.

## Installation

```bash
git clone https://github.com/NIVANorge/FjordSim.jl.git
cd FjordSim.jl
julia --project -e 'using Pkg; Pkg.instantiate()'
```

`instantiate` resolves a `Manifest.toml` from `Project.toml`, then installs and precompiles the dependencies.
To track the repository from a project of your own rather than working inside it,
`add https://github.com/NIVANorge/FjordSim.jl.git`.

Julia **1.12 or newer** is required: the `-m FjordSim` command form used throughout this README is
Julia 1.12's package entry point.

A GPU is not strictly required — `architecture = :auto` falls back to the CPU — but running a
simulation without one is impractically slow. Preparing the input data is fine on CPU.

Data lives under `~/FjordSim_data/<fjord>/` and results under `~/FjordSim_results/<fjord>/`.
Set `FJORDSIM_DATA_ROOT` or `FJORDSIM_RESULTS_ROOT` to change the default directories,
or set the destination in a config file.

## Quick start: Oslofjorden

[`examples/oslofjorden.jl`](examples/oslofjorden.jl) is the reference setup: Geonorge bathymetry,
an initial state from one NorKyst-800m daily mean, NVE rivers, hourly NorKyst-800m open-boundary
data along the southern edge, and NORA3 atmosphere. It names no forcing, so `download_forcing` and
`prepare_forcing` report that and do nothing.
Please note: It may take several hours to download the data.

`add_rivers` also needs a free NVE HydAPI key — sign up at https://hydapi.nve.no/Users and set it
as `NVE_API_KEY` in your environment before running that step.

Preparing the data and running it is nine commands. Only `prepare_bathymetry` is
required of every setup; the others belong to optional components (see
[What a run is made of](#what-a-run-is-made-of)):

```bash
julia --project -m FjordSim prepare_bathymetry     --config examples/oslofjorden.jl  # required
julia --project -m FjordSim download_initial_state --config examples/oslofjorden.jl  # optional: initial state
julia --project -m FjordSim prepare_initial_state  --config examples/oslofjorden.jl  # optional: initial state
julia --project -m FjordSim add_rivers             --config examples/oslofjorden.jl  # optional: rivers
julia --project -m FjordSim download_boundary_data --config examples/oslofjorden.jl  # optional: open boundaries
julia --project -m FjordSim prepare_boundary_data  --config examples/oslofjorden.jl  # optional: open boundaries
julia --project -m FjordSim download_atmosphere    --config examples/oslofjorden.jl  # optional: atmosphere
julia --project -m FjordSim prepare_atmosphere     --config examples/oslofjorden.jl  # optional: atmosphere
julia --project -m FjordSim run_simulation         --config examples/oslofjorden.jl  # the run itself
```

Each subcommand is named after the function it calls, so the same steps run from the REPL:

```julia
using FjordSim
config = fjord_config("examples/oslofjorden.jl")
prepare_bathymetry(config)
download_initial_state(config)
prepare_initial_state(config)
add_rivers(config)
download_boundary_data(config)
prepare_boundary_data(config)
download_atmosphere(config)
prepare_atmosphere(config)
run_simulation(config)
```

`build_simulation(config)` returns the assembled simulation without starting it, which is the entry
point for the REPL and the debugger. Snapshots and a per-launch log land in
`~/FjordSim_results/oslofjorden/`; `julia --project -m FjordSim --help` lists every subcommand.

What to re-run after changing the run window or the grid, how to archive a run before a fresh
rerun, and how to resume are in [`docs/running.md`](docs/running.md); running on a cloud GPU is in
[`gcp/README.md`](gcp/README.md).

## What a run is made of

A run reads prepared files only, each written by its own step onto — or along — the simulation grid.
Only the grid and the bathymetry are required. Every other component is a `FjordConfig` field
defaulting to `nothing`; a step whose config is `nothing` is a no-op rather than an error, and the
run does without that component. So a setup names what it wants and skips the rest.

| Component | `FjordConfig` field | Steps | Without it |
|---|---|---|---|
| Bathymetry | `bathymetry_config` | `prepare_bathymetry` | — (required) |
| Initial state | `initial_state_config` | `download_initial_state`, `prepare_initial_state` | no `FromInitialState` initial conditions |
| Rivers | `river_config` | `add_rivers` | no river inflow |
| Open boundaries | `boundary_data_config` | `download_boundary_data`, `prepare_boundary_data` | every lateral edge is a closed wall |
| Atmosphere | `atmosphere_config` | `download_atmosphere`, `prepare_atmosphere` | no prescribed atmosphere or radiation |
| Forcing | `forcing_config` | `download_forcing`, `prepare_forcing` | no interior forcing terms, no `FromForcing` initial conditions |
| Model | `model_config` | `run_simulation` | data preparation only; required with a run config |
| Biogeochemistry | `biogeochemistry_config` | `run_simulation` | no biogeochemical model |
| Boundary conditions | `boundary_condition_configs` | `run_simulation` | data preparation only; required with a run config, `()` for the model's defaults |
| Run | `run_config` | `run_simulation`, `validate_simulation` | data preparation only |

The model, its biogeochemistry, its boundary conditions and the run are separate components, side
by side with the data: `model_config` says what is built, `run_config` how it is run.

1. **Bathymetry** — the domain's depth and land/sea mask. `prepare_bathymetry` regrids a bathymetry
   source onto the simulation grid (built from the grid config's bounds and resolution) and writes the
   result as a NetCDF file. That file is what actually defines the model's `ImmersedBoundaryGrid`. The
   built-in source is the public [Geonorge](https://www.geonorge.no/) Sjøkart Dybdedata dataset.

2. **Initial state** — one record of a larger-scale ocean model at the run's `start_date`, for a
   realistic start without the whole-window download the forcing (below) needs.
   `download_initial_state` fetches a single day (built in: the NorKyst-800m daily mean of the day
   `start_date` falls in) and `prepare_initial_state` regrids it onto the simulation grid;
   `FromInitialState()` initial conditions read it. The file records the `start_date` it was prepared for, and a run starting at
   another one refuses it, so changing `start_date` means re-running both steps.

3. **Rivers** — `add_rivers` writes each river's daily discharge and temperature on the surface
   columns it enters through (built in: NVE's HydAPI), which the run adds as a volume source. It is its
   own config and file, independent of the forcing: it needs only the processed bathymetry and the
   run config's run window, `start_date` to `start_date + stop_time`, which its daily time
   axis spans.

4. **Open-boundary data** — the exterior state along the open lateral edges, which is what makes them
   open. `prepare_boundary_data` regrids the *hourly* NorKyst collection onto the boundary rows,
   including surface elevation and depth-averaged velocities, so the tide survives. The boundary
   condition itself is Flather on the transport plus Orlanski radiation with inflow nudging on
   everything else. Any subset of the four edges can be open — a fjord names one, a region in the
   open ocean names all four, and all of them share one prepared file.

5. **Atmosphere** — the surface forcing a run needs but a regional ocean model does not carry.
   `prepare_atmosphere` regrids a reanalysis product (built in: MET Norway's
   [NORA3](https://thredds.met.no/thredds/projects/nora3.html)) onto its own regular longitude/latitude
   grid, independent of the ocean grid since NumericalEarth interpolates between them, producing the
   eight variables a `PrescribedAtmosphere`/`PrescribedRadiation` pair needs.

6. **Forcing** — interior fields from a larger-scale ocean model. `prepare_forcing` regrids them
   (built in: NorKyst-800m daily means of temperature, salinity and velocity) onto the simulation
   grid, each with a rate field `λ` beside it. The run adds every variable to the model as a forcing
   term: where `|λ| < 1` it relaxes the model field towards the file's value at rate `λ` (s⁻¹) — the
   interior nudging ROMS does towards its climatology file — and where `|λ| > 1` it reads the value as
   a tracer flux into the cell, a source or sink. The built-in `prepare_forcing` writes `λ = 0`
   everywhere, so today these terms are zero and the file serves as the `FromForcing` initial state
   — the ocean starts from a realistic state rather than a uniform water column.

The initial state, rivers, open-boundary data and forcing go onto the bathymetry's grid, so their
steps need the processed bathymetry: that is where their land mask comes from. The atmosphere is
independent. `FromForcing` initial conditions need a prepared forcing file and `FromInitialState` a
prepared initial state; a setup with neither starts from constants, functions, or `FromResults` from
an earlier run.

### A new fjord

No setup is part of the package: every one is a standalone `.jl` file whose last expression is its
`FjordConfig`, passed by path to `--config`. A new fjord is a new file, written from scratch or built
on another with `base = include("oslofjorden.jl")`, as the variants in `examples/` are:

```bash
julia --project -m FjordSim run_simulation --config examples/oslofjorden_npzd.jl
```

## Some `run_simulation` details

`run_simulation` has no command-line options: everything is stated in the setup's
`RunConfig`. These fields control the run itself.

| Field | Meaning |
|---|---|
| `start_date` | the calendar instant model time zero stands for |
| `stop_time` | how long one run window is |
| `loops` | how many times to repeat that window, carrying the ocean state over |
| `initial_conditions` | where the first state comes from |
| `pickup` | resume from the newest checkpoint instead of starting fresh |


### Checkpointing

```julia
CheckpointWriter(interval = 30days, cleanup = true)
```

It writes the coupled model's full prognostic state as JLD2 — a few hundred MB each on a real grid.
The filename carries the iteration, so every fire writes a new file and `cleanup = true` is what
keeps only the newest; `false` accumulates one per interval. Only one checkpointing writer is
allowed.

### Continuing a run

Set `pickup = true` and the newest checkpoint is restored — clock, time
step, tracers, velocities, free surface, CATKE's diffusivities and the Adams-Bashforth tendencies — so
this is an *exact* continuation, and it supersedes `initial_conditions` (they are set at build time
and then overwritten). `pickup` without a `CheckpointWriter` is a configuration error, not a run that
fails later; `pickup` with no checkpoint on disk warns and starts from the beginning. The reverse
is refused too: a fresh run (`pickup = false`) into a `results_root` that already holds checkpoints,
since the two runs' checkpoints would mix.

The resumed run gets its own run tag, so its snapshots start at the resume point and the records
before it stay in the previous launch's file.

### Looping

`loops > 1` runs the same window repeatedly, saving the ocean state.
Each iteration resets all clocks to zero, saves the time step reached by the master, and writes its own file (`snapshots_ocean_<tag>_loop02.nc`) so that loops can be compared without overwriting each other.

### Starting from forcing, an initial state or a previous run

`initial_conditions` takes one of four shapes:

```julia
initial_conditions = (T = 5.0, S = 33.0)                       # constants, functions or fields
initial_conditions = FromForcing()                             # prepared forcing at `start_date`
initial_conditions = FromForcing(DateTime(2020, 3, 1))         # ...or at a named date
initial_conditions = FromInitialState()                        # the prepared initial state at `start_date`
initial_conditions = FromResults("snapshots_ocean_20260814T143012.nc")   # a previous run's last record
initial_conditions = FromResults("snapshots_ocean_20260814T143012.nc", DateTime(2020, 6, 1))
```

`FromInitialState()` reads the day the setup's `initial_state_config` downloaded and prepared:

```julia
initial_state_config = NorKystInitialStateConfig(
    data_root = data_root,
    output_directory = "norkyst_initial",
    parameters = ["temperature", "salinity", "u_eastward", "v_northward"],
),
```

All are a **lossy warm start**: the tracers from the config plus `u` and `v` are set, while the free surface, the turbulence state and the timestepper tendencies keep their defaults, so the first hours are a barotropic adjustment.
That is the difference from `pickup`, which is exact — the two are complementary, not alternatives.
