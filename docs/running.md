# Running a fjord

How to take a setup from nothing to a scored run, and the rules that keep its data and results
consistent. The quick start and the `run_simulation` fields are in the `README.md`; running on a
cloud GPU is in `gcp/README.md`.

`--config` is the only option of every subcommand, and it is required. It takes the path to a setup:
a `.jl` file whose last expression is a `FjordConfig`, such as `examples/oslofjorden/oslofjorden.jl`. No setup
is part of the package. Everything else is a field of that config.

## The steps

In order. Each step needs the ones above it that its setup configures; a step the setup does not
configure is a no-op.

| Step | Needs | Writes (under the fjord's data directory) | Notes |
|---|---|---|---|
| `prepare_bathymetry` | — | `bathymetry.nc` | downloads the Geonorge FileGDB (4.5 GB) on first use |
| `download_forcing` | — | raw NorKyst months | network |
| `prepare_forcing` | bathymetry, download | `forcing.nc` | pure regrid |
| `add_rivers` | bathymetry | `rivers_nve.nc` | needs `NVE_API_KEY` (free: https://hydapi.nve.no/Users) |
| `download_boundary_data` | — | hourly NorKyst band along the open edges | network |
| `prepare_boundary_data` | bathymetry, download | `boundaries.nc` | pure regrid |
| `download_atmosphere` | — | raw NORA3 months | slow: ~10000 OPeNDAP reads per simulated year |
| `prepare_atmosphere` | download | `atmosphere.nc` | pure regrid, independent of the ocean grid |
| `run_simulation` | every prepared file the setup names | results under `results_root` | needs a GPU in practice |
| `validate_simulation` | a finished run | `<results_root>/validation/` | only Kartverket water level is fetched automatically |

```bash
julia --project -m FjordSim prepare_bathymetry --config examples/oslofjorden/oslofjorden.jl
julia --project -m FjordSim run_simulation --config examples/oslofjorden/oslofjorden_npzd.jl
julia --project -m FjordSim validate_simulation --config examples/oslofjorden/oslofjorden.jl
```

The same without `-m`, Julia 1.12's package entry point:

```bash
julia --project -e 'using FjordSim; FjordSim.main(ARGS)' -- prepare_bathymetry --config examples/oslofjorden/oslofjorden.jl
```

`scripts/surface_videos.jl` turns one file of a `SurfaceSnapshotWriter` (and, with `--sea-ice`, the
matching `SeaIceSnapshotWriter` file) into one MP4 per variable; its header has the usage.
`scripts/section_videos.jl` makes one MP4 of a `SnapshotWriter`'s records (all its parts, e.g. a
year of daily snapshots) along a section that follows the deepest water between given points, with
the velocity along the section and, from the surface and sea-ice snapshot files if given, strips of
the surface elevation and the ice and snow.

### From the REPL, and debugging a step

Each subcommand is a method on the generic function of the same name, taking the whole setup:

```julia
using Pkg; Pkg.activate(".")
using FjordSim
config = fjord_config("examples/oslofjorden/oslofjorden.jl")
prepare_bathymetry(config)

using Debugger
@enter add_rivers(config)
```

`build_simulation(config)` returns the assembled `Simulation` without starting it — the entry point
for inspecting a model or stepping into it.

## What to re-run after a config change

The prepared files are regridded onto the grid and padded to their dataset's period *at prepare
time*, so some config edits silently do nothing until the files are rewritten.

**The run's `start_date` or `stop_date`.** Nothing, as long as every prepared file still spans the
new window: each dataset is downloaded and prepared for its own config's `start_date` to `stop_date`,
independent of the run's, and `build_simulation` checks the run window lies inside each one. The
exception is the initial state, whose config's `start_date` has to be the run's: `FromInitialState`
refuses a file prepared for another one.

**A dataset config's `start_date` or `stop_date`.** That dataset's download fetches every month the
period touches, its end exclusive (`coverage_months`), skipping the months already present, and its
prepare step pads to the period's ends; the river file's daily axis *is* the river config's period.
Re-run that dataset's download and prepare steps — for the rivers, `add_rivers` does both:

```bash
julia --project -m FjordSim download_atmosphere --config examples/oslofjorden/oslofjorden.jl
julia --project -m FjordSim prepare_atmosphere  --config examples/oslofjorden/oslofjorden.jl
```

For the initial state, `download_initial_state` and `prepare_initial_state`.

**`z_faces`, the grid `size`, or any `bathymetry_config` smoothing knob.** `z_faces` is written into
`bathymetry.nc` and `simulation_grid` reads the *file*, not the config; the land mask the forcing and
boundary prepares build their masks from comes from the same file. Re-run, in this order, with
`prepare_forcing` and `prepare_initial_state` after `prepare_bathymetry` for a setup that names
them:

```bash
julia --project -m FjordSim prepare_bathymetry  --config examples/oslofjorden/oslofjorden.jl
julia --project -m FjordSim add_rivers          --config examples/oslofjorden/oslofjorden.jl
julia --project -m FjordSim prepare_boundary_data  --config examples/oslofjorden/oslofjorden.jl
```

`prepare_atmosphere` is not in this list: the atmosphere is 2D on its own lon/lat grid. Existing
checkpoints are invalid after either change.

No download step is affected by a grid or bathymetry change, and `loops` re-runs nothing.

## One run per results directory

Checkpoints carry no run tag, so a results directory holds exactly one resumable run: `pickup`
resumes the highest-numbered checkpoint it finds, whichever launch wrote it. A fresh run into a used
directory would mix with the old one, and analyses that read every stamp would splice them together.

- **Every run uses the setup's own config and `results_root`.** A rerun, a test or a changed
  parameter is an edit to the setup, not a new config file.
- **A finished run is archived**, so the undated name is free for the next one. Rename `<name>` to
  `<name>_<YYYY-MM-DD>`, the UTC date of that run's last output (`-2` if taken) — locally and, if
  the run went through GCP, in the bucket. A `--resume` is not a fresh run and archives nothing.

Locally:

```bash
mv ~/FjordSim_results/oslofjorden_validation ~/FjordSim_results/oslofjorden_validation_2026-10-16
```

In the bucket: `gcp/README.md`, "When a run finishes".

The rule is enforced twice: `run_simulation` refuses a fresh run (`pickup = false`) into a
`results_root` that already holds checkpoints, and `gcp/fjordsim-gcp run` refuses a fresh
`run_simulation` while the bucket's `results/<name>/` does.

## Resuming

Set `pickup = true` in the config (or pass `--resume` to `gcp/fjordsim-gcp run`, which sets it) and
the newest checkpoint is restored exactly. It needs a `CheckpointWriter` among the writers. See
"Continuing a run" in the `README.md`.
