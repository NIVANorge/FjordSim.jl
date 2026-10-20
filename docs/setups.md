# Adding a simulation setup

> Read this before writing a new setup or changing an existing one's parameters. The worked setups
> are in `examples/`; the rationale behind each of their values is in their own comments.

## What a setup is

A setup is a Julia script whose last expression is a `FjordConfig`. It is not part of the package:
it is passed by path to every step, and `fjord_config` evaluates it afresh each time.

```julia
using FjordSim

data_root = fjord_data_root("hardangerfjorden")

FjordConfig(
    grid_config = EvenGrid(...),
    bathymetry_config = DybdedataConfig(data_root = data_root, ...),
    # every other field is optional
)
```

```bash
julia --project -m FjordSim prepare_bathymetry --config path/to/hardangerfjorden.jl
```

Two ways to write one:

- **From scratch**, or by copying the example closest to the new domain: `examples/oslofjorden.jl`
  (inner fjord, open southern edge), `examples/drammensfjorden.jl` (a sub-domain sharing Oslofjord's
  downloads), `examples/lista.jl` (open coast, three open edges).
- **As a variant of another setup**: `base = include("oslofjorden.jl")` resolves beside the
  including file and returns that setup's config, whose fields the variant then changes or wraps
  (`examples/oslofjorden_validation.jl`, `examples/oslofjorden_npzd.jl`). The included script's
  top-level names (`data_root`, `FT`) become the variant's globals too, so assign the variant's own
  after the `include`. `gcp/fjordsim-gcp` uploads every `.jl` file beside the config, so the
  include works on the VM as well.

Name the file after the fjord, lowercase, the same name `fjord_data_root` gets: `gcp/fjordsim-gcp`
takes the basename as the fjord whose data and results directories it stages.

## What it must and may contain

Only the grid and the bathymetry are required. Every other field defaults to `nothing`, and every
step whose config is `nothing` is a no-op:

| Field | Built-in config | Steps | Prepared file | Needed for |
|---|---|---|---|---|
| `grid_config` | `EvenGrid` | — | — | always |
| `bathymetry_config` | `DybdedataConfig` (Geonorge) | `prepare_bathymetry` | `bathymetry.nc` | always: it *is* the model grid |
| `forcing_config` | `NorKystConfig`, `NorKystHindcastConfig` | `download_forcing`, `prepare_forcing` | `forcing.nc` | `FromForcing` initial conditions |
| `river_config` | `NVERiversConfig` | `add_rivers` | `rivers_nve.nc` | river inflow; `RiverMouthMixing` |
| `boundary_config` | `NorKystBoundariesConfig`, `NorKystHindcastBoundariesConfig` | `download_boundaries`, `prepare_boundaries` | `boundaries.nc` | any open edge; `OpenLateralBoundaryFromData` |
| `atmosphere_config` | `NORA3Config` | `download_atmosphere`, `prepare_atmosphere` | `atmosphere.nc` | the prescribed atmosphere and radiation the surface fluxes are computed from |
| `simulation_config` | `SimulationConfig` | `run_simulation`, `validate_simulation` | results | running at all |

NorKyst's operational archive (`NorKystConfig`, `NorKystBoundariesConfig`) begins on 2017-02-20; a
window before that needs the hindcast configs, which cover 2012 onwards. Each field is a type
parameter of `FjordConfig`, so it is chosen when the config is built and cannot be set to `nothing`
on an existing instance afterwards.

A new data source is a new config type, not an edit to the pipeline: `docs/adding-a-source.md`.

## Paths

- **Roots.** Build the data root with `fjord_data_root("<fjord>")` and the simulation config's
  `results_root` with `fjord_results_root("<fjord>")`, never from `homedir()`. Both read
  `FJORDSIM_DATA_ROOT` / `FJORDSIM_RESULTS_ROOT`, which is how one setup runs unchanged on a laptop
  and on the cloud VM. Results are never under the data root.
- **File names** (`output_file`, `plot_file`, `output_directory`, `geodatabase_file`) are relative
  to the config's own `data_root`; an absolute path relocates just that one.
- **Sharing a download** with another fjord is an absolute path into that fjord's root, built from
  `fjord_data_root("<other>")` so a relocation moves both — `examples/drammensfjorden.jl` shares
  Oslofjord's FileGDB and NORA3 months this way. Share only from a domain that *contains* the new one, for the same years: a month missing from the shared directory is
  downloaded subset to the smaller box and then shrinks the larger setup's own prepare step.
- **Never shared:** prepared files (their names carry no fjord, so two setups would overwrite each
  other), the NVE caches (they record the domain they were fetched for under a fixed name), and the
  open-boundary band (it is cut along that setup's own edges).
- One results directory holds one run (`docs/running.md`, "One run per results directory"); a
  variant that is a different experiment gets its own `fjord_results_root` name.

## The grid

`EvenGrid` takes `size`, `halo`, `longitude`, `latitude`, `z_faces` and `float_type` (default
`Float32`).

- **Float type.** Define `FT = Float32` (or `Float64`) once, give it to `EvenGrid(float_type = FT)`
  and to every model component that takes one — `WENO(FT)`, `WENOVectorInvariant(FT)`,
  `CATKEVerticalDiffusivity(FT; ...)`, `SeawaterBuoyancy(FT, ...)`, `TEOS10EquationOfState(FT)`,
  `HydrostaticSphericalCoriolis(FT)`, `FreezingLimitedOceanTemperature(FT)`. Never set
  `Oceananigans.defaults.FloatType`: it is process-wide.
- **Box.** A closed edge should be land along its whole length; an edge cutting open water must be
  opened (`boundary_config.open_edges`) or the domain moved. Check the plotted bathymetry after the
  first `prepare_bathymetry`.
- **`z_faces`** run bottom to top, end at `0.0`, and number `size[3] + 1`.
  - Fit the deepest face to the deepest cell `prepare_bathymetry` writes (read it from
    `bathymetry.nc`), not to a round number: deeper wastes levels and shortens the barotropic
    substep (`sqrt(g·Lz)`), shallower silently truncates the basin.
  - Avoid thick bottom layers. A geometric stretch from a ~1 m surface cell, with no two adjacent
    layers differing by more than a factor ~1.33, keeps a column's floor in a thin layer; a 25 m
    bottom cell holds a site's whole near-bottom water in one cell, and that is where tracer
    extremes collect.
  - More levels give the seabed more faces to cross, so more laterally isolated bottom cells.
- **Changing `size`, `z_faces` or the box is a data change**: `simulation_grid` reads
  `bathymetry.nc`, not the config. Re-run the prepare steps `docs/running.md` lists.

## The bathymetry

`DybdedataConfig` requires `data_root`, `output_file` and `plot_file`. Its defaults leave every
optional smoothing stage off:

| Field | Default | What it does |
|---|---|---|
| `geodatabase_file` | `GEONORGE_DYBDEDATA_GDB` | the national FileGDB, downloaded (4.5 GB) into `data_root` if missing |
| `raw_resolution_factor` | `4` | native gridding resolution, as a multiple of the model resolution |
| `padding_cells` | `0` | model-cell widths of soundings read around the box |
| `include_contours`, `contour_stride` | `true`, `1` | use depth contours besides soundings |
| `max_sounding_distance` | `Inf` | turn wet cells farther than this many metres from any sounding into land |
| `interpolation_passes` | `1` | passed to `NumericalEarth.regrid_bathymetry` |
| `major_basins` | `1` | keep this many largest connected water bodies (drops lakes) |
| `minimum_depth` | `0.0` | deepen every wet cell shallower than this, m |
| `open_boundary_land_cells` | `0` | flood land within this many rows of an open edge |
| `max_island_cells` | `0` | remove interior land components up to this size |
| `close_narrow_passages` | `false` | turn one-cell-wide passages into land where both sides stay connected |
| `spike_ratio` | `0.0` | replace a cell shallower than this fraction of its neighbours' median depth by that median |
| `minimum_cell_fraction` | `0.0` | snap partial bottom cells; set it to `0.2`, `PartialCellBottom`'s `minimum_fractional_cell_height` |
| `max_slope_factor` | `0.0` | limit the bottom slope between neighbours |
| `preserve_sills` | `false` | restore sills the regrid averaged away |
| `geonorge_cache`, `regrid_cache` | `true`, `true` | reuse intermediate files |

`DybdedataConfig`'s docstring states the failure each stage prevents and `docs/architecture.md`
(Bathymetry) why their order is fixed. Check that `max_island_cells` stays below the smallest real
island, and that a known sill comes out at its real depth.

## The simulation config

`SimulationConfig` and every config it nests have **no defaults**: each field is a scientific choice
about one fjord, stated in the setup, so a field added later breaks every setup until each names it.
`README.md` lists the fields; the examples state all of them. Coefficients that depend on the grid
spacing have to be worked out for the new grid rather than copied:

- **Explicit diffusion stability.** `AdaptiveTimeStep` measures only the advective CFL. A biharmonic
  `ν₄` needs `Δt ≤ Δx⁴ / (32 ν₄)` and damps the 2Δx mode with e-folding time `Δx⁴ / (16 ν₄)`; a
  harmonic `ν` (including `BoundarySponge`'s `viscosity` and `diffusivity`) needs
  `Δt ≤ Δx² / (4 ν)`. Keep both limits well above `max_time_step`, or the diffusion, not the CFL,
  sets the time step. Taking a value from a setup with a different cell size means scaling it by
  `(Δx_new / Δx_old)⁴` for biharmonic and `(Δx_new / Δx_old)²` for harmonic coefficients.
- **`tracer_advection`** is one scheme for every tracer, not a `NamedTuple`: tracers a `NamedTuple`
  omits, such as CATKE's `e`, fall back to an unbounded `Centered()` scheme.
- **`CATKEVerticalDiffusivity`'s `minimum_tke`** sets the background `κ = 0.098·e_min/N` wherever the
  prognostic TKE is below it, which in a stratified fjord is most of the column, so it is effectively
  the vertical-mixing knob.
- **River inflow.** `NVERiversConfig` requires `data_root`, `outlets`, `years` and
  `maximum_inflow_velocity`. `maximum_inflow_velocity` bounds the vertical Courant number
  `w·Δt/Δz` through the top cell: at 1e-3 m/s, a 1 m top cell and a 180 s `max_time_step` it is
  0.18. `minimum_discharge` defaults to `Inf`, meaning no discovery: `outlets` must then state every
  river and its mouth. A discharge in m³/s instead discovers every mouth at or above it from NVE's
  river network, and `outlets` become overrides that attach gauges.
- **Open boundaries.** `OpenLateralBoundaryFromData`'s `inflow_timescale` relaxes by `Δt/τ` per step;
  compare that with the advective rate into the boundary cell, `u/Δx`, or the nudging does nothing.
  Too strong an inflow nudge reflects. A `BoundarySponge` damps what the open edges radiate.
- **Initial conditions.** `FromForcing()` needs a `forcing_config`; `FromResults(path)` continues
  from an earlier run's snapshot; a `NamedTuple` of constants or functions needs nothing.
- **Run window.** `start_date` and `stop_time` decide what the prepare steps write, since they pad
  every time axis to the window — changing either is a data change too.

## Bringing a new setup up

1. Load it: `fjord_config("path/to/it.jl")` must return a `FjordConfig`.
2. `prepare_bathymetry`, then inspect `bathymetry.png` and `bathymetry.nc`: closed edges are land,
   the deepest cell matches `z_faces`, sills and channels survive.
3. The other prepare steps the setup names, in the order `docs/running.md` lists.
4. A short run first (days, not the full window), on the GPU it will run on, before the long one.
