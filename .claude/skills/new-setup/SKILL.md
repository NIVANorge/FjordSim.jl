---
name: new-setup
description: Checklist for writing a new FjordSim simulation setup (a new fjord, domain or experiment config file) or a variant of an existing one — choosing the domain, grid and z_faces, data sources, data roots, and the model coefficients that depend on the grid. Use when asked to set up, create or configure a new fjord simulation, a new domain, or a new experiment based on an existing setup.
---

# Writing a new FjordSim setup

The rules and defaults are in `docs/setups.md`; read it first. The worked setups are in `examples/`.
This is the order to work through.

1. **Pin down what the run is for.** Ask the user for the domain (box), the run window, which
   components it needs (forcing, initial state, rivers, open boundaries, atmosphere, a run at all)
   and what it will be judged on. A variant of an existing setup is a new file that does
   `base = include("<setup>.jl")`, not an edit to the original; ask before changing an existing setup.
2. **Check the data covers it.** The run window against the source archives (NorKyst's operational
   archive starts 2017-02-20; earlier needs the hindcast configs), and the box against Geonorge's
   coverage. A shared download only from a domain that contains the new one.
3. **Write the file** from the closest example: `fjord_data_root`/`fjord_results_root` with the
   file's own basename, one `FT` passed everywhere, every `CoupledHydrostaticModel` and
   `SimulationConfig` field stated, and the `boundary_condition_configs` tuple. Scale every
   grid-spacing-dependent coefficient to the new Δx with the formulas in `docs/setups.md` rather than
   copying it, and state in a comment the number that justifies each value.
4. **Load it**: `julia --project -e 'using FjordSim; fjord_config("<file>.jl")'` must return a
   `FjordConfig`.
5. **Bathymetry first.** Confirm with the user before running `prepare_bathymetry` (it may download
   the 4.5 GB FileGDB). Then check `bathymetry.png`/`bathymetry.nc`: closed edges are land, the
   deepest cell matches the deepest face of `z_faces` (refit and re-run if not), sills and channels
   survive.
6. **Hand over to the `run-simulation` skill** for the remaining prepare steps and a short test run
   before the full window.
