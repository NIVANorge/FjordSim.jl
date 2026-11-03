# TODO

## `examples/oslofjorden/oslofjorden.jl` at 160 m, 2024-10 to 2026-10

- [ ] Bring it up: every prepare step (nothing is prepared at 160 m), then check the deepest cell
      against `z_faces`' −467 m (the box's deepest sounding is 470 m) and whether
      `max_island_cells = 6`, chosen on 193 m cells, still deletes only one-cell ridges.
- [ ] Check the eight gauge overrides have daily values on HydAPI over 2024–2026; they were chosen for
      2020.
- [ ] Observations for 2024–2026: NIVA CTD exports for `validate_simulation` (the cached ones are
      2014-15); Kartverket water level it fetches itself.
- [ ] `oslofjorden_npzd.jl` now runs 2024–2026 with nutrient profiles measured in 2020.

## `examples/mjosa/mjosa.jl`

- [ ] Gudbrandsdalslågen, three quarters of the inflow, has no 2025 water temperature after August
      (Hovdefossen `2.656.0` ends 2025-08-19), so it carries Glomma's at Glomstadfoss. Replace the
      proxy if NVE fills the series. Glomma reads down to −1.3 °C in winter; `add_rivers` raises every
      river temperature below 0 °C to 0 °C.
- [ ] Ice: the setup runs `SlabSeaIce` (no drift) under NORA3 snowfall, which froze the Mjøsa grid
      on CPU from a near-freezing surface layer under 2025-01-01 NORA3. Still to do:
  - [ ] Score the ice against NVE's freeze-up and break-up dates.
- [ ] Score the run against NIVA's Mjøsa monitoring profiles (Skreia, Brøttum, Furnes) in Vannmiljø.
      Check the Vorma outlet's mean against the ~344 m³/s at Ertesekken, η against the
      120.3–123.0 m observed at Hamar in 2025, and that T stays smooth in the outlet's columns.

## Design

`docs/dispatch-map.md` lists every supertype's hooks and which built-in types implement them. These
are the places where a new component still cannot plug in through dispatch alone.

- [ ] State the longitude–latitude grid assumption in a type. A grid config only promises
      `domain_grid` and `simulation_grid`, but `regrid_forcing`, `regrid_initial_state`,
      `regrid_boundary_data`, `place_rivers` and the prepared files' layout (`src/Regridding.jl`)
      read the target grid's coordinates with `λnodes`/`φnodes`. A rectilinear or single-column grid
      config therefore loads and then fails inside the regridding. Either give those work functions
      a method per grid kind, or add a grid-config hook returning each point's longitude and
      latitude.
- [ ] Turn validation's station-writer contract into hooks. `station_output` opts a writer in, but
      `score_run` (`model_series` in `src/Validation/Validation.jl`) then reads the writer's
      `variables`, `stations` and `output_file` fields and expects one file per station, named by
      `station_tag`. A writer that stores its stations another way, such as `SectionWriter` with many
      profiles in one file, is not scored.
