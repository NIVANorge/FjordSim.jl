# Known issues

## `examples/oslofjorden/oslofjorden.jl` (160 m, 2024-10 to 2026-10)

- `max_island_cells = 6` and the check of the deepest cell against `z_faces`' −467 m (the box's
  deepest sounding is 470 m) were measured on 193 m cells, not on the 160 m grid.
- The eight gauge overrides were chosen for 2020. Whether each has daily values on HydAPI over
  2024–2026 is unchecked.
- `validate_simulation` scores only Kartverket water level, which it fetches itself. The NIVA casts
  of June 2024 – December 2025 are exported to `observations/niva_ctd_2024-2025`, but it reads them
  only from a `CsvObservations` passed to it. AquaMonitor holds no 2026 CTD casts yet, only bottle
  chemistry (checked 2026-11-08).
- `CastCorrectedInitialState` is untested in a run. It leaves 121 surface cell pairs between cast
  and Norkyst columns, 6.3 g/kg apart (median), at the Drøbak and Svelvik sills and around inner-fjord
  pockets it does not reach; the cast profiles add 3200 small density inversions (≤ 0.43 kg/m³,
  against 1.31 without them).
- The run of stamp 20261107T210205 started before `CastCorrectedInitialState`: its sill-basin water
  starts 1–3 g/kg too salty and about 2 °C too cold, and is renewed only slowly (Ep1 at 120 m:
  34.1 to 33.9 g/kg from October to January, the casts 32.9).
- The initial state has T = S = 0 in 20 wet cells: the south-west corner of the open edge, i 1–2,
  j 1, the top ten levels. Those nodes lie 535–617 m west of the downloaded Norkyst subset (see
  "Downloads" below), where the interpolation returns 0; no level of their columns is covered, so
  `TargetFill` leaves them.
- In November 2024 – January 2025 the top 5 m are 3.9 g/kg too salty against the NIVA casts, and
  the 25 m − 1 m salinity difference 2.9 g/kg too weak. Norkyst has the same error (+3.8 g/kg); the
  cause, river input, mixing or the boundary, is unknown. In the Inner Oslofjord the top 5 m are
  also 1.2 °C too warm, where Norkyst is within 0.1 °C.
- The non-tidal water level has 60 % of the observed variance at all three gauges, as in Norkyst,
  whose surge the open boundary carries in; the set-down of 9 December 2024 reached −0.25 m against
  −0.6 m observed. M4 is half the observed (0.6–1.6 cm against 1.2–2.5) and 31–46° early; the other
  constituents are within 0.5 cm.
- Single cells reach 29.6 °C (23 November 2024), with at most two above 25 °C on any day of
  October 2024 – January 2025; which cells, and why, is unexamined.
- The `:hydrography` station writer's list is the 2014–15 programme's: OF-1, OF-5 and SF-1 have no
  casts in 2024–2025, while the programme's SF-3, VAL-1, BO-1, I-1, I-4, KF-1 and SAN-5 have no
  hourly output and are scored from the daily snapshots.
- Holmenfjorden (Cj1) is about 25 m deep within a cell of the station against casts to 56 m.
- Below 104 m the Drøbak Sound south of the sill (Im2) is a pocket closed off from the open edge, so
  `TargetFill` fills it from above: 9.9 °C at 150 m on 2024-10-01 against 7.1 in the cast of
  2024-09-26 (7.8 with the earlier same-depth fill). `TargetFill` also doubles the cells denser than
  the one below them, 4746 to 11 278, the largest 0.72 to 1.31 kg/m³.
- `oslofjorden_npzd.jl` runs 2024–2026 with nutrient profiles measured in 2020.
- No sea ice. `SlabSeaIce` (no drift) grew runaway ice on this setup, so it runs
  `FreezingLimitedOceanTemperature`, which clamps every level at freezing and discards the heat: no
  ice insulates the surface and no melt freshens it. In January 2025 the surface is 0.3–2.6 g/kg
  saltier than the slab-ice run where that run sat at freezing, and 1729 surface cells sit at
  freezing on 6 Jan against its 473. The slab-ice runs are archived as `oslofjorden_2026-11-06`,
  `oslofjorden_2026-11-07` and `oslofjorden_2026-11-07-2`:
  - The cause is WENO's fallback at a boundary. Oceananigans' default `WENO()` drops to order 3
    and then to `Centered(order=2)` on faces too close to land, the seabed or the surface to fit
    its stencil. A centred face exports the mean of the two cells, so flow out of a cell colder than
    its neighbour removes heat the cell does not have, and the frazil reset turns the deficit into
    ice. At Larvik this happened on horizontal faces beside land: 305 m of ice by April 2025 at
    (17, 61) in the Numedalslågen channel, fixed by upwinding coastal columns (`narrow_cells`),
    after which Larvik peaked at 0.19 m.
  - At Svelvik it happens on the face under the top cell, which `RiverMouthAdvection` leaves to
    WENO. Ice in (147, 458) and (148, 458), 4 m and 2 m deep beside the channel, grew from 16 to
    93 cm on 10–16 January 2025. A 1 m fresh top layer (S 2–3) sits at freezing over water 6–7 K
    warmer (S 22–28), and the converging flow sinks out of it at 2–6·10⁻⁴ m/s. Replayed on the
    daily snapshots of 11–13 January, every downwelling face under a top cell carried exactly the
    two-cell mean, and every top cell at freezing in i 146–148, j 452–461 was cooled by 7–121 K/day
    while horizontal advection mostly warmed it. With `WENO(FT; minimum_buffer_upwind_order = 1)`,
    which falls back to first-order upwind instead, 14 of those 16 are no longer cooled. That is
    untested in a run. Under the clamp the same sink is discarded instead of growing ice. The ice
    also outlived the warm water because NumericalEarth's default `ThreeEquationHeatFlux` uses a
    fixed friction velocity of 0.002 m/s.
- The same fallback takes two bottom cells to −1.90 °C in October: (346, 39, 24) at 64 m, a corner,
  and (344, 35, 27) at 49 m, a dead end draining through its top face. The coastal rule misses
  them: their surface columns have four wet neighbours, and it upwinds only horizontal faces.

## `examples/mjosa/mjosa.jl`

- Gudbrandsdalslågen, three quarters of the inflow, has no 2025 water temperature after August
  (Hovdefossen `2.656.0` ends 2025-08-19), so it carries Glomma's at Glomstadfoss. Glomma reads
  down to −1.3 °C in winter; `add_rivers` raises every river temperature below 0 °C to 0 °C.
- The setup runs `SlabSeaIce` (no drift) under NORA3 snowfall, which froze the Mjøsa grid on CPU
  from a near-freezing surface layer under 2025-01-01 NORA3. The ice is not scored against NVE's
  freeze-up and break-up dates.
- The run is not scored against NIVA's Mjøsa monitoring profiles (Skreia, Brøttum, Furnes) in
  Vannmiljø. Unchecked: the Vorma outlet's mean against the ~344 m³/s at Ertesekken, η against the
  120.3–123.0 m observed at Hamar in 2025, and whether T stays smooth in the outlet's columns.

## Downloads

- The Norkyst initial-state and forcing downloads (`NorKystV3InitialStateConfig`,
  `NorKystInitialStateConfig`, `NorKystConfig`, `NorKystHindcastConfig`) take the bounding box of the
  source cells whose centres lie inside the target's longitude–latitude box, with no margin. On the
  rotated Norkyst grid the box's corners can lie outside it, and a wet target cell there
  interpolates beyond the subset to 0 (Oslofjord, above). The boundary-data configs have a `margin`
  for this.

## Design

`docs/dispatch-map.md` lists every supertype's hooks and which built-in types implement them. These
are the places where a new component still cannot plug in through dispatch alone.

- The longitude–latitude grid assumption is in no type. A grid config only promises `domain_grid`
  and `simulation_grid`, but `regrid_forcing`, `regrid_initial_state`, `regrid_boundary_data`,
  `place_rivers` and the prepared files' layout (`src/Regridding.jl`) read the target grid's
  coordinates with `λnodes`/`φnodes`. A rectilinear or single-column grid config therefore loads
  and then fails inside the regridding. The fix is either a method of those work functions per grid
  kind, or a grid-config hook returning each point's longitude and latitude.
- Validation's station-writer contract is not hooks. `station_output` opts a writer in, but
  `score_run` (`model_series` in `src/Validation/Validation.jl`) then reads the writer's
  `variables`, `stations` and `output_file` fields and expects one file per station, named by
  `station_tag`. A writer that stores its stations another way, such as `SectionWriter` with many
  profiles in one file, is not scored.
