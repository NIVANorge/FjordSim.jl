# Known issues

## `examples/oslofjorden/oslofjorden.jl` (160 m, 2024-10 to 2026-10)

- `max_island_cells = 6` and the check of the deepest cell against `z_faces`' −467 m (the box's
  deepest sounding is 470 m) were measured on 193 m cells, not on the 160 m grid.
- The eight gauge overrides were chosen for 2020. Whether each has daily values on HydAPI over
  2024–2026 is unchecked.
- `validate_simulation` scores only Kartverket water level, which it fetches itself. The NIVA casts
  are exported to `observations/niva_ctd_2024-2025` (AquaMonitor, June 2024 – December 2025) and
  `observations/niva_ctd_odm2_2024-2026` (NIVA's ODM2 database, October 2024 – September 2026), IMR's
  ØkoKyst and Tisler casts to `observations/imr_ctd_2024-2026` (Vannmiljø, to March 2026), but it
  reads them only from a `CsvObservations` passed to it.
- `CastCorrectedInitialState` is untested in a run. It leaves 121 surface cell pairs between cast
  and Norkyst columns, 6.3 g/kg apart (median), at the Drøbak and Svelvik sills and around inner-fjord
  pockets it does not reach; the cast profiles add 3200 small density inversions (≤ 0.43 kg/m³,
  against 1.31 without them).
- The run of stamp 20261107T210205 started before `CastCorrectedInitialState`: its sill-basin water
  starts 1–3 g/kg too salty and about 2 °C too cold, and is renewed only slowly (Ep1 at 120 m:
  34.1 to 33.9 g/kg from October to January, the casts 32.9). Over June 2025 – May 2026 the Inner
  Oslofjord below 50 m is −0.28 g/kg, but Ep1 at 120 m is still 1.4–2.4 °C too cold and D-2 at 100 m
  0.9 g/kg too salty in March 2026 (3.2 at the start). Ep1 at 120 m is never renewed in the run: by
  August 2026 its start-state salt has drained to 32.5 g/kg (casts 32.9) at 6.5 °C (casts 8.9), and
  D-2 at 100 m is still 0.5 g/kg too salty in September 2026.
- In November 2024 – January 2025 the top 5 m are 3.9 g/kg too salty against the NIVA casts, and
  the 25 m − 1 m salinity difference 2.9 g/kg too weak. Norkyst has the same error (+3.8 g/kg); the
  cause, river input, mixing or the boundary, is unknown. IMR's NorFjords 160 m hindcast (ROMS, GLS
  mixing) is saltier still: on the casts all three reach, the November top 5 m are +4.5 g/kg in
  FjordSim, +5.2 in Norkyst and +6.7 in NorFjords, so FjordSim's mixing is the least likely cause.
  In the Inner Oslofjord the top 5 m are also 1.2 °C too warm, where Norkyst is within 0.1 °C.
- From February to April 2025 the water at 20–50 m is 1.0–1.5 °C too cold. Norkyst and NorFjords
  are 0.7–1.3 °C too cold there too, NorFjords under AROME and GLS mixing, so it does not point at
  NORA3 or CATKE.
- The deep-water renewal of February 2026 is missed, and from March to May 2026 the water below
  50 m is 2.0–2.5 °C too cold. The casts show warm, salty water entering the inner basins (Fl1 at
  100 m 32.8 to 33.5 g/kg, Ep1 at 120 m to 9.1 °C); FjordSim's basins keep freshening and cool to
  5.6–6.7 °C. Outside the Drøbak sill (Im2, 100 m) the casts read 35.0 g/kg, FjordSim 34.1–34.3,
  Norkyst 34.2–34.4 and NorFjords 34.8–35.0. At the fjord mouth in March (IMR's Torbjørnskjær and
  Færder, 30–75 m) the casts read 7.9–8.3 °C and 34.5–34.9 g/kg, FjordSim 4.7–5.0 °C and
  32.6–33.6 g/kg, Norkyst 4.4–5.3 °C; in February the casts there held a cold brackish layer at 30 m
  the models lack. Norkyst and NorFjords are 1.8–2.2 and 1.5–1.9 °C too cold below 50 m too. Norkyst
  is the open boundary's source, so its Skagerrak water at 30–100 m is the first suspect. The model
  renews Vestfjorden and Steilene instead in early August 2026 (Fl1 at 100 m 5.7 to 7.9 °C), after
  an inflow in July 2026 that raises the volume-mean salinity 1.1 g/kg in three weeks; the casts
  show no inflow then.
- Glomma's fresh water stays too close to its mouths. In the top 1.5 m over November 2024 –
  September 2025, I-4 off Fredrikstad is 3.4 g/kg too fresh (casts 3.4 g/kg), while I-1, 8 km
  farther south in the Hvaler skerries, is 7.9 g/kg too salty (casts 12.2) and S-9 4.8 (casts 19.6).
  Norkyst is 3.1 and 2.8 g/kg too salty at I-1 and S-9, NorFjords 5.1 and 3.1. At 5 m all three are
  about 7 g/kg too salty at I-1 outside February–April, so there the error is shared. Unchecked: the
  2/3–1/3 split between Østerelva and Vesterelva, and how open the skerries' channels are on the
  160 m grid.
- The non-tidal water level is too weak at all three gauges: as daily means over November 2024 –
  September 2025 it has 66–67 % of the observed standard deviation, Norkyst 61–62 % and NorFjords,
  with another open boundary, 87–88 %. The missing surge is in Norkyst's sea level, which the open
  boundary carries in; the set-down of 9 December 2024 reached −0.25 m against −0.6 m observed. M4
  is half the observed (0.6–1.6 cm against 1.2–2.5) and 31–46° early; the other constituents are
  within 0.5 cm. Over June 2025 – May 2026 the daily means carry 79 % (Norkyst 71–73 %, NorFjords
  92 %), over the whole run 74 % (Norkyst 67–69 %).
- The summer heat stays in the top metre. In the heat wave of 15–29 July 2025 up to 48 088 top
  cells exceed 25 °C (21 July), the warmest 29–33.3 °C in columns two to six cells (to 9 m) deep,
  several of them river water (1–6 g/kg). The casts of 21–22 July read 24 °C uniform to 2.5 m;
  FjordSim is 0.9–1.3 °C warmer at 0.5 m and 3.1–3.5 °C colder at 2.5 m (Gk1, Bl4), and its 5–20 m
  layer is 0.8–1.8 °C too cold in June–September, more than Norkyst's except in July. In July
  Norkyst and NorFjords are as cold or colder there (−1.9 and −3.7 °C against −1.8) without the hot
  top metre, so that layer's cold is not only heat held at the surface. Before July
  at most seven cells are above 25 °C on a day; the 29.6 °C cell of 23 November 2024 is unexamined.
  July 2026 repeats it in FjordSim alone: the inner-fjord casts hold 17–19 °C down to 4.5 m where
  FjordSim has 10–15 °C, and the top 5 m are 3.0 °C too cold (Norkyst −0.8); the warmest cell is
  35.3 °C on 17 July 2026.
- The `:hydrography` station writer's list is the 2014–15 programme's: SF-1 has no casts in
  2024–2026 and OF-1 and OF-5 only IMR's monthly ones (Torbjørnskjær, Breiangen vest), while the
  programme's SF-3, VAL-1, BO-1, I-1, I-4, KF-1 and SAN-5 have no hourly output and are scored from
  the daily snapshots.
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
  freezing on 6 Jan against its 473. In January–February 2026 the top 5 m are 2.3 and 1.5 °C too
  warm against the casts, where Norkyst is within 0.8 °C, while up to 46 905 cells are below −1 °C
  (15 February); whether the clamp causes it is unchecked. The slab-ice runs are archived as `oslofjorden_2026-11-06`,
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
  (346, 39, 24) is at −1.9 °C again on nine days of July–September 2025. In May–September 2026
  single cells reach −1.9 °C on 28 days, most often (208, 632, 32), (251, 649, 32) and (346, 39, 24).
- Salinity overshoots: the saltiest cell is above 36 g/kg on 157 days of October 2024 – September
  2025 (37.9 on 11 January, 37.3 in August) and on every day of October 2025 – September 2026 (38.0
  on 2 January 2026), against 35.1 in the initial state. In the Glomma spring flood of 22–24 May
  2025, up to 1821 cells in the Hvaler skerries (i 277–365, j 113–181), in columns one to three
  cells (1.6–3 m) deep under a 13–25 g/kg top cell, reach up to 37.6 g/kg. From 16 to 22 May one to
  five bottom cells there sit at freezing, e.g. (347, 166, 37) on 20 May: 35.2 g/kg and −1.2 °C
  under 30.1 g/kg and 10.9 °C, then 25.3 g/kg and 16.7 °C. A cell saltier and colder than all
  around it is what the centred fallback at its top face leaves when water leaves through that
  face; unconfirmed, as the snapshots hold no vertical velocity. Where the January and August
  maxima and those of 2026 sit is unexamined.
- (213, 418, 9), a 182 m pit in Breiangen with no wet neighbour at its level, climbs from 35.9 g/kg
  on 1 May 2025 to 36.05 at the end of July under 34.4–34.7 in the cell above; from June it is the
  saltiest cell on most days.

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
