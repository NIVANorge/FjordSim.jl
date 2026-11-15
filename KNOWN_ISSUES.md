# Known issues

## `examples/oslofjorden/oslofjorden.jl` (160 m, 2024-10 to 2026-10)

Numbers are from the run of stamp 20261112T205524, scored to the end of April 2026, unless another
run is named. "Before" is the run of stamp 20261107T210205, which started from Norkyst's basin
water, flooded islands up to six cells, and ran WENO's default boundary fallback with no shortwave
penetration. Against the NIVA casts after the first month, S RMSE is 1.90 g/kg (1.81 before) and T
RMSE 1.58 °C (1.45): the overshoots below are gone, but the water mixes more.

- Gauge gaps, filled from the nearest date: Drammenselva's discharge (Mjøndalen bru) on 27 July –
  13 August and 20 August – 12 September 2025, Glomma's temperature (Sarpfossen) from 18 August
  2026 and Sandvikselva's from 4 August 2026 to the end. Lierelva (Oppsal) reads 0 m³/s from 6 May
  2026 to the end and Mossevassdraget (Moss dam) from 7 May to early September 2026, both
  uncontrolled values.
- IMR's Missingene is 110 m deep within a cell of the station against casts to 350 m, Tisler-5 31 m
  against 120 m (120 m water is a cell away), and Færder lies inside the southern boundary sponge.
- The Drammensfjord loses salt at 30–100 m. From November 2024 to March 2026 D-2 falls from 28.4 to
  24.7 g/kg at 30 m (casts 29.6–30.6), from 30.8 to 26.7 at 50 m (casts 31.0–31.4) and from 31.6 to
  29.8 at 100 m (casts 31.6). Before, it lost 1.1–1.5 g/kg at those depths. The Svelvik sill, a
  one-cell channel 12–14 m deep, is the same in both runs. First suspect: the first-order upwind
  fallback (`minimum_buffer_upwind_order = 1`), the one change that adds mixing; the run also
  changed the initial state and the shortwave, so it is unconfirmed. The setup now keeps the
  fallback for T and e only, untested in a run; S's centred fallback is expected to bring back its
  single cells above 36 g/kg (157 days of October 2024 – September 2025 before, up to 37.9).
- `CastCorrectedInitialState` leaves 121 surface cell pairs between cast and Norkyst columns, 6.3
  g/kg apart (median), at the Drøbak and Svelvik sills and around inner-fjord pockets it does not
  reach; the cast profiles add 3200 small density inversions (≤ 0.43 kg/m³, against 1.31 without
  them).
- The surface is too salty and the halocline too weak. In November 2024 – January 2025 the top 5 m
  are 4.0 g/kg too salty against the NIVA casts (3.9 before), and the 25 m − 1 m salinity difference
  2.8 g/kg too weak (1.4 before); over November 2024 – April 2026 the top 5 m are +3.4 g/kg (+2.5)
  and the difference 5.7 g/kg against 9.3 observed (7.6). The inner fjord's surface misses the
  fresh events the run before caught. Norkyst has the same November–January error (+3.8 g/kg), and
  IMR's NorFjords 160 m hindcast (ROMS, GLS mixing) is saltier still: on the casts all three reach,
  the November top 5 m were +4.5 g/kg in FjordSim before, +5.2 in Norkyst and +6.7 in NorFjords. In
  the Inner Oslofjord the top 5 m are also 1.5 °C too warm (1.1 before), where Norkyst is within
  0.1 °C.
- From February to April 2025 the water at 20–50 m is 1.8–1.9 °C too cold (1.2–1.5 before). Norkyst
  and NorFjords are 0.7–1.3 °C too cold there too, NorFjords under AROME and GLS mixing, so it does
  not point at NORA3 or CATKE.
- The deep-water renewal of February 2026 is missed, and in March–April 2026 the water below 50 m is
  2.5 °C too cold (2.1–2.2 before). The casts show warm, salty water entering the inner basins: Fl1
  at 100 m reads 33.3–33.5 g/kg against FjordSim's 32.2–32.6, Ep1 at 120 m 33.0 g/kg and 9.1 °C
  against 31.9 g/kg and 5.9 °C. Outside the Drøbak sill (Im2, 100 m) the casts read 34.9 g/kg,
  FjordSim 33.8–34.5, Norkyst 34.2–34.4 and NorFjords 34.8–35.0. At the fjord mouth in March (IMR's
  Torbjørnskjær and Færder, 30–75 m) the casts read 6.7–8.3 °C and 34.0–34.9 g/kg, FjordSim
  4.8–5.2 °C and 32.5–33.6 g/kg, Norkyst 4.4–5.3 °C; in February the casts there held a cold
  brackish layer at 30 m the models lack. Norkyst and NorFjords are 1.8–2.2 and 1.5–1.9 °C too cold
  below 50 m too. Norkyst is the open boundary's source, so its Skagerrak water at 30–100 m is the
  first suspect. Before, the model renewed Vestfjorden and Steilene instead in early August 2026
  (Fl1 at 100 m 5.7 to 7.9 °C), after an inflow in July 2026 that raised the volume-mean salinity
  1.1 g/kg in three weeks; the casts show no inflow then.
- Glomma's fresh water stays too close to its mouths. In the top 1.5 m over November 2024 –
  September 2025, I-4 off Fredrikstad is 3.2 g/kg too fresh (casts 3.2 g/kg), while I-1, 8 km
  farther south in the Hvaler skerries, is 6.3 g/kg too salty (5.4 before; casts 11.4) and S-9 5.5
  (4.4; casts 19.5). Norkyst is 3.1 and 2.8 g/kg too salty at I-1 and S-9, NorFjords 5.1 and 3.1.
  At 5 m all three are about 7 g/kg too salty at I-1 outside February–April, so there the error is
  shared. Unchecked: the 2/3–1/3 split between Østerelva and Vesterelva, and how open the skerries'
  channels are on the 160 m grid.
- The non-tidal water level is too weak at all three gauges: as daily means over November 2024 –
  September 2025 it had 66–67 % of the observed standard deviation before, Norkyst 61–62 % and
  NorFjords, with another open boundary, 87–88 %. The missing surge is in Norkyst's sea level, which
  the open boundary carries in; the set-down of 9 December 2024 reached −0.25 m against −0.6 m
  observed. M4 is half the observed (0.7–1.8 cm against 1.2–2.3) and 31–47° early; the other
  constituents are within 0.5 cm, and all are within 0.1 cm of the run before. Over June 2025 – May
  2026 the daily means carried 79 % before (Norkyst 71–73 %, NorFjords 92 %), over the whole run
  74 % (Norkyst 67–69 %).
- The summer heat sits too high in the column, though less than before. In July 2025 the warmest
  cell is 28.4 °C (21 July) and up to 6023 cells exceed 25 °C, against 33.3 °C and 48 088 before. The casts of
  21–22 July read 24 °C uniform to 2.5 m; FjordSim has 22.9–24.2 °C at 0.5 m and 22.3–23.4 °C at
  2.5 m (Gk1, Bl4; 20.5–21.0 at 2.5 m before). The 5–20 m layer is 0.4–1.3 °C too cold in
  June–September (0.8–1.8 before). In July Norkyst and NorFjords are as cold or colder there (−1.9
  and −3.7 °C against −1.3), so that layer's cold is not only heat held at the surface. Below it,
  20–50 m is now too warm from July to November 2025, by 0.5–1.8 °C (−0.4 to +0.7 before). Before,
  July 2026 repeated the top-metre heat in
  FjordSim alone: the inner-fjord casts held 17–19 °C down to 4.5 m where FjordSim had 10–15 °C, and
  the top 5 m were 3.0 °C too cold (Norkyst −0.8); the warmest cell was 35.3 °C on 17 July 2026.
- Holmenfjorden (Cj1): the station writer places it in a column 7.9 m deep, 221 m from the station,
  against casts to 62 m; 54 m water lies 437 m away.
- Im2's cast replaces the Drøbak Sound's pocket below 104 m alone, under Norkyst's water, which at
  99 m is denser than the cast's top cell in 863 of its 903 columns, by at most 0.13 kg/m³.
  `TargetFill` doubles the cells denser than the one below them, 4746 to 11 278, the largest 0.72 to
  1.31 kg/m³.
- `oslofjorden_npzd.jl` runs 2024–2026 with nutrient profiles measured in 2020.
- No sea ice. `SlabSeaIce` (no drift) grew runaway ice on this setup, so it runs
  `FreezingLimitedOceanTemperature`, which clamps every level at freezing and discards the heat: no
  ice insulates the surface and no melt freshens it. In January 2025 the surface was 0.3–2.6 g/kg
  saltier than the slab-ice run where that run sat at freezing, and 1729 surface cells sat at
  freezing on 6 Jan against its 473. In January–February 2026 the top 5 m are 2.4 and 1.4 °C too
  warm against the casts, where Norkyst is within 0.8 °C, while up to 46 546 cells are below −1 °C
  in February; whether the clamp causes it is unchecked. The slab-ice runs are archived as
  `oslofjorden_2026-11-06`, `oslofjorden_2026-11-07` and `oslofjorden_2026-11-07-2`:
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
    which falls back to first-order upwind instead, 14 of those 16 are no longer cooled. The setup
    runs it for T, under the clamp, so whether it stops the slab ice's growth is untested. The ice also
    outlived the warm water because NumericalEarth's default `ThreeEquationHeatFlux` uses a fixed
    friction velocity of 0.002 m/s.
- (213, 418, 9), a 182 m pit in Breiangen with no wet neighbour at its level, gains salt steadily:
  from 35.2 g/kg in June 2025 to 36.06 at the end of April 2026. From June 2025 it is the saltiest cell, and from 10 April 2026 the only one above 36 g/kg.

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
