# Known issues

## `examples/oslofjorden/oslofjorden.jl` (160 m, 2024-10 to 2026-10)

Numbers are scored against the NIVA casts from November 2024 to September 2026; October 2024 is the
adjustment from the initial state. The salinity RMSE is 1.59 g/kg and the temperature RMSE 1.52 °C.
On the samples NorKyst v3 160 m has, FjordSim's are 2.02 g/kg and 1.42 °C, NorKyst's 2.23 g/kg and
1.34 °C.

- Gauge gaps, filled from the nearest date: Drammenselva's discharge (Mjøndalen bru) on 27 July –
  13 August and 20 August – 12 September 2025, Glomma's temperature (Sarpfossen) from 18 August
  2026 and Sandvikselva's from 4 August 2026 to the end. Lierelva (Oppsal) reads 0 m³/s from 6 May
  2026 to the end and Mossevassdraget (Moss dam) from 7 May to early September 2026, both
  uncontrolled values.
- IMR's Missingene is 110 m deep within a cell of the station against casts to 350 m, Tisler-5 31 m
  against 120 m (120 m water is a cell away), and Færder lies inside the southern boundary sponge.
- The Drammensfjord's deep water freshens. D-2 starts on the casts and falls from 29.7 to 28.0 g/kg
  at 30 m by March 2026 (casts 30.1–30.4), from 31.2 to 30.0 at 50 m by September 2026 (casts
  31.0–31.4) and from 31.6 to 31.0 at 100 m (casts 31.6). The Svelvik sill is a one-cell channel
  12–14 m deep. Whether mixing below the halocline or too little salt water over the sill causes it
  is unknown.
- The February 2026 inflow to the Inner Oslofjord is missed. The casts show salty, warm water at Dk1
  (60 m) and Fl1 (100 m) from February; on 9 February FjordSim has 32.64 g/kg at Dk1 against 33.45
  observed (NorKyst v3 160 m 33.76, NorFjords 34.10), and from March to July both stations are
  2.0–3.2 °C too cold. Bunnefjorden is never renewed: at Ep1, 120 m, FjordSim falls from 32.9 to
  32.1 g/kg and from 8.4 to 6.0 °C, while the casts stay at 32.6–33.0 g/kg and 7.7–9.1 °C. Outside
  the Drøbak sill (Im2, 100 m) the February–June 2026 casts average 34.67 g/kg, FjordSim 34.10 and
  NorKyst v3 160 m 34.15. At sill depth outside (Im2, 18–28 m) FjordSim's density is within
  0.1 kg/m³ of the casts' in January–February 2026, and its excess over Fl1's deep water is larger
  (0.98 kg/m³ on 9 February against 0.61 observed). The water is diluted inside the sill: on 9
  February the casts at Hm4, 30–80 m, are as dense as Im2 at sill depth, FjordSim 0.85 kg/m³ lighter.
  The cause is convection where the inflow lands. The 3D snapshots show the crest's bottom cell as
  dense as the water outside on the inflow days. The water then flows off the step into the next,
  deeper column, lands on lighter water, and CATKE overturns that column to the bottom. On each pulse
  (4–14 January, 28 January–2 February, 10–15 February) the TKE in Hm4 below 25 m rises from about
  1e-7 to 1e-3 m²/s² with no shear, and the column from 36 to 99 m turns uniform within a day; Fl1
  below 50 m follows. The casts show the inflow as a bottom current under a stratified column.
  NorKyst and NorFjords follow the seabed with their levels, and NorKyst v3 160 m's inner basins are
  only 50–70 m deep. Finer z-levels would still overturn at every step; a downslope-flow
  parameterization would target it. The bathymetry adds a deeper route: the sill's east channel is
  21 m (about 20 m real), but the west channel, which the Drøbak jetty closes to about 6 m, is open at
  23.5 m. In March–May 2026 the water below 50 m is 1.9–2.3 °C too cold (NorKyst v3 160 m
  1.5–2.1 °C); NorKyst v3 800 m, the boundary source, has the same error near the boundary.
- `CastCorrectedInitialState` leaves 121 surface cell pairs between cast and Norkyst columns, 6.3
  g/kg apart (median), at the Drøbak and Svelvik sills and around inner-fjord pockets it does not
  reach; the cast profiles add 3200 small density inversions (≤ 0.43 kg/m³, against 1.31 without
  them).
- The surface is too salty and the halocline too weak. In November 2024 – January 2025 the top 5 m
  are 3.5 g/kg too salty, as in NorKyst v3 160 m on the same samples (NorFjords +5.4). Over the
  whole period the top 5 m are 2.1 g/kg too salty and S(25 m) − S(1 m) is 7.0 g/kg against 8.4
  observed. In November 2024 – January 2025 the Inner Oslofjord's top 5 m are also 0.9 °C too warm,
  where NorKyst v3 160 m is within 0.1 °C.
- From February to April 2025 the water at 20–50 m is 1.3–1.5 °C too cold. NorKyst v3 160 m and
  NorFjords are 1.1 and 0.9 °C too cold there too, NorFjords under AROME and GLS mixing, so it does
  not point at NORA3 or CATKE.
- In January–February 2026 the top 5 m are 1.8 and 1.2 °C too warm (NorKyst v3 160 m +0.1 and
  −1.3 °C, NorFjords +1.3 and +1.6 °C).
- Glomma's fresh water stays too close to its mouths. In the top 1.5 m over November 2024 –
  September 2025, I-4 off Fredrikstad is 3.4 g/kg too fresh (casts 3.4 g/kg), while I-1, 8 km
  farther south in the Hvaler skerries, is 6.7 g/kg too salty (casts 12.2) and S-9 4.3 (casts
  19.6). NorKyst v3 800 m is 3.1 and 2.8 g/kg too salty at I-1 and S-9, NorFjords 5.1 and 3.1. At
  5 m all three are 6–7 g/kg too salty at I-1 outside February–April, so there the error is shared.
  Unchecked: the 2/3–1/3 split between Østerelva and Vesterelva, and how open the skerries'
  channels are on the 160 m grid.
- The non-tidal water level is too weak at all three gauges: its daily means have 74 % of the
  observed standard deviation, NorKyst v3 67–69 % and NorFjords, with another open boundary, 89 %.
  The missing surge is in NorKyst v3 800 m's sea level, which the open boundary carries in; the
  set-down of 9 December 2024 reached −0.25 m against −0.6 m observed. M4 is 0.7–1.8 cm against
  1.2–2.3 cm observed and 32–48° early; the other constituents are within 0.3 cm.
- In summer the 5–20 m layer is too cold, by 0.3–1.1 °C in June–September 2025; in July NorKyst v3
  160 m and NorFjords are colder there (−2.0 and −3.5 °C against −1.1 on the same samples). The
  20–50 m layer is 0.5–1.1 °C too warm in September–November 2025, where NorKyst v3 160 m and
  NorFjords are within 0.2 °C. In July 2026 the top 5 m are 1.8 °C too cold (NorKyst v3 160 m
  −0.8 °C).
- Holmenfjorden (Cj1): the station writer places it in a column 7.9 m deep, 221 m from the station,
  against casts to 62 m; 54 m water lies 437 m away.
- Im2's cast replaces the Drøbak Sound's pocket below 104 m alone, under Norkyst's water, which at
  99 m is denser than the cast's top cell in 863 of its 903 columns, by at most 0.13 kg/m³.
  `TargetFill` doubles the cells denser than the one below them, 4746 to 11 278, the largest 0.72 to
  1.31 kg/m³.
- `oslofjorden_npzd.jl` runs 2024–2026 with nutrient profiles measured in 2020.

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

## Performance

- Output writing leaves the GPU idle about a quarter of the time (`oslofjorden.jl` on an NVIDIA
  L4). The pauses, about 1 s each, match the 10-minute `FieldStationWriter` and the hourly
  writers; this is by timing, not profiled.

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
