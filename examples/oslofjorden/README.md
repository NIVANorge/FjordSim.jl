# Oslofjorden

`oslofjorden.jl` runs the Oslofjord from 1 October 2024 to 1 October 2026, two whole winters, on a
424 × 681 × 39 grid of 160 m cells covering 10.00–11.20°E, 58.95–59.93°N:

- Geonorge Sjøkart Dybdedata bathymetry;
- an operational Norkyst-v3 initial state (the 160 m Oslofjord nest, filled out by the 800 m
  reference run), with NIVA casts in the sill basins;
- NVE rivers;
- hourly operational Norkyst-v3 800 m on the open southern edge, with its tide corrected to Viker
  and its upper layer's salinity to NIVA's FerryBox;
- NORA3 atmosphere with snowfall;
- water held at or above freezing, with no ice.

It writes daily 3D snapshots, an hourly 4 m level, hourly sea level at the three tide gauges and an
hourly section along the Drammensfjord. `oslofjorden_npzd.jl` adds an
NPZD model to it. Problems with the simulation are in [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md).

## Preparing and running it

```bash
julia --project -m FjordSim prepare_bathymetry     --config examples/oslofjorden/oslofjorden.jl
julia --project -m FjordSim download_observations  --config examples/oslofjorden/oslofjorden.jl  # VANNMILJO_API_KEY
julia --project -m FjordSim download_initial_state --config examples/oslofjorden/oslofjorden.jl  # VANNMILJO_API_KEY
julia --project -m FjordSim prepare_initial_state  --config examples/oslofjorden/oslofjorden.jl
julia --project -m FjordSim add_rivers             --config examples/oslofjorden/oslofjorden.jl  # NVE_API_KEY
julia --project -m FjordSim download_boundary_data --config examples/oslofjorden/oslofjorden.jl
julia --project scripts/measure_boundary_tide.jl examples/oslofjorden/oslofjorden.jl Viker
    # paste its gains into `VIKER_TIDE_GAINS`, and put the FerryBox files in place (below)
julia --project -m FjordSim prepare_boundary_data  --config examples/oslofjorden/oslofjorden.jl
julia --project -m FjordSim download_atmosphere    --config examples/oslofjorden/oslofjorden.jl
julia --project -m FjordSim prepare_atmosphere     --config examples/oslofjorden/oslofjorden.jl
julia --project -m FjordSim run_simulation         --config examples/oslofjorden/oslofjorden.jl  # a GPU
julia --project -m FjordSim validate_simulation    --config examples/oslofjorden/oslofjorden.jl
```

The FerryBox files come from NIVA's FerryBox API (https://ferrybox-api.niva.no, which needs a NIVA
service account) and go in `<data root>/observations/niva_ferrybox_2024-2026/`: one per ship and
month, `<ship>_<yyyymm>.csv`, with columns `time` (UTC), `lon`, `lat`, `temperature` and
`salinity`, one row per minute, only values that passed NIVA's QC and only while the ship moves.

## The domain

The box is the one MET Norway's ROMS model FjordOs CL covers (Røed et al., *A high-resolution,
curvilinear ROMS model for the Oslofjord*, METreport 4/2016; Hjelmervik et al., *Evaluation of the
FjordOs-model*, METreport 11/2017), and most station positions come from those reports. Beyond the
inner fjord it reaches east to 11.20°E for the Hvaler archipelago, both arms of Glomma's delta and
S-9; west to 10.00°E for Numedalslågen and LA-1; and south to 58.95°N, so LA-1 lies 7.7 km inside
the southern edge rather than in its 3.2 km sponge.

At 160 m the grid is finer than NorKyst-800 everywhere but about twice as coarse as FjordOs CL in
the sounds: the Drøbak Sound is 6–12 cells wide, and Svelvik, 180 m wide and 11 m deep, is one cell.
The sill-controlled tidal jets of METreport 4/2016 and the Drammensfjord's exchange through Svelvik
are out of reach whatever is tuned.

## Scoring it

`validate_simulation` fetches Kartverket's water level at the `:tidegauge` writer's stations itself.
The NIVA and IMR casts come from Vannmiljø, fetched by `download_observations`, one source per
contractor so each cast is scored once; each station takes the Vannmiljø location nearest it. The
casts are scored from the daily 3D snapshots. The Drammensfjord section is not an observation
programme: it shows the basin's vertical mixing and its renewal over the sill. Its stations are
every second cell of the deepest path through the 193 m grid, found by Dijkstra over wet cells with
a step cost of distance / depth².

## Why the values are what they are

Measured on this box at 193 m against the NIVA casts of 2014–15, unless stated.

### Grid and bathymetry

- **Levels.** A geometric stretch (ratio 1.25, capped at 33.5 m) from a 1 m surface cell to −467 m,
  with 41–132 m in eighteen 5.06 m levels. On the stretch's five 11–27 m levels there, WENO's
  implicit vertical diffusion (~|w|Δz) mixed the sill basins at 5–10× CATKE's κ; the split cut the
  Drammensfjord's salt-budget diffusivity below 41 m from 4.5e-5–3.8e-4 to 1.7e-5–7.0e-5 m²/s. The
  deepest sounding is 470 m; only 13 columns reach below −467 m, by at most 2.4 m. Re-check the
  deepest cell in `bathymetry.nc` after a change to the gridding or the box: too shallow truncates
  the basin silently, too deep shortens the barotropic substep for water that is not there.
- **Basin depths.** The deepest cell within 2 km of each station, in m:

  |                                        | Ep1 | Fl1 | Im2 | Hm4 | Svelvik |
  |----------------------------------------|----:|----:|----:|----:|--------:|
  | soundings                              | 154 | 166 | 211 | 119 | 11      |
  | inverse-distance gridding, factor 2    | 114 | 138 | 176 | 104 | 8.8     |
  | ... slope limit 0.25 and snap          | 111 | 114 | 166 |  71 | 5.5     |
  | linear, factor 4, 0.5, preserve_sills  | 153 | 153 | 209 |  99 | 10.8    |

  A ~97 m native grid holds Svelvik at 8.5 m, a ~48 m one at 11.3 m, which `preserve_sills` carries
  onto the model grid. The slope limit 0.25 cost 34 m at Fl1 and 24 m at Hm4, 0.5 costs 8 and 4.
- **The Drøbak jetty** (1874–79, Søndre Kaholmen to Hurum) is capped at 6 m along the line of 1–2 m
  soundings. Without it the sound's deepest crossing is 23.5 m through the channel it closes; with
  it the sill is the 21 m east channel (real ~20 m).
- **Islands** are all kept: flooding them added no open channel and removed Kaholmen.
- **`max_sounding_distance`** removes two places the FileGDB has no data for, which the gridding
  filled with water from soundings far away: north of Drammen (a 7 m sea) and the Swedish side
  south of Strömstad (100 m). Every real wet cell lies within 1.5 km of a sounding.

### Initial state

The 160 m nest reaches every sill basin and the fjord north of about 59.5°N; south of about 59.3°N
the state is the 800 m run's. The nest's own parent starts only on 2024-11-23; on 2024-10-01 the
seam between nest and 800 m run is 1.0 g/kg and 0.15 °C (median) at the surface.

The nest's sill basins are too shallow (z-levels end at 65–100 m against casts to 98–166 m), and the
water carried down below them is wrong: Ep1 at 120 m had 33.90 g/kg and 6.18 °C against 32.91 and
8.45 in the cast, D-2 at 100 m 33.53 g/kg against 31.65. `CastCorrectedInitialState` replaces the
basins' columns with the NIVA casts nearest the start (2024-09-20 and -26), which the next weeks'
casts match within 0.1 g/kg and 0.2 °C below 50 m. Im2's cast replaces only the Drøbak Sound's
pocket below 104 m: over whole columns it would leave 12.8 g/kg beside Norkyst's 24.5 at the surface.

### Rivers

Every mouth of 0.5 m³/s or more is discovered from NVE's ELVIS river network and sized by its REGINE
catchment, so the setup states no river coordinate; on a smaller box that was 21 mouths carrying
99.6 % of the inflow. A mouth runs at its catchment's normal unless an override attaches a gauge.
Glomma's two mouths share Solbergfoss's discharge 2/3 to Østerelva, 1/3 to Vesterelva, an
assumption: NVE describes the division and publishes no fraction. Akerselva's and Lysakerelva's
temperature gauges read 20–21 °C on average, a sun-exposed sensor, so those rivers carry the mean
of the gauged ones. `maximum_inflow_velocity` (1 mm/s) is a vertical Courant number of 0.18 on the
1 m top cell at the 3-minute time step.

### Open boundary

- **Tide.** Norkyst-v3's boundary M2 is 17.4 cm against 12.3 cm observed at Viker, and nearly uniform
  along the edge (17.2–17.8 cm, phase within 1.3°), and the model carries it unchanged to Viker. So
  one gain per constituent, measured at the boundary against Viker, corrects the whole edge, as
  FjordOs did (METreport 4/2016 §4.3).
- **Salinity.** Norkyst-v3 800 m lacks most of the brackish surface layer at the mouth. At 4.6 m the
  edge is 4.4 g/kg saltier than MS Color Fantasy's FerryBox where it crosses (10.66°E), by 2.1–6.1 in
  every month, and 6.8 saltier than MS Color Hybrid's (11.02°E). Against IMR's casts at Færder,
  1.8 km inside the edge, it is 4.7 g/kg too salty at the surface, 3.7 at 5 m, 2.5 at 10 m, 0.6 at
  20 m and right from 30 m down; temperature is within 1 °C. The ferry's excess matches Færder's over
  the top 5 m (4.2 g/kg), so the correction is full at the surface and tapers to none at 25 m.
  Corrected, the edge at Færder is within 0.6 g/kg of the casts at every depth, and the ships'
  excess at 4.6 m is 0.7 and 1.3 g/kg.
- **Closed walls.** The western wall's 14 wet cells lie inside the southern edge's sponge; the eastern
  wall's 26 are the Idefjord approach, a sill fjord. Opening either would make the downloaded band the
  whole domain, about 110 GB over two years against 5 GB.
- **Nudging.** 3 hours on inflow: a day relaxed 200× more weakly than the advection into the boundary
  cell at ~10 s steps.

### Model

- **Freezing.** In the mild winter of 2014–15 surface cells sat at freezing on 115 of 641 days, up
  to 5.6 % of the surface, all brackish. A non-drifting slab ice turned every heat deficit the
  advection left into ice that stayed (305 m in a one-cell channel at Larvik), so the water is only
  clamped at freezing.
- **Sponge.** Without it the boundary row carried 0.4 m/s velocity std and 2 m/s peaks against 0.13
  and 0.43 inside, and salinity roughness that ran away after 50 days. Viscous rather than a tracer
  relaxation, so it does not fight the boundary's own nudging. 20 m²/s is under the explicit limit
  (320 s on 160 m cells).
- **CATKE.** Below 6–9 m the prognostic TKE sits under the floor in most cells, so `minimum_tke` is
  the vertical mixing there. 7e-6 left the top 5 m too salty and too weakly stratified (25 m − 1 m
  salinity 7.3 at Ep1 against 11.0 observed); 7e-7 gave 9.3, at the cost of overheating 2–3 m deep
  cells. `Cᵇ` stays at 0.28: it scales the near-bottom mixing length, which Stigebrandt credited for
  the basins' high diffusivity.
- **Biharmonic.** ν₄ = 9.4e3 m⁴/s damps the 2Δx mode in 72 min and leaves an 8 km eddy alone; its
  explicit limit is 2180 s, which `AdaptiveTimeStep` does not watch. κ = 2e3 overshot at fronts,
  leaving 169 cells below 0 g/kg.
- **Tracer advection.** At river mouths WENO across the fresh front made water saltier than any
  source, so T and S are upwinded within three cells (WENO5's stencil) of them. Coastal cells are not
  upwinded: first-order upwind diffuses 8–40 m²/s at the surface. Plain WENO on S took cells beside
  Numedalslågen below TEOS10's −32 floor, so S is bounded. Near land and the seabed WENO falls back to
  `Centered(2)`, which leaves single cells at −1.9 °C or above 36 g/kg; T's fallback ends in
  first-order upwind instead, in the top four cells only, since below them it mixed heat into the sill
  basins (κ_T below 52 m at Bn1 2.4e-4 m²/s against 0.86e-4). S keeps the centred fallback: upwind
  faces mixed its halocline. Below the top four cells T's centred fallback leaves single cells up to
  25 °C in Glomma's salt wedge at 9 m.
- **Shortwave.** With all of it in the 1 m top cell, summer cells reached 29–33 °C and 5–20 m was
  0.8–1.8 °C too cold. The blue-green band's 3 m comes from the fjord's 4–6 m Secchi depths.
- **Downslope flow.** Each inflow pulse over the Drøbak sill in winter 2026 landed on lighter water
  in the next column and CATKE overturned it, so the February renewal arrived 0.85 kg/m³ too light at
  Hm4. A friction timescale of 1000 s gives 0.5 m/s for the sill water's 0.4 kg/m³ excess on the
  sound's 0.13 slope, which the scheme's cap limits to 0.22 m/s at the 3-minute step.

### Output

The daily 3D snapshots are 225 MB a record, 164 GB over the run; three-hourly would be 1.3 TB. They
are what the casts are scored on. The hourly 4 m level, about 100 GB, is the level the FerryBox
measures, and the surface and current maps are drawn from it. The tide gauges are hourly, six
samples in a period of MS4 (6.1 h), the shortest constituent analysed. The Drammensfjord section
is hourly for the tidal exchange over the Svelvik sill. There is no `StationWriter`: each of its
stations is an output writer of its own, written every interval, and the casts need no hourly series.

## Data the setup uses

### Rivers

NVE's gauges have gaps, which `add_rivers` fills from the nearest date: Drammenselva's discharge
(Mjøndalen bru) on 27 July – 13 August and 20 August – 12 September 2025, Glomma's temperature
(Sarpfossen) from 18 August 2026 and Sandvikselva's from 4 August 2026 to the end. Lierelva (Oppsal)
reads 0 m³/s from 6 May 2026 to the end and Mossevassdraget (Moss dam) from 7 May to early September
2026, both uncontrolled values (not yet quality-checked by NVE).

### Observation stations

- IMR's Missingene is 110 m deep within a cell of the station against casts to 350 m, Tisler-5 31 m
  against 120 m (120 m water is a cell away), and Færder lies inside the southern boundary sponge.
- Holmenfjorden (Cj1): the model column nearest it is 7.9 m deep, 221 m from the station, against
  casts to 62 m; 54 m water lies 437 m away.
- Vannmiljø places Ø-1 and Cj1 about 600 m from the positions in the setup.

### Biogeochemistry

`oslofjorden_npzd.jl` runs 2024–2026 with nutrient profiles measured in 2020.
