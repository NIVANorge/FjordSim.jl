# TODO

Code changes found while validating `examples/oslofjorden_validation.jl` against the NIVA CTD casts,
NorKyst v3 and FjordOs. The problem statements quote the first run (days 0–298, 2014-04 to 2015-01).
The results under each task, and the sections from "Sill basins are not renewed" on, quote the full
640-day run with every fix below in it (stamp 20261014T142318, 2014-04-01 to 2015-12-31). Its report
and scripts are in `~/FjordSim_results/oslofjorden_validation/validation/report_full/`.

## Rivers as a mass source, not a relaxation

`add_rivers` writes rivers as a relaxation of T and S towards the river values over the plume
levels, with `λ = Q / V` (`src/Rivers/Rivers.jl`). The tracer dilution is right, but no volume
enters the domain:

- the river adds no volume, so sea level does not rise and there is no outflow driven by the river
  itself;
- the outflowing surface layer does not entrain water from below, so the estuarine circulation it
  should drive is missing;
- in the sill basins, without that circulation there is no deep inflow of saltier water to replace
  what the surface layer carries out. In the model, basin water is only diluted from above and is
  never renewed.

That fits what the run shows in the sill basins. Observed water at D-2 50 m stays at 31.0–31.4 g/kg,
while the model's falls from 32 to 27.7. In the Drammensfjord, FjordSim's upper 20 m is within
about 1 g/kg of the observations, but the deep basin freshens steadily. In the open fjord, the
surface is 2–4 g/kg saltier than NorKyst.

- [x] Implement rivers as a volume source.
      - Oceananigans allows exactly one volume hook without type piracy, `forcing.η` (m/s), so the
        water enters through the surface face of its columns; `LwSrc`-with-`Qshape` and `LuvSrc`
        would both need overrides of `compute_w_from_continuity!` or the immersed masking.
      - The tracer half is `w (c_r − c) / Δz` in the top cell, replacing the top cell's own value
        that the surface inflow carries in — NumericalEarth's correction, as a forcing.
      - Each river is spread over connected columns until its peak discharge enters no faster than
        `maximum_inflow_velocity` (1 mm/s in the setups): Glomma's eastern arm takes 47 columns.
      - Rivers moved to a file of their own (`rivers_nve.nc`, ~3 MB, daily over the run window)
        instead of a patched copy of `forcing.nc`; `FromForcing` no longer starts river cells at
        S = 0. Relaxation, plume depths and `OF800RiversConfig` are gone; optional
        `RiverMouthMixing` reuses NumericalEarth's mouth diffusivity.
      - "Rivers as a volume source" in `docs/architecture.md` has the details.
- [x] Check that the river cell's mass and salt budgets close, so that `∫S dV` changes only by
      boundary fluxes. Do it with a test on a small offline grid.
      - "river budget" in `test/runtests.jl`: on a closed basin `Σ η A` rises by `Q t` to ~5e-9, and
        a fresh river leaves the linear-free-surface salt content `Σ S V + Σ S_N η A` unchanged to
        ~2e-8.
- [x] Then re-check the D-2 and D-3 basin salinity trend and the near-surface salinity bias.
      - Near-surface bias fixed: the top 5 m is +0.1 g/kg on average (RMSE 4.2, FjordOs 5.8), and the
        25 m − 1 m salinity difference is unbiased overall. The Drammensfjord's is right (RMSE 1.1,
        FjordOs 4.5).
      - The basin trend is slower but not gone: D-2 at 50 m goes 32.1 → 30.0 g/kg over 640 days
        against an observed 31.0–31.4, and does not level off. See "Sill basins are not renewed".

## Vertical mixing is too strong

This probably applies to every setup and example: `oslofjorden` and `drammensfjorden` use
`CATKEVerticalDiffusivity(FT; minimum_tke = 7e-6)`. `examples/lista.jl` has moved to 7e-7. The
`oslofjorden()` comment shows that this TKE floor, not the prognostic TKE, sets κ below about 6 m.
It argues that 7e-6 is at the low end of the measured basin diffusivities, but the run contradicts
that:

- Ep1 at 50 m freshens from 32.8 to 30.2 g/kg and warms from 6.3 to 10.5 °C. NorKyst stays near
  33.3 g/kg and observations at 33.0–33.3.
- At Ep1 the salinity difference between the surface and 25 m is about 3.5 g/kg, against NorKyst's
  7. The water columns become more uniform over time, and summer heat reaches deeper than observed.
- In D-3, the 30 g/kg isohaline sinks from 14 m to below the 66 m bottom by October.

Tasks:
- [ ] Run a short sensitivity test, 60–90 days from the same start, with a lower `minimum_tke`.
      Measure the 50 m salinity trend at D-2 and Ep1.
      - 13 Oct, after the volume-source rivers and the new bathymetry: `minimum_tke = 7e-7` for 27
        days (`oslofjorden_validation_tke7e-7_30d`, cut short when GCE deleted the VM), against the
        first 27 days of the full 7e-6 run. Over the first 26 days, the 50 m trend per 30 days is
        D-2 +0.39 against +0.28 (the old run's −0.31), and Ep1 0.00 against −0.14 (old −0.55). The
        new rivers stopped the 50 m freshening at both settings, and the floor differs by about 0.1,
        too little to tell apart in 26 days.
      - The only casts inside that window are Ep1 and Dk1 on 9 April (day 8). Both settings are too
        salty in the top 5 m and too weakly stratified, and 7e-7 is closer. Salinity at 25 m minus
        1 m is 11.0 observed at Ep1, against 9.3 at 7e-7, 7.3 at 7e-6 and 6.8 in the old run. Dk1
        has 9.4 against 7.4, 6.2 and 5.8. At 1 m Ep1 is 20.1 observed, against 22.9 and 24.7. Below
        20 m all runs are about 0.3–1.0 too salty and about 2 °C too cold (Ep1 6.0–6.6 °C against
        8.2–8.5). The floor does not change that, so it is the starting state.
      - NorKyst's surface-to-25 m difference of about 7 at Ep1 is weaker than that cast, so it is
        not a reliable target for stratification in the inner fjord. Judge against the casts.
      - With 7e-7 the thin surface layer in shallow bays reached 22–23 °C by day 27, against about
        16–18 °C at 7e-6. No observations cover those bays.
      - Next: `examples/oslofjorden_validation.jl` now sets 7e-7 itself, and the full validation is
        to be rerun with it. The 7e-6 run it replaces is archived as
        `oslofjorden_validation_2026-10-13` (27 days). Check it against the casts on 8 May (day 37,
        inner fjord) and on 16 June and 5 July (days 76 and 95, the first at D-2 and D-3), and then
        against the rest of the casts.
      - Done: the full 640-day run at 7e-7. Against 216 casts, the salinity RMSE is 1.73 g/kg
        (FjordOs 2.79) and the temperature RMSE 1.22 °C (FjordOs 2.34). Temperature is unbiased at
        every depth band, and the summer heat no longer reaches too deep in the Inner Oslofjord.
        Salinity is 0.4–1.4 g/kg too fresh below 5 m, and 1.2–2.8 too fresh at 5–20 m in every
        region: the brackish layer is a few metres too deep.
      - Costs of the lower floor: the 2–3 m deep cells overheat (see "Hot shallow cells"), and
        Åsgårdstrand's monthly variance is larger than observed in most months. Its monthly means
        are still within 0.6 °C in 17 of 21 months.
      - Keep 7e-7 in the example. Whether `oslofjorden()` should follow is the open rationale task
        below.
- [ ] Measure numerical mixing from WENO plus `PartialCellBottom` z-levels on steep slopes, for
      example from the tracer variance budget, and compare it with CATKE's κ.
- [ ] Revisit `HorizontalScalarBiharmonicDiffusivity`'s κ. NorKyst uses a 10 m² s⁻¹ harmonic κ on
      tracers.
- [ ] Update the rationale comments in the setups once this has been measured.

## Bathymetry is smoothed too much

The model basins are 15–40 m too shallow. The model column is the deepest wet cell within 6 km of
the station, from `bathymetry.nc`; the real column is the deepest NIVA cast there:

| Station | Deepest model cell (m) | Deepest cast (m) |
|---|---|---|
| Ep1 | 111 | 149 |
| Fl1 | 132 | 165 |
| Im2 | 166 | 201 |
| D-2 | 100 | 115 |

The channels are wrong too: the Svelvik sill is 5.5 m against 11 m real, and the Drøbak Sound is
at most 131 m deep against 208 m. The comment in the example says the slope limit
(`max_slope_factor = 0.25`) costs only 0.6 m at the Drøbak sill. That measures one stage, so the
loss must come from the others.

- [x] Run the `smooth_bathymetry_gaps!` pipeline stage by stage and record, after each stage, the
      basin maxima and the sill depths along the deepest path.
      - Most of the loss came *before* the pipeline, in the gridding of the soundings: `invdistnn`
        with `smoothing=0.2` in degrees weighted the 16 nearest soundings almost equally, a plain
        mean. That took Ep1 from 154 m (deepest sounding within 2 km) to 114 m and the Drøbak Sound
        from 211 m to 176 m.
      - Then the regrid (Fl1 −10 m, Hm4 −16 m, Svelvik 8.8 → 6.6 m), and `max_slope_factor = 0.25`
        (Fl1 −13 m, Hm4 −17 m). The snap costs at most a fifth of a layer, 4 m at Im2.
      - Table and numbers: "Bathymetry gridding and sills" in `docs/setups.md`.
- [x] Preserve the extremes where it matters: keep basin maxima when regridding (not only the cell
      mean), and keep the deepest connection across a sill rather than averaging it with the flanks.
      - Soundings are now gridded linearly on their Delaunay triangulation (`GEONORGE_GRIDDING_ALGORITHM`).
      - The new `preserve_sills` stage (`preserve_sill_depths`) restores sills along one path per
        lost connection, from the native grid's merge tree.
      - `oslofjorden()` now uses `raw_resolution_factor = 4` and `max_slope_factor = 0.5`, and all
        three setups use `preserve_sills = true`.
      - On the validation box: Ep1 153 m, Fl1 153, Im2 209, Svelvik 10.8 m, Drøbak 23.5 m.
      - Drøbak is still deeper than the real ~20 m. The stage only deepens, so it cannot correct a
        sill the regrid made too deep.
- [x] Below 132 m the faces are 33.5 m apart, so check whether the coarse deep `z_faces` add to the
      loss in basin depth.
      - Only through the snap, which costs at most 0.2 × 33.5 = 6.7 m.
      - The real `z_faces` problem was at the bottom: the faces had been fitted to the averaged
        Hvaler trench (395 m), but the soundings reach 463–470 m. Both Oslofjord grids now go to
        −467 m (26 levels).
- [ ] Re-run `prepare_bathymetry`, `prepare_forcing`, `add_rivers` and `prepare_boundaries` for
      `oslofjorden`, `drammensfjorden` and `examples/oslofjorden_validation.jl`. Existing checkpoints
      are invalid: the `z_faces` and the land mask have changed.
      - Done for `examples/oslofjorden_validation.jl`. It ran 640 days on the new bathymetry
        (`max_slope_factor = 0.5`, `preserve_sills`) without a failure, so the stability concern is
        settled for that box.
      - Model columns at the CTD stations are now within 15 m of the casts (Ep1 150 against 149 m,
        D-2 116 against 115), except Im2 (165 against 201 m) in the narrow Drøbak Sound. Svelvik is
        10.8 m, Drøbak 23.5 m.
      - Still to do for `oslofjorden` and `drammensfjorden`.
      - `examples/lista.jl`'s `z_faces` gate needs re-measuring on the linear gridding.
- [ ] Fill single-cell pits on open edges. The new validation bathymetry has five on the south edge:
      wet edge cells whose east, west and inward neighbours are all below the seabed, so their only
      side neighbour is the open boundary. On the 13 Oct full run, the two thinnest went to the
      freezing point and were held there by the limiter: (255, 1, 17), a 3.6 m partial cell, at
      day 25, and (277, 1, 16), 2.1 m, at day 27. Neither the boundary data nor the rivers caused
      it; the boundary feeds 2.4–19 °C there. Unbounded WENO for T undershoots in such a small cell,
      walled on three sides beside the edge. It is one cell at a time, inside the boundary sponge.
      - Add a stage that raises such an open-edge column to the depth of its deepest neighbour, then
        re-run the four prepares above. Bounded WENO for T would hide it without removing the pits.
      - The other three pits are i = 76 and 84 at k = 11 (6 m of water) and i = 2 at k = 15 (8 m).
      - Full run: the domain minimum is below −1.5 °C on 591 of 640 days, always on the southern row.
        The cells are (84, 1, 11), (163, 1, 8), (144, 1, 9), (92, 1, 11), (76, 1, 11),
        (277, 1, 16), (94, 1, 11) and (255, 1, 17). So i = 92, 94, 144 and 163 are pits too, and the
        stage has to find them by the rule rather than from a list.

## Open-boundary tide

NorKyst v3's boundary `zeta` gives M2 17.7 cm and S2 1.4 cm, against 12.4 cm and 3.1 cm observed at
Viker. FjordSim inherits this error. The M2 current at Km1 is 5.6 cm/s, against 1.8 cm/s observed.

Full run, harmonic fits from 2014-05-01 to 2015-12-31, amplitudes in cm:

| | Boundary (NorKyst) | Viker model / obs | Oscarsborg model / obs | Oslo model / obs |
|---|---|---|---|---|
| M2 | 17.7 | 17.3 / 12.3 | 20.3 / 14.3 | 21.5 / 14.9 |
| S2 | 1.4 | 1.4 / 2.9 (34° late) | 1.7 / 3.4 (22° late) | 1.8 / 3.5 (22° late) |
| N2 | 3.4 | 3.3 / 3.0 | 3.9 / 3.5 | 4.1 / 3.6 |
| O1 | 1.1 | 1.2 / 1.8 | 1.2 / 1.9 | 1.2 / 1.9 |
| K1 | 1.0 | 1.0 / 0.4 | 1.1 / 0.6 | 1.2 / 0.7 |

- The error comes in at the boundary; the interior is right. The model's Viker tide is the
  boundary's. From Viker to Oscarsborg, M2 amplifies by 1.17 in the model and 1.16 observed, and
  S2 by 1.21 and 1.17. The model/observed ratio is the same at every gauge (M2 1.41–1.44, S2
  0.48–0.51), so one correction per constituent fixes them all.
- The boundary tide is nearly uniform along the 60 km southern edge (M2 17.2–17.8 cm, phase
  within 1.3°). One correction set therefore serves the whole edge.
- With the new bathymetry, the depth-mean M2 current at Km1 is 3.2 cm/s against 1.8 observed.
  The boundary's 1.4× accounts for most of the excess.
- The +0.19–0.22 m mean offset at the gauges is NorKyst's own mean `zeta` at the boundary
  (+0.14–0.17 m), a datum difference rather than a tide error. The non-tidal residual is good
  (r = 0.90), but its standard deviation is 14–16 % too small and the largest surges are
  18–25 % too low. That is worth a look at the boundary's subtidal `zeta` too, separately.

Ways to correct it, in order of preference:

1. **Rescale NorKyst's tide per constituent, after the regrid.**
   - Add a post-pass to `prepare_boundaries`, through a new boundary-config hook with a no-op
     fallback.
   - Fit the main constituents (M2, S2, N2, K2, K1, O1) to `eta`, `ubar` and `vbar` at every
     boundary node over the prepared window. The file spans two years, enough to separate S2 from
     K2.
   - Replace each fitted constituent with `G_k` times itself, where `G_k` is a complex gain: an
     amplitude factor and a phase shift. Everything else is left unchanged: the mean, the surges
     and the subtidal exchange.
   - The velocities must get the same `G_k` as the sea level. Flather (`GravityWaveRadiation`)
     combines the exterior transport and the exterior elevation, so correcting only `eta` would
     leave the two inconsistent and set the wrong incoming wave.
   - Gains from this run, observed/model at Viker: M2 ×0.71 (+3°), S2 ×2.1 (−34°), N2 ×0.91,
     O1 ×1.5, K1 ×0.4. State them in the setup as constants.
   - `harmonic_fit` lives in `Validation`, which is included after `Boundaries`, so it has to move
     down (to `Utils`) for the prepare step to use it.
   - Cost: the code plus a re-run of `prepare_boundaries`; no new download. The gains belong to
     the NorKyst v3 hindcast in 2014–15. The operational archive `oslofjorden()` reads needs its
     own check.
2. **Force the tide from a tidal atlas** (TPXO10-atlas or FES2022).
   - Remove NorKyst's tide with a 25 h low-pass filter, which keeps only the changes slower than
     a day.
   - Add the atlas tide in elevation and transport, with nodal corrections: the 18.6-year
     modulation of each constituent's amplitude and phase.
   - It is independent of NorKyst's tide and covers every setup and period without gauges.
   - But both atlases need registration, the nodal-correction code (or a dependency) is new, and
     global atlases are weakest in the Skagerrak, where M2 is small. It would still need checking
     against Viker.
   - Worth it only if several setups or periods need the tide right.
3. **Impose the observed tide.**
   - Low-pass NorKyst as in 2, then add the tide fitted to the Viker record, plus Helgeroa's for
     the western end of the edge if it is fetched.
   - The velocities still need NorKyst's tidal currents scaled to match, so in practice this is 1
     with target values in place of gains, and 1 is simpler.

- [x] Implement option 1, in `examples/oslofjorden_validation.jl` only, since the gains belong to
      this dataset and box.
      - `TideCorrectedBoundaries` wraps the `NorKystHindcastBoundariesConfig` and forwards its
        hooks to it. Its `prepare_boundaries` method calls `correct_boundary_tide!` on the written
        file, so nothing in `src/` changed. The gains are `VIKER_TIDE_GAINS`, for M2, S2, N2, K2,
        O1 and K1.
      - Checked on a copy of the prepared file. At the node south of Viker, each corrected
        constituent has exactly its gain and shift: M2 17.7 → 12.6 cm, S2 1.4 → 2.9 cm. The mean,
        the other seven constituents and the residual are unchanged to four decimals. The 3D `v`
        moves by the `vbar` increment at every wet level (to 3e-8). A second call is refused (the
        `tide_correction` attribute). It takes 41 s on the CPU.
- [ ] Re-run `prepare_boundaries` for the example, then a 60-day run. Check M2 and S2 at the three
      gauges, and the Km1 current (expected about 2.3 cm/s from the M2 factor alone). The interior
      M2 phase error (11° ahead at Oscarsborg and Oslo, 3° at Viker) will remain.

## Sill basins are not renewed

The main error left after the fixes above. Inside the sills, the deep water is freshened from above
and never renewed, steadily over both years. Outside the sills (OF-1 at 200 m, OF-5 and Im2 at
100 m) there is no trend, and the salinity stays within about 1 g/kg of the casts. Values in g/kg:

| Station, depth | Apr 2014 | Jan 2015 | Dec 2015 | Observed |
|---|---|---|---|---|
| Ep1 (Bunnefjorden), 120 m | 34.24 | 31.94 | 30.86 | 33.6–33.7 (2014; no 2015 casts) |
| D-2 (Drammensfjord), 50 m | 32.06 | 31.63 | 30.02 | 31.0–31.4 (both years) |
| D-3 (Drammensfjord), 80 m | 33.50 | 31.89 | 30.52 | – |

- Ep1's 33 g/kg water leaves the basin by about day 270 and never comes back. The rate slows
  from −2.3 g/kg in the first nine months to −1.1 in 2015, but it does not stop.
- The Drammensfjord starts 1.2–1.4 g/kg too salty (the NorKyst initial state) and freshens past
  the observations in spring 2015. It ends 0.9 too fresh, and its bias below 50 m goes from
  +0.86 (2014) to −0.44 (2015).
- Vestfjorden (Fl1, Dk1) gets renewal pulses (around days 300 and 600) and freshens between them.
- FjordOs holds its basin water steady, but 0.5–2 g/kg too salty.
- The sill depths and the stratification above them are now right, so neither explains this.

- [ ] Measure the dense inflow over the Drøbak and Svelvik sills: transport and salinity across a
      section at each sill from the daily output, compared with the basin's loss of salt. Im2's
      deep water is 0.4 g/kg too fresh, so the inflowing water may be too light as well as too
      rare.
- [ ] Measure the numerical mixing on the basin slopes (the WENO/`PartialCellBottom` task above),
      in Bunnefjorden first.
- [ ] Check the starting state. NorKyst's basin water on 9 April is up to 1 g/kg off and 1–3 °C too
      cold. If the adjustment from it drives the drift, a few months' spin-up from a cast-based
      initial state in the basins would show it.

## Intermediate water too warm in the middle fjord

In summer and autumn, FjordSim carries the surface heat to 100–150 m at OF-5 (Breiangen) and Im2
(Drøbak Sound). It reaches about 10 °C at 100 m, against 7–9.5 °C in the casts, in both years.
OF-1 at 200 m near the boundary follows the casts, so the heat comes from inside the domain. FjordOs
keeps that water 1.5–2.5 °C too cold.

- [ ] Find out whether it is local vertical mixing or warm water advected in at mid-depth: compare
      κ from CATKE with the heat budget at OF-5 between June and October.

## Surface salinity near the large rivers

The top 5 m is too salty off Glomma, at I-1 (+5.1 g/kg) and S-9 (+3.4). It is too fresh at LA-1
beside Numedalslågen (−3.5, and too fresh at every depth there, −1.5 to −2.2) and in the Oslo
harbour (−2 to −3.1).

- [ ] Before changing physics, check each river's discharge and how it is spread over its inflow
      columns against the casts' surface salinity. Glomma's eastern arm takes 47 columns, which may
      put its fresh water east of where I-1 and S-9 see it. Check the gauge substitutions
      (Akerselva, Lysakerelva) and the catchment-normal fallbacks.

## Negative salinity in a Glomma delta column

On 152 of 640 days some cell is below −0.1 g/kg, and on 35 days below −0.5; the minimum is −1.35
(31 October 2014). In every snapshot checked it is the same column, (260, 132) at 10.89 °E,
59.19 °N: 5.5 m deep, walled by land on two sides, about 10 cells downstream of Glomma Vesterelva's
inflow columns. The undershoot is at k = 23, 4.6 m down. It comes in the flood periods (June 2014,
October–November 2014, June and September 2015). Bounds-preserving WENO bounds only the advective
reconstruction, so another term drives this cell below zero.

- [ ] Reproduce it on a short run over June 2014 and find the term: the river tracer source, the
      vertical diffusion of fresh water, or the free-surface inflow interacting with advection in a
      nearly enclosed cell.

## Hot shallow cells

With the TKE floor at 7e-7, the warmest cell reaches 33.8 °C (27 July 2014) and 29.4 °C in 2015.
More than 10 000 cells exceed 25 °C on 10 days, all in the July 2014 heat wave, in 2–3 m deep
columns. No observations cover those bays.

- [ ] Add mixing in very shallow columns (like `RiverMouthMixing`), or set a minimum depth for wet
      cells, rather than raising the TKE floor everywhere.

## Validation

- [ ] `validate_simulation` reads only the newest run stamp (`station_files` takes
      `first(candidates)`). A resumed run splits its output across stamps, so it should concatenate
      the segments in time order.
- [ ] Add the NIVA CTD casts to the example as `CsvObservations`. The files in
      `observations/niva_ctd/` already use that schema.
- [ ] Add CTD-profile scoring to `src/Validation`:
      - pair each cast with the nearest wet column at matching depth;
      - compute bias and RMSE by depth band;
      - draw profile and basin-drift figures.

      The one-off scripts `ctd_validation.jl` and `ctd_fjordos.jl` in
      `<results_root>/validation/report_figures/scripts/` show the method. Comparing against FjordOs
      (MET THREDDS, ROMS s-levels) could be an optional reference source.

      The full run's scripts in `validation/report_full/scripts/` are the most recent version.
      `cells.jl` picks the column, `extract.jl` reads the daily snapshots in one pass,
      `fjordos_fetch.jl` gets FjordOs, and `ctd.jl` scores. Lessons for the module:
      - Score on the samples every model reaches, and per region and per year as well as overall.
        The drift in the sill basins only shows in a by-year or by-date view.
      - Pick the scored column by depth match, not by distance alone. AquaMonitor's OF-1 lies in
        70–100 m of water while its casts reach 440 m, so OF-1 is scored at METreport 11/2017
        Table 2's position (now in the example, 7fe5642).
      - A FjordOs reader must screen for bad records. Its Bn1 column is about 4 g/kg from top to
        bottom on 2014-12-08, below zero at the surface. The report keeps Bn1 out of the totals.
      - THREDDS takes about 5 minutes per station column over 21 months. Cache per station and
        fetch several stations in parallel processes.
- [ ] Add station-scale water-level scoring of the non-tidal residual per season. The
      surge-amplitude deficit (14–16 % in standard deviation) is invisible in the overall r = 0.90.
- [ ] The Inner Oslofjord programme (AquaMonitor project 6302) has no 2015 casts in Nivabasen, so the
      inner basins' 2015 drift is unscored. Ask NIVA whether the 2015 Indre Oslofjord casts exist
      elsewhere.

## Run infrastructure

- [ ] The run kept all 91 checkpoints, about 62 GB, even though `cleanup = true` is set. Find out
      why, and make `gcp/fjordsim-gcp`'s results sync, which uses `rsync` without delete, stop
      re-uploading pruned checkpoints.
      - The full run did the same: 186 checkpoints, 130 GB of its 182 GB in the bucket. The sync
        also crashed once (`FileNotFoundError`) on a checkpoint pruned between listing and upload,
        so `cleanup` does delete locally. The bucket keeps them because the sync never deletes.
- [ ] The periodic results sync uploads files that are still being written. The tide-gauge JLD2
      files get a record every 10 model minutes, about every 1.4 s of wall time, so their bucket
      copies are cut off mid-write and unreadable until the final sync after the run ends.
      Twenty-six copies of `Oslo` taken over SSH all failed. The `Not enough bytes in stream` errors
      in the run log are the same thing. Copy each growing file to a snapshot first, or exclude
      station files from the periodic sync and upload them at exit.
- [ ] `fjordsim-gcp logs` runs `tail -f` and never returns. Add a `--lines N` mode that prints and
      exits.
