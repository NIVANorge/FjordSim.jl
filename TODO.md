# TODO

Code changes found while validating `examples/oslofjorden_validation.jl` against the NIVA CTD casts,
NorKyst v3 and FjordOs (days 0–298 of the run, 2014-04 to 2015-01). The numbers are from that
assessment.

## Rivers as a mass source, not a relaxation

`add_rivers` writes rivers as a relaxation of T and S towards the river values over the plume
levels, with `λ = Q / V` (`src/Forcing/rivers.jl`). The tracer dilution is right, but no volume
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
- [ ] Then re-check the D-2 and D-3 basin salinity trend and the near-surface salinity bias.

## Vertical mixing is too strong

This probably applies to every setup and example: `oslofjorden`, `drammensfjorden` and
`examples/lista.jl` all use `CATKEVerticalDiffusivity(FT; minimum_tke = 7e-6)`. The
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
      - Short GPU run first: `max_slope_factor = 0.5` has not been run since the stages that fixed
        the earlier instabilities went in.
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

## Open-boundary tide

NorKyst v3's boundary `zeta` gives M2 17.7 cm and S2 1.4 cm, against 12.4 cm and 3.1 cm observed at
Viker. FjordSim inherits this error. The M2 current at Km1 is 5.6 cm/s, against 1.8 cm/s observed.

- [ ] Add a way to correct the tidal part of the boundary sea level: either scale the harmonics to a
      reference gauge (as FjordOs did) or force from a tidal atlas. This should be a boundary-config
      hook, not an edit to the generic pipeline.

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

## Run infrastructure

- [ ] The run kept all 91 checkpoints, about 62 GB, even though `cleanup = true` is set. Find out
      why, and make `gcp/fjordsim-gcp`'s results sync, which uses `rsync` without delete, stop
      re-uploading pruned checkpoints.
