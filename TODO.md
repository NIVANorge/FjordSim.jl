# TODO

Code changes found while validating `examples/oslofjorden_validation.jl` against the NIVA CTD casts,
NorKyst v3 and FjordOs. Unless a section says otherwise, numbers quote the full 640-day run (stamp
20261014T142318, 2014-04-01 to 2015-12-31), archived as `oslofjorden_validation_2026-10-16`; its
report and scripts are in that results directory under `validation/report_full/`.

## Vertical mixing is too strong

This probably applies to every setup and example: `oslofjorden` and `drammensfjorden` use
`CATKEVerticalDiffusivity(FT; minimum_tke = 7e-6)`. `examples/lista.jl` has moved to 7e-7. The
`oslofjorden()` comment shows that this TKE floor, not the prognostic TKE, sets κ below about 6 m.
It argues that 7e-6 is at the low end of the measured basin diffusivities, but the first run at 7e-6
(days 0–298, 2014-04 to 2015-01) contradicts that:

- Ep1 at 50 m freshens from 32.8 to 30.2 g/kg and warms from 6.3 to 10.5 °C. NorKyst stays near
  33.3 g/kg and observations at 33.0–33.3.
- At Ep1 the salinity difference between the surface and 25 m is about 3.5 g/kg, against NorKyst's
  7. The water columns become more uniform over time, and summer heat reaches deeper than observed.
- In D-3, the 30 g/kg isohaline sinks from 14 m to below the 66 m bottom by October.

Tasks:
- [ ] Measure numerical mixing from WENO plus `PartialCellBottom` z-levels on steep slopes, for
      example from the tracer variance budget, and compare it with CATKE's κ.
- [ ] Revisit `HorizontalScalarBiharmonicDiffusivity`'s κ. NorKyst uses a 10 m² s⁻¹ harmonic κ on
      tracers.
- [ ] Decide whether `oslofjorden()` and `drammensfjorden()` follow the example to 7e-7, and update
      the rationale comments in the setups. The 640-day run at 7e-7 scores salinity RMSE 1.73 g/kg
      and temperature RMSE 1.22 °C against 216 casts (FjordOs 2.79 and 2.34), unbiased in
      temperature, but 0.4–1.4 g/kg too fresh below 5 m and overheats the 2–3 m deep cells (see
      "Hot shallow cells").

## Bathymetry

The gridding and sill fixes are in (linear gridding of the soundings, `preserve_sills`,
`raw_resolution_factor = 4`, `max_slope_factor = 0.5`, faces to −467 m); see "Bathymetry gridding and
sills" in `docs/setups.md`. On the validation box the model columns at the CTD stations are within
15 m of the casts except Im2 (165 against 201 m) in the narrow Drøbak Sound, and Drøbak is 23.5 m
against the real ~20 m.

- [ ] Re-run `prepare_bathymetry`, `prepare_forcing`, `add_rivers` and `prepare_boundaries` for
      `oslofjorden` and `drammensfjorden` on the linear gridding and `preserve_sills` (done for the
      validation example). Their existing checkpoints are invalid.
- [ ] Re-measure `examples/lista.jl`'s `z_faces` gate on the linear gridding.
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

The correction is implemented in `examples/oslofjorden_validation.jl` only: `TideCorrectedBoundaries`
wraps the hindcast boundary config and rescales each constituent by `VIKER_TIDE_GAINS` after the
regrid (`correct_boundary_tide!`), velocities by the same complex gain as the sea level.

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
