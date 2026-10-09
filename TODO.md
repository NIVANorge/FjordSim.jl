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

- [ ] Implement rivers the way ROMS does (`LuvSrc`/`LwSrc`): a volume transport `Q` through a cell
      face or into a cell, carrying the river's T and S. That needs a volume source in the free
      surface and the continuity equation. Find out what Oceananigans' `HydrostaticFreeSurfaceModel`
      and `SplitExplicitFreeSurface` allow before designing the hook. The tracer half already has an
      unused entry point: `ForcingFromFile` reads `λ > 1` as an x-flux and `λ < -1` as a y-flux
      (`forcing_term_x_flux`/`forcing_term_y_flux` in `src/Forcing/Forcing.jl`).
- [ ] Check that the river cell's mass and salt budgets close, so that `∫S dV` changes only by
      boundary fluxes. Do it with a test on a small offline grid.
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

- [ ] Run the `smooth_bathymetry_gaps!` pipeline stage by stage and record, after each stage, the
      basin maxima and the sill depths along the deepest path.
      - The stages are the regrid at `raw_resolution_factor = 2`, `interpolation_passes`,
        `close_narrow_passages`, `spike_ratio`, `minimum_cell_fraction` and `max_slope_factor`.
      - The sill depths can be computed with the maximin flood-fill used for the validation report.
      - Then find which stage removes the depth.
- [ ] Preserve the extremes where it matters: keep basin maxima when regridding (not only the cell
      mean), and keep the deepest connection across a sill rather than averaging it with the flanks.
- [ ] Below 132 m the faces are 33.5 m apart, so check whether the coarse deep `z_faces` add to the
      loss in basin depth.
- [ ] Any change here means re-running `prepare_bathymetry`, `prepare_forcing`, `add_rivers` and
      `prepare_boundaries`, and invalidates existing checkpoints.

## Open-boundary tide

NorKyst v3's boundary `zeta` gives M2 17.7 cm and S2 1.4 cm, against 12.4 cm and 3.1 cm observed at
Viker. FjordSim inherits this error. The M2 current at Km1 is 5.6 cm/s, against 1.8 cm/s observed.

- [ ] Add a way to correct the tidal part of the boundary sea level: either scale the harmonics to a
      reference gauge (as FjordOs did) or force from a tidal atlas. This should be a boundary-config
      hook, not an edit to the generic pipeline.

## Station output

- [ ] In `StationWriter`, interpolate u and v to cell centres. Today u is sampled on the west face
      and v on the south face, so a closed coastal face gives zero; Slagen's v is 0 at every level.
- [ ] The `model_depth` attribute is the depth of a level face, not the water depth `h`. Either write
      `h` as well or rename the attribute.
- [ ] In `examples/oslofjorden_validation.jl`, move these stations to AquaMonitor's positions:
      - Dk1 → 10.569384 E, 59.814999 N (now about 10 km off)
      - Ep1 → 10.723783 E, 59.786301 N (2.8 km off)
      - OF-1 → 10.6652 E, 59.0360 N (5 km off)

      Also check `INNER_OSLOFJORD_STATIONS` (marked APPROXIMATE) against
      `observations/niva_ctd/stations.csv`.

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
