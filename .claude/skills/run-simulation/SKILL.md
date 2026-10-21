---
name: run-simulation
description: Checklist for running FjordSim pipeline steps (prepare_*, download_*, add_rivers, run_simulation, validate_simulation) locally or on the GCP VM, re-running prepares after a config change, archiving a previous run before a fresh rerun, and resuming. Use when asked to run, rerun, resume or validate a fjord simulation or prepare its data.
---

# Running a FjordSim simulation

The rules and commands are in `docs/running.md` (local) and `gcp/README.md` (cloud). Read the
relevant one first; this is the order to work through.

1. **Ask before launching data steps.** `prepare_*` and `download_*` take minutes to days and can
   download hundreds of GB. Confirm with the user which steps to run and where (laptop or VM).
2. **Find out what is stale.** Compare the config against what was prepared: a changed
   `start_date`/`stop_time`, or `z_faces`/`size`/a bathymetry knob, needs the prepare steps listed in
   `docs/running.md` ("What to re-run after a config change"), in that order.
3. **Check the results directory.** If `results_root` (locally, and `results/<name>/` in the bucket
   for a GCP run) already holds checkpoints:
   - a resume → `pickup = true` / `--resume`, nothing to archive;
   - anything else (a rerun, a test, a changed parameter — all with the same config file, edited,
     never a new one) → the previous run should have been archived as `<name>_<YYYY-MM-DD>` when it
     finished; if it was not, archive it first (`gcp/README.md`, "When a run finishes"); confirm with
     the user before deleting anything in the bucket.
4. **Run.**
   - Locally: `julia --project -m FjordSim <step> --config <setup.jl>`.
   - On GCP: check the VM is up, then
     `gcp/fjordsim-gcp run --config <setup.jl> --steps <a,b,...> [--gpu] [--resume] [--run-id <id>]`.
     Use `--dry-run` first when unsure what a launch will do. `add_rivers` needs `NVE_API_KEY`.
5. **Watch it.** `gcp/fjordsim-gcp status`, `gcp/fjordsim-gcp logs <run-id>`; stop any monitoring
   loop you started once it's no longer needed.
6. **After the run.** `gcp/fjordsim-gcp pull-results <name>`, archive the run (step 3), then
   `julia --project -m FjordSim validate_simulation --config <config>`. Record findings and open
   problems in `TODO.md`; remove the items the run completed.
