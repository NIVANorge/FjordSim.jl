# Running FjordSim on GCP

`run_simulation` needs a GPU and runs on the GPU VM. Data preparation — every `download_*` and
`prepare_*` step and `add_rivers` — runs on the CPU-only preprocessing VMs, not locally and not on
the GPU VM, which would bill a GPU the prepare steps never use. `download_atmosphere` alone is close
to 10000 OPeNDAP reads per simulated year, and the Geonorge bathymetry source is a 4.5 GB FileGDB.

This directory runs any pipeline step of any setup on shared GCE VMs — a GPU VM for the
simulation, and optionally CPU-only preprocessing VMs for the downloads and regrids — and moves data
between your machine and a Cloud Storage bucket. Nothing about a fjord or a step is baked into the
image: which setup and which steps are arguments to one launch.

## The deployment model

`fjordsim-gcp` creates, starts, stops and deletes the VMs itself with `gcloud compute instances`,
and everything after `up` goes over SSH with sudo:

- VMs are **long-lived and shared**, not one per run: one GPU VM, `VM_NAME`, and optionally the
  preprocessing VMs `PREP_VM_PREFIX-1` ... `-N`. Every subcommand that acts on a VM takes
  `--vm N|NAME` to pick one; the default is the GPU VM. Any other name is another GPU VM, for runs in
  parallel or on another GPU type: `up --vm NAME [--machine-type T] [--accelerator A]` creates it.
- There is **no startup script** — `run` ships a script over `scp` and detaches it with `nohup`.
- Nothing self-deletes. `down` shuts a VM down *from inside* (`sudo poweroff`); `delete` removes a
  VM and its disk. Both act on the one VM `--vm` names, or with `--prep` on every preprocessing VM.

`up` creates a missing VM from `config.env`, trying each zone in `ZONES` once, in order, and moving
on only when a zone has no capacity or no such machine type; it fails if none has. A stopped VM is
started in the zone it lives in, which every later call looks up rather than configures. Every VM
gets Ubuntu 22.04, a `pd-balanced` boot disk, OS Login and the `cloud-platform` scope, running as
`SERVICE_ACCOUNT` (the project's default compute service account if empty).

| VM | Created with |
|---|---|
| GPU, `up [--vm NAME]` | `GPU_MACHINE_TYPE`, `GPU_ACCELERATOR` if set, `GPU_DISK_SIZE`; `--machine-type` and `--accelerator` replace the first two for a VM being created |
| preprocessing, `up --prep N` | `PREP_MACHINE_TYPE`, `PREP_DISK_SIZE`; `-1` ... `-N`. Without `N`, `up --prep` starts every existing one. |

All are on-demand (standard) VMs. Size the GPU disk for the data directory plus the output written
between two syncs ("Watching a run"), and a preprocessing disk for one fjord's data directory — the
Geonorge FileGDB alone is 4.5 GB.

### What it needs

| Who | Needs |
|---|---|
| The operator | `roles/compute.instanceAdmin.v1` on the project (it includes OS Login with sudo); `iam.serviceAccounts.actAs` on the VMs' service account; object admin on the bucket; Artifact Registry push on the image's repository |
| The VMs' service account | Read and write on the bucket; `roles/artifactregistry.reader` on the image's repository. It must be in the VMs' own project where the organisation enforces `iam.disableCrossProjectServiceAccountUsage`; the bucket and the repository may be in others. |
| The project | The Compute Engine API enabled; GPU quota for the machine type in the zones' region, and in `GPUS_ALL_REGIONS`; a way in on port 22: the default network's `default-allow-ssh` firewall rule, or IAP TCP forwarding from `35.235.240.0/20` |

> **A 403 from `list` proves nothing.** `gcloud storage buckets list` and
> `gcloud artifacts repositories list` are *project-level* calls, and they are commonly denied to
> accounts that nonetheless have full access to a specific bucket and can push to a specific
> repository. Test what you will actually use: `gcloud storage ls "$BUCKET/"` for the bucket, and a
> real push (`gcp/fjordsim-gcp build`) for the registry. For a project or an instance,
> `testIamPermissions` on that resource answers exactly.

## How it fits together

```
your machine                    gs://BUCKET                     a VM (shared, long-lived)
────────────                    ───────────                     ─────────────────────────
~/FjordSim_data/<fjord>  ──push──►  data/<fjord>/  ──stage in──►  /mnt/stage/data/<fjord>
                                                                        │  bind mount
                                    configs/<run-id>/*.jl ──────►  container: FJORDSIM_DATA_ROOT=/data
                                                                        │
~/FjordSim_results/<fjord> ◄─pull─  results/<fjord>/ ◄─stage out─  /mnt/stage/results/<fjord>
```

The package never learns about Cloud Storage. `FJORDSIM_DATA_ROOT` and `FJORDSIM_RESULTS_ROOT` name
the parent of the per-fjord directories (see `fjord_data_root` in `src/Utils.jl`), the VM
bind-mounts staged directories onto them, and every path in every setup follows. Staging runs on the
*host*, which has `gcloud` installed and the VM's own service account — so the image carries no
cloud SDK and no credentials.

Object storage is not used as a filesystem. `add_rivers` copies a 360 MB file, `Checkpointer` lists
a directory, and the NetCDF writers append — all of which a FUSE mount does badly. Copy in, run,
copy out.

The remote run script syncs results and raw downloads to the bucket every ten minutes, so an
interrupted job loses minutes of work, not all of it ("Resuming after an interruption").

## One-time setup

### 1. Fill in `gcp/config.env`

```bash
cp gcp/config.env.example gcp/config.env
```

All of them facts about the deployment:

| Variable | What it is |
|---|---|
| `PROJECT_ID` | the project the VMs live in |
| `BUCKET` | `gs://…` — staged inputs, results and uploaded configs |
| `IMAGE_URI` | the full image URI; everything else about the registry follows from it |
| `VM_NAME` | the GPU instance |
| `PREP_VM_PREFIX` | the preprocessing VMs are `PREP_VM_PREFIX-1` ... `-N`; `--vm N` names the Nth. Empty without them. |
| `ZONES`, `GPU_MACHINE_TYPE`, `GPU_ACCELERATOR`, `GPU_DISK_SIZE`, `PREP_MACHINE_TYPE`, `PREP_DISK_SIZE`, `SERVICE_ACCOUNT` | what `up` creates the VMs with ("The deployment model") |

`IMAGE_URI` is written out in full, so the repository may be in any project or region.

Every variable is also overridable per invocation, so a one-off needs no edit at all:

```bash
BUCKET=gs://other-bucket gcp/fjordsim-gcp run --config examples/oslofjorden/oslofjorden.jl --steps run_simulation --gpu
```

### 2. Build and push the image

```bash
source gcp/config.env
gcloud auth configure-docker "${IMAGE_URI%%/*}"
gcp/fjordsim-gcp build
```

The image bakes the precompile cache for the whole Oceananigans / NumericalEarth / CairoMakie
stack. That build is slow once so that no VM pays 10-20 minutes of precompilation on boot, while
being billed for a GPU that is doing nothing.

## The loop

```bash
gcp/fjordsim-gcp up                                  # create or start the VM
gcp/fjordsim-gcp bootstrap                           # driver, Docker, toolkit, image pull
gcp/fjordsim-gcp run --config examples/drammensfjorden/drammensfjorden.jl --steps run_simulation --gpu
gcp/fjordsim-gcp logs <run-id>
gcp/fjordsim-gcp pull-results drammensfjorden
gcp/fjordsim-gcp down                                # stop billing; the disk survives
gcp/fjordsim-gcp delete                              # when done for good: the disk goes too
```

`up` is idempotent: it creates the VM if it is gone, starts it if it is stopped, and exits cleanly
if it is already running. A *stopped* VM can only be started in the zone it already lives in, so a
stockout there fails `up` even when another zone has capacity — `delete` it, and the next `up`
creates it wherever there is room.

`bootstrap` is idempotent too. On a plain Ubuntu 22.04 image it installs the NVIDIA driver (~5 min,
then a reboot it drives itself), Docker and the container toolkit; on a second call, or after
`down`/`up` or host maintenance, every stage is a no-op. After a `delete` it starts over.

### Prepare on the preprocessing VMs

Always pass `--vm N`: without it `run` goes to the GPU VM. The VMs are cheaper than the GPU VM and
there can be several, so the slow downloads of one setup run side by side, each on its own VM:

```bash
gcp/fjordsim-gcp up --prep 3                         # creates or starts -1, -2 and -3
for n in 1 2 3; do gcp/fjordsim-gcp bootstrap --vm $n; done

gcp/fjordsim-gcp run --vm 1 --config examples/moro/moro.jl \
    --steps prepare_bathymetry,download_initial_state,prepare_initial_state,add_rivers
gcp/fjordsim-gcp run --vm 2 --config examples/moro/moro.jl --steps download_boundary_data
gcp/fjordsim-gcp run --vm 3 --config examples/moro/moro.jl --steps download_atmosphere,prepare_atmosphere
# once vm 1 and vm 2 have succeeded (gcp/fjordsim-gcp status):
gcp/fjordsim-gcp run --vm 2 --config examples/moro/moro.jl --steps prepare_boundary_data

gcp/fjordsim-gcp down --prep                         # stops all of them
gcp/fjordsim-gcp delete --prep                       # deletes all of them, disks included
```

With a single preprocessing VM, the same steps go to `--vm 1` one job after another.

Prepare steps take no `--gpu`: `architecture = :auto` resolves to `CPU()`, so the same config serves
both VMs without an edit. Anything a prep job writes under the fjord's data directory is synced
back to the bucket when it finishes, so the next job finds it. Raw downloads in subdirectories are
also synced every ten minutes, once a file is complete, so an interrupted download resumes from the
bucket; prepared top-level files wait for the final sync. Every job other than a plain
`run_simulation` stages the whole data directory, and a job that succeeds deletes it from the VM's
disk once it is in the bucket, so the disk needs room for one job, not every job so far. Extra
fjords from `--stage` are neither synced back nor deleted. A prep job's results directory is not
synced: only a fresh `run_simulation` or a `--resume` owns what is there, and a prep job would
otherwise upload whatever an earlier run left on the VM's disk, back under a name that run was
archived from.

`add_rivers` needs `NVE_API_KEY` in your own environment — the same variable a local run uses. `run`
writes it to a mode-600 file on the VM rather than passing it as an argument, which would be
world-readable in `ps`. Free key: https://hydapi.nve.no/Users

`download_observations` with a Vannmiljø source — and a `download_initial_state` whose setup takes
casts from one — needs `VANNMILJO_API_KEY` the same way; `run` passes it whenever it is set.

A setup that reads another fjord's downloads by absolute path (`docs/setups.md`, Paths) needs that
fjord staged too:

```bash
gcp/fjordsim-gcp run --vm 1 --config path/to/myfjord.jl --steps prepare_bathymetry \
    --stage oslofjorden --stage-all
```

- **A job reads what is in the bucket when it starts**, and writes back when it finishes. A step
  that needs another's output (`docs/running.md`, "The steps": every ocean prepare needs
  `bathymetry.nc`) must start after the job producing it has succeeded, not merely started.
- **Parallel jobs on one fjord are safe** because a job uploads only the files it wrote — those
  changed after its stage-in — never its staged copies of the rest, which another VM may have
  replaced in the bucket meanwhile.
- **One job per VM at a time.** A job that succeeds deletes the fjord's data directory from its
  VM's disk, under any other job there for the same fjord.
- `bootstrap` installs Docker only — the driver stages are skipped on a VM without a GPU — so it
  takes a minute or two.
- `logs RUN_ID` finds the run's VM by itself (`run` records it in `runs/<run-id>/vm`).
- `down --prep` stops every preprocessing VM, including one in the middle of a job; a single one
  can be powered off from inside with `down --vm N`.

### A one-off experiment

`--config` is a `.jl` file, uploaded to the bucket together with every `.jl` file beside it and
fetched by the VM, so a config may `include` a sibling — as `examples/oslofjorden/oslofjorden_npzd.jl` includes
`oslofjorden.jl` — and trying a variant never means rebuilding the image. Name the file after the
fjord — the basename is what roots the staged directories.

```bash
gcp/fjordsim-gcp run --config examples/oslofjorden/oslofjorden_npzd.jl --steps run_simulation --gpu
```

### A setup with its own environment

A setup that loads packages FjordSim's environment does not have — a biogeochemistry package, YAML
— cannot run in the plain image. `build --environment DIR` builds a second image on top of it
(`gcp/Dockerfile.environment`) whose Julia environment is the one in `DIR`, and pushes that instead:

```bash
export IMAGE_URI=REGION-docker.pkg.dev/PROJECT/REPO/fjordsim-mysetup:latest
gcp/fjordsim-gcp build --environment path/to/mysetup_env
gcp/fjordsim-gcp run --vm 1 --config path/to/mysetup_env/mysetup.jl --steps prepare_bathymetry,add_rivers
```

- **Its own `IMAGE_URI`**, for the build and for every `run` against it: `run` pulls `IMAGE_URI` at
  the start of each job, so the two images live side by side and a plain setup keeps using the
  plain one.
- **`DIR` needs a `Project.toml` and a `Manifest.toml`**, instantiated. The build copies the git
  repository holding `DIR`, so a package the environment names by a relative `path` source inside
  that repository (`MyBiogeochemistry = {path = ".."}`) comes along, data files and all; a `path` source
  outside it does not.
- **FjordSim is this checkout's**, developed into the environment from the image's `/workspace`
  whatever source the environment names for it, so every step runs against the FjordSim the image
  was built from.
- **Rebuild after changing** FjordSim, a package the environment loads or a file such a package
  reads. The config file itself is uploaded with each run, as for any setup.

Every step and `--resume` run in the setup's environment. The GPU preflight runs in FjordSim's own,
which is what has CUDA as a direct dependency.

### Resuming after an interruption

Host maintenance can stop the GPU VM in the middle of a long run. Recovery is `up`, `bootstrap` (a
no-op, since the disk survived), then the same run id with `--resume`:

```bash
gcp/fjordsim-gcp run --config examples/oslofjorden/oslofjorden.jl --steps run_simulation --gpu --run-id oslo-2020
# after an interruption:
gcp/fjordsim-gcp up && gcp/fjordsim-gcp bootstrap
gcp/fjordsim-gcp run --config examples/oslofjorden/oslofjorden.jl --steps run_simulation --gpu --run-id oslo-2020 --resume
```

`--resume` stages the newest checkpoint in `results/<fjord>/` back down and sets `pickup = true`
(`gcp/resume.jl`), so the run continues from it rather than restarting. Nothing else is staged: the
resumed writers start files of their own, and the finished output stays in the bucket. The sync
keeps only the two newest checkpoints in the bucket — the newest, and one in case its upload is
incomplete — since every one of `mjosa`'s is 580 MB. This only works if the setup names
a `CheckpointWriter` — `examples/oslofjorden/oslofjorden.jl` and `examples/drammensfjorden/drammensfjorden.jl` both do.

Checkpoints carry no run tag (`docs/running.md`, "One run per results directory"), so
`fjordsim-gcp` refuses a run whose id is already marked running, and a fresh `run_simulation` while
`results/<fjord>/` in the bucket holds checkpoints — archive the earlier run first ("When a run
finishes"). A fresh run starts on the VM from an empty results directory.

### Watching a run

```bash
gcp/fjordsim-gcp status            # every run marker in the bucket
gcp/fjordsim-gcp logs oslo-2020    # tail on the VM, bucket copy once it is gone
gcp/fjordsim-gcp ssh               # the VM is a normal box; docker ps, nvidia-smi, etc.
```

Results sync up every 10 minutes while the run is in flight, so `pull-results` mid-run gives you
something to plot. A writer with `split_times` starts a new file part at each of them, and every
part whose successor exists is finished: the sync uploads it and deletes it from the VM's disk, so
the disk holds only the parts being written. A year of output in one file can fill a disk —
`mjosa`'s 3-hourly 3D snapshots reached 110 GB in 185 days — so give a long run's writers
`split_times` (monthly in `examples/mjosa/mjosa.jl`).

### When a run finishes

Archive it: move `results/<fjord>/` in the bucket to `results/<fjord>_<YYYY-MM-DD>/`, the UTC date of
the run's last output (`-2` if taken). The undated directory is then free for the next run — a rerun,
a test or a changed parameter alike — which uses the same config file, edited if need be, rather than
a new one. The config each run used is kept under `configs/<run-id>/`.

An object rename is a copy, so copy, check, then delete:

```bash
source gcp/config.env
from="$BUCKET/results/oslofjorden_validation"
to="$BUCKET/results/oslofjorden_validation_2026-10-16"
gcloud storage cp -r "$from/*" "$to/"
diff <(gcloud storage ls -L "$from/**" | grep -E 'Content-Length|Hash \(crc32c\)') \
     <(gcloud storage ls -L "$to/**"   | grep -E 'Content-Length|Hash \(crc32c\)') \
    && gcloud storage rm -r "$from"
```

A local copy from `pull-results` is renamed the same way (`docs/running.md`, "One run per results
directory").

## Cost

- **Nothing stops on its own.** A finished job leaves its VM up and billing; `down` (or
  `down --prep`) is a manual step.
- **A stopped VM's disk is still billed**, at the same rate as a running one's. `delete` ends that.
- Keep the bucket in the VMs' region — that traffic repeats all run long, and within a region
  staging is free. An image in another region is pulled across regions, once per VM, since later
  pulls reuse its cached layers. Pulling results to your laptop is egress — a few hundred MB per run.
- The stage-in filter is the other lever: a simulation pulls ~900 MB instead of ~5.5 GB because it
  never touches the raw sources.

## Reproducibility

`Manifest.toml` is git-ignored but *is* copied into the image, so the container resolves to exactly
the package versions on the machine that built it. The trade-off is that an image cannot be rebuilt
from a clean clone.

## Notes

- The `Scratch.jl` bathymetry cache (`src/Bathymetry/geonorge.jl`) lives in the Julia depot inside
  the image, so it is rebuilt on each VM. `raw_directory` is an ordinary config field; point it
  under `data_root` in a config file if you would rather it were staged and reused.
- The files here: `fjordsim-gcp` (the driver), `load_config.sh` (reads `config.env`, which
  `config.env.example` is the template for), `bootstrap.sh` and `remote-run.sh` (both rendered
  with `envsubst` and run on the VM), `Dockerfile` / `build_image.sh` (the image),
  `Dockerfile.environment` (a setup's environment on top of it), `preflight.jl` (fails a GPU run
  that silently fell back to CPU) and `resume.jl` (sets `pickup = true`).
