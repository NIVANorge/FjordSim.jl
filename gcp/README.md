# Running FjordSim on GCP

`run_simulation` needs a GPU. Everything else — the downloads and the regrids — runs anywhere, but
`download_atmosphere` is close to 10000 OPeNDAP reads per simulated year and the Geonorge
bathymetry source is a 4.5 GB FileGDB, so both halves are often better off in the cloud.

This directory runs any pipeline step of any setup on shared GCE VMs — a GPU VM for the
simulation, and optionally CPU-only preprocessing VMs for the downloads and regrids — and moves data
between your machine and a Cloud Storage bucket. Nothing about a fjord or a step is baked into the
image: which setup and which steps are arguments to one launch.

## The deployment model

The VMs come up one of two ways, picked by `gcp/config.env`:

- **Direct mode** — no launch job configured. `fjordsim-gcp` creates, starts, stops and deletes the
  VMs itself with `gcloud compute instances`: for an operator with Compute Engine rights in the
  project, such as their own.
- **Job mode** — launch jobs configured. For a project where the operator has no Compute Engine
  admin rights: someone who has them (the *admin*) sets up Cloud Run jobs that create or start
  fixed VMs, and the operator executes those jobs. `fjordsim-gcp` then never creates, starts, stops
  or deletes an instance itself, and never sets instance metadata.

Either way, everything after `up` goes over SSH with sudo, and is the same in both modes:

- VMs are **long-lived and shared**, not one per run: one GPU VM, `VM_NAME`, and optionally the
  preprocessing VMs `PREP_VM_PREFIX-1` ... `-N`. Every subcommand that acts on a VM takes
  `--vm N|NAME` to pick one; the default is the GPU VM.
- There is **no startup script** — `run` ships a script over `scp` and detaches it with `nohup`.
- Nothing self-deletes. `down` shuts a VM down *from inside* (`sudo poweroff`), which needs no
  `compute.instances.stop`; `down --prep` stops every preprocessing VM; `delete` removes a VM and
  its disk, where the operator may delete instances.
- A VM's zone is **discovered, never configured**: whoever created it was free to pick whichever
  zone had capacity.

### Direct mode

`up` creates a missing VM from `config.env`, trying each zone in `ZONES` in order and moving on only
when a zone has no capacity or no such machine type, and starts a stopped one in the zone it lives
in. Every VM gets Ubuntu 22.04, a `pd-balanced` boot disk, OS Login and the `cloud-platform` scope,
running as `SERVICE_ACCOUNT` (the project's default compute service account if empty).

| VM | Created with |
|---|---|
| GPU, `up` | `GPU_MACHINE_TYPE`, `GPU_ACCELERATOR` if set, `GPU_DISK_SIZE`; Spot, with `STOP` as the termination action: a preemption stops it and keeps its disk |
| GPU, `up --standard` | the same, on-demand |
| preprocessing, `up --prep N` | `PREP_MACHINE_TYPE`, `PREP_DISK_SIZE`, on-demand; `-1` ... `-N`. Without `N`, `up --prep` starts every existing one. |

A VM keeps the provisioning model it was created with: `delete` it to switch between Spot and
on-demand. Size the GPU disk for the data directory plus the output written between two syncs
("Watching a run"), and a preprocessing disk for one fjord's data directory — the Geonorge FileGDB
alone is 4.5 GB.

The operator needs `roles/compute.instanceAdmin.v1` on the project (it includes OS Login with
sudo), `iam.serviceAccounts.actAs` on the VMs' service account, the bucket and registry access in
"What the operator needs" below, and a way in on port 22: the default network's
`default-allow-ssh` firewall rule, or IAP TCP forwarding from `35.235.240.0/20`. The VMs' service
account needs what "What the admin provides" lists for it.

### Job mode

#### What the admin provides

| Piece | What it must do |
|---|---|
| GPU launch jobs, `LAUNCH_JOB_SPOT` and `LAUNCH_JOB_STANDARD` | Create the instance `VM_NAME` if it is gone, start it if it is stopped, and succeed if it is running. One NVIDIA GPU, Ubuntu 22.04, OS Login, and a boot disk sized as in direct mode. `up` executes the first, `up --standard` the second — conventionally Spot and on-demand. |
| Preprocessing launch job, `LAUNCH_JOB_PREP` (optional) | Task *i* (from 1) creates or starts `PREP_VM_PREFIX-i`. No GPU, Ubuntu 22.04, OS Login, and a boot disk sized as in direct mode. |
| Preprocessing stop job, `STOP_JOB_PREP` (optional) | Stops every preprocessing VM. |
| The VMs' service account | Read and write on the bucket; `roles/artifactregistry.reader` on the image's repository; the `cloud-platform` scope. |

#### What the operator needs

| Permission | On |
|---|---|
| `run.jobs.run` | every launch and stop job |
| `run.executions.get` | the GPU launch jobs: `up` waits for the execution to finish |
| `compute.instances.osAdminLogin` — SSH **with sudo** | every VM |
| `compute.instances.get`, `.list` | the project: zone discovery and status |
| `roles/storage.admin` (or object admin) | the bucket |
| `iam.serviceAccounts.actAs` | the VMs' service account, which OS Login requires to SSH into a VM that runs as one |
| Artifact Registry push | the image's repository |

Optional, and what lacking each costs:

| Permission | Without it |
|---|---|
| `run.jobs.runWithOverrides` and `run.jobs.get` on `LAUNCH_JOB_PREP` | `up --prep N` cannot set the number of VMs; `up --prep` starts the job's own task count. gcloud reads the job to validate any override, hence `get`. |
| `run.executions.get` on the preprocessing jobs | none: `up --prep` and `down --prep` launch `--async` and wait on the VMs' status instead |
| `compute.instances.delete` | `delete` fails, and a stopped VM's disk stays billed until the admin deletes it |
| `secretmanager.*` | none: the NVE key comes from the operator's environment |

> **A 403 from `list` proves nothing.** `gcloud storage buckets list` and
> `gcloud artifacts repositories list` are *project-level* calls, and they are commonly denied to
> accounts that nonetheless have full access to a specific bucket and can push to a specific
> repository. Test what you will actually use: `gcloud storage ls "$BUCKET/"` for the bucket, and a
> real push (`gcp/fjordsim-gcp build`) for the registry. For a job or an instance,
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

**A VM's disk may not be durable.** A launch job that creates its VM with
`--instance-termination-action=DELETE` — usual for Spot, and for a `--max-run-duration` cap — loses
the instance *and its boot disk* to a preemption or to that cap. A guest `poweroff`, or a direct-mode
Spot VM's preemption, is different — that only stops it, and the disk survives. This is why the remote run script syncs results and raw
downloads to the bucket every ten minutes, and why `bootstrap` has to run again after a preemption
but not after a `down`/`up`.

## One-time setup

### 1. Fill in `gcp/config.env`

```bash
cp gcp/config.env.example gcp/config.env
```

All of them facts about the deployment:

| Variable | What it is |
|---|---|
| `PROJECT_ID` | the project the VMs (and jobs) live in |
| `BUCKET` | `gs://…` — staged inputs, results and uploaded configs |
| `IMAGE_URI` | the full image URI; everything else about the registry follows from it |
| `VM_NAME` | the GPU instance |
| `PREP_VM_PREFIX` | the preprocessing VMs are `PREP_VM_PREFIX-1` ... `-N`; `--vm N` names the Nth. Empty without them. |
| `JOB_REGION`, `LAUNCH_JOB_SPOT`, `LAUNCH_JOB_STANDARD`, `LAUNCH_JOB_PREP`, `STOP_JOB_PREP` | job mode: the jobs and their region, named by the admin. All empty for direct mode. |
| `ZONES`, `GPU_MACHINE_TYPE`, `GPU_ACCELERATOR`, `GPU_DISK_SIZE`, `PREP_MACHINE_TYPE`, `PREP_DISK_SIZE`, `SERVICE_ACCOUNT` | direct mode: what `up` creates the VMs with ("Direct mode") |

Access to Artifact Registry is granted per *repository*, so the one you can push to may sit in
another project or region entirely. `IMAGE_URI` is therefore written out in full rather than
assembled from a project and a region, and both `bootstrap` and the remote run script take the
Docker auth host from it — which is all it takes to make the split work. The VMs' service account
then needs `roles/artifactregistry.reader` on *that* repository.

Deliberately absent: `ZONE` (discovered from the running instance), and anything about the NVE key
(`run` takes it from your environment). In job mode the machine, disk and service account belong to
the launch jobs, and `bootstrap` reads the service account from the metadata server when it needs
to name one.

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

`up` is idempotent and does the right thing in every state: it creates the VM if it is gone, starts
it if it is stopped, and exits cleanly if it is already running. It waits until the VM is running,
which can take long when the zones are short of GPUs. A *stopped* VM can only be started in the
zone it already lives in, so a stockout there fails `up` even when another zone has capacity —
`delete` it, and the next `up` creates it wherever there is room.

`bootstrap` is idempotent too. On a plain Ubuntu 22.04 image it installs the NVIDIA driver (~5 min,
then a reboot it drives itself), Docker and the container toolkit; on a second call, or after
`down`/`up`, every stage is a no-op. After a `delete`, or a job-mode preemption that deleted the
disk, it starts over.

### Prepare locally, simulate in the cloud

Preparing on a laptop is fine and often simplest; only the simulation really needs the cloud.

```bash
julia --project -m FjordSim prepare_bathymetry --config examples/drammensfjorden/drammensfjorden.jl
julia --project -m FjordSim add_rivers         --config examples/drammensfjorden/drammensfjorden.jl

gcp/fjordsim-gcp push-data drammensfjorden          # only the top-level *.nc
gcp/fjordsim-gcp run --config examples/drammensfjorden/drammensfjorden.jl --steps run_simulation --gpu
gcp/fjordsim-gcp pull-results drammensfjorden
```

`push-data` uploads only the prepared NetCDFs, because that is all a simulation opens. `--all` adds
the raw sources (`norkyst_hourly/`, `norkyst_initial/`, `nora3/`, the FileGDB), which only a
`prepare_*` step needs.

### Prepare in the cloud

Better for the slow steps, and it keeps the raw sources off your disk entirely.

```bash
gcp/fjordsim-gcp run --config examples/oslofjorden/oslofjorden.jl \
    --steps prepare_bathymetry,download_initial_state,prepare_initial_state,add_rivers,download_boundary_data,prepare_boundary_data
gcp/fjordsim-gcp run --config examples/oslofjorden/oslofjorden.jl \
    --steps download_atmosphere,prepare_atmosphere
gcp/fjordsim-gcp run --config examples/oslofjorden/oslofjorden.jl --steps run_simulation --gpu
```

Omitting `--gpu` hides the GPU from the container, and `architecture = :auto` resolves to `CPU()` —
the same config runs either way without an edit. Anything a prep job writes under the fjord's data directory is synced
back to the bucket when it finishes, so the next job finds it. Raw downloads in subdirectories are
also synced every ten minutes, once a file is complete, so a preempted download resumes from the
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

A setup that reads another fjord's downloads by absolute path (`docs/setups.md`, Paths) needs that
fjord staged too:

```bash
gcp/fjordsim-gcp run --config path/to/myfjord.jl --steps prepare_bathymetry \
    --stage oslofjorden --stage-all
```

### Prepare on the preprocessing VMs

The CPU-only VMs are cheaper than the GPU VM and there can be several, so the slow downloads of one
setup can run side by side, each on its own VM:

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
- **In job mode the number of VMs is the job's task count.** `up --prep N` overrides it, which
  needs `run.jobs.runWithOverrides` and `run.jobs.get` on the job ("What the operator needs");
  without them, `up --prep` starts the job's own count. Both `up --prep` and `down --prep` launch
  the job `--async` and wait for the VMs' status, so they need no access to its executions.
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
gcp/fjordsim-gcp run --config path/to/mysetup_env/mysetup.jl --steps prepare_bathymetry,add_rivers
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

### Resuming after a preemption

A Spot VM will interrupt a long run. Recovery is `up`, `bootstrap` (a no-op if the disk survived),
then the same run id with `--resume`:

```bash
gcp/fjordsim-gcp run --config examples/oslofjorden/oslofjorden.jl --steps run_simulation --gpu --run-id oslo-2020
# after a preemption:
gcp/fjordsim-gcp up && gcp/fjordsim-gcp bootstrap
gcp/fjordsim-gcp run --config examples/oslofjorden/oslofjorden.jl --steps run_simulation --gpu --run-id oslo-2020 --resume
```

`--resume` stages the newest checkpoint in `results/<fjord>/` back down and sets `pickup = true`
(`gcp/resume.jl`), so the run continues from it rather than restarting. Nothing else is staged: the
resumed writers start files of their own, and the finished output stays in the bucket. The sync
keeps only the two newest checkpoints in the bucket — the newest, and one in case its upload is
incomplete — since every one of `mjosa`'s is 580 MB. This only works if the setup names
a `CheckpointWriter` — `examples/oslofjorden/oslofjorden.jl` and `examples/drammensfjorden/drammensfjorden.jl` both do.

Checkpoint filenames carry no run tag and are scoped to `results_root`, so **two live runs must not
share a fjord**. With a single shared VM that is easy to do by accident. `fjordsim-gcp` refuses to
start a run whose id is already marked running.

Nor may a fresh run follow an earlier one into the same directory: the results sync never deletes, so
the earlier run's checkpoints would stay beside the new one's and a later `--resume` would pick up the
highest-numbered of them, likely the earlier run's. A fresh `run_simulation` (no `--resume`) is
therefore refused while `results/<fjord>/` in the bucket holds checkpoints — archive the earlier run
first ("When a run finishes", below) — and starts on the VM from an empty results
directory, clearing what an earlier run left on its disk.

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
  `down --prep`) is a manual step. A launch job's `--max-run-duration`, if it sets one, is the only
  backstop, and with `DELETE` as the termination action it deletes rather than stops.
- **A stopped VM's disk is still billed**, at the same rate as a running one's. `delete` ends that;
  in job mode it needs `compute.instances.delete`, which the admin may not have granted.
- Spot is roughly a 60-70% discount and is survivable here because of checkpoint resume.
- Keep the bucket in the VMs' region — that traffic repeats all run long, and within a region
  staging is free. An image in another region is pulled across regions, once per VM, since later
  pulls reuse its cached layers. Pulling results to your laptop is egress — a few hundred MB per run.
- The stage-in filter is the other lever: a simulation pulls ~900 MB instead of ~5.5 GB because it
  never touches the raw sources.

## Reproducibility

`Manifest.toml` is git-ignored but *is* copied into the image, so the container resolves to exactly
the package versions on the machine that built it. The trade-off is that an image cannot be rebuilt
from a clean clone. If you want CI builds later, drop `/Manifest*.toml` from `.gitignore` and commit
it — nothing else has to change.

## Notes

- The `Scratch.jl` bathymetry cache (`src/Bathymetry/geonorge.jl`) lives in the Julia depot inside
  the image, so it is rebuilt on each VM. `raw_directory` is an ordinary config field; point it
  under `data_root` in a config file if you would rather it were staged and reused.
- `--dry-run` prints the launch command and the fully rendered remote script without touching the
  VM. It is the fastest way to see what a launch will actually do.
- The files here: `fjordsim-gcp` (the driver), `bootstrap.sh` and `remote-run.sh` (both rendered
  with `envsubst` and run on the VM), `Dockerfile` / `build_image.sh` (the image),
  `Dockerfile.environment` (a setup's environment on top of it), `preflight.jl` (fails a GPU run
  that silently fell back to CPU) and `resume.jl` (sets `pickup = true`).
