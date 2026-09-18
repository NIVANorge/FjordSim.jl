# Running FjordSim on GCP

`run_simulation` needs a GPU. Everything else — the downloads and the regrids — runs anywhere, but
`download_atmosphere` is close to 10000 OPeNDAP reads per simulated year and the Geonorge
bathymetry source is a 4.5 GB FileGDB, so both halves are often better off in the cloud.

This directory runs any pipeline step of any setup on a shared GCE VM, and moves data between your
machine and a Cloud Storage bucket. Nothing about a fjord or a step is baked into the image: which
setup and which steps are arguments to one launch.

## The shape of this deployment

**You cannot create, start, stop or delete instances, and you cannot set instance metadata.** Two
Cloud Run jobs, owned by an admin, create or start one fixed VM on your behalf. After that you have
SSH with sudo, and that is the whole interface. Everything here is built around those facts:

- There is **one** VM, named by `VM_NAME` (`fjordsim-gpu`), not one per run.
- There is **no startup script** — `run` ships a script over `scp` and detaches it with `nohup`.
- Nothing self-deletes. `down` shuts the VM down *from inside* (`sudo poweroff`), which is the only
  stop available without `compute.instances.stop`.
- The VM's zone is **discovered, never configured**: the launcher tries `europe-west4-{a,b,c}` and
  takes whichever has L4 capacity.

Measured access for `shamil.iakubov@niva.no` in `nivatest-1` (`testIamPermissions`, 2026-09-16):

| Have | Scope |
|---|---|
| `run.jobs.run`, `run.jobs.get`, `run.executions.*`, `run.tasks.*` | the two launcher jobs |
| `compute.instances.osAdminLogin` / `.osLogin` — SSH **with sudo** | the instance `fjordsim-gpu` only |
| `compute.instances.get` / `.list` | project-wide |
| `roles/storage.admin` | `gs://fjordsim` |
| `iam.serviceAccounts.actAs` | `fjordsim-sa@nivatest-1.iam.gserviceaccount.com` |
| AR push (`uploadArtifacts`, `tags.create/update`) | `europe-west1-docker.pkg.dev/niva-cd/vertex-images` |

| Do **not** have | Consequence |
|---|---|
| `compute.instances.create/start/stop/delete/setMetadata` | no per-run VMs, no startup script, no self-delete, no real stop |
| `run.jobs.runWithOverrides` | the jobs' retry counts and zone list are fixed |
| `secretmanager.*` | the NVE key comes from your own environment, not Secret Manager |

> **A 403 from `list` proves nothing.** `gcloud storage buckets list` and
> `gcloud artifacts repositories list` are *project-level* calls, and they are commonly denied to
> accounts that nonetheless have full access to a specific bucket and can push to a specific
> repository. Test what you will actually use: `gcloud storage ls "$BUCKET/"` for the bucket, and a
> real push (`gcp/fjordsim-gcp build`) for the registry.

## How it fits together

```
your machine                    gs://BUCKET                     the VM (one, shared, long-lived)
────────────                    ───────────                     ────────────────────────────────
~/FjordSim_data/<fjord>  ──push──►  data/<fjord>/  ──stage in──►  /mnt/stage/data/<fjord>
                                                                        │  bind mount
                                    configs/<name>.jl ──────────►  container: FJORDSIM_DATA_ROOT=/data
                                                                        │
~/FjordSim_results/<fjord> ◄─pull─  results/<fjord>/ ◄─stage out─  /mnt/stage/results/<fjord>
```

The package never learns about Cloud Storage. `FJORDSIM_DATA_ROOT` and `FJORDSIM_RESULTS_ROOT` name
the parent of the per-fjord directories (see `fjord_data_root` in `src/Configs.jl`), the VM
bind-mounts staged directories onto them, and every path in every setup follows. Staging runs on the
*host*, which has `gcloud` installed and the VM's own service account — so the image carries no
cloud SDK and no credentials.

Object storage is not used as a filesystem. `add_rivers` copies a 360 MB file, `Checkpointer` lists
a directory, and the NetCDF writers append — all of which a FUSE mount does badly. Copy in, run,
copy out.

**Nothing on the VM is durable.** The launcher creates it as a Spot instance with
`--instance-termination-action=DELETE` and `--max-run-duration=24h`, so a preemption or 24 hours of
uptime deletes the instance *and its boot disk*. A guest `poweroff` is different — that only stops
it, and the disk survives. This is why the remote run script syncs results to the bucket every ten
minutes, and why `bootstrap` has to run again after a preemption but not after a `down`/`up`.

## One-time setup

### 1. Fill in `gcp/config.env`

```bash
cp gcp/config.env.example gcp/config.env
```

Six values, all of them facts about the deployment:

| Variable | What it is |
|---|---|
| `PROJECT_ID` | the project everything is billed to |
| `BUCKET` | `gs://…` — staged inputs, results and uploaded configs |
| `IMAGE_URI` | the full image URI; everything else about the registry follows from it |
| `JOB_REGION` | where the launcher jobs live — its *only* use is `gcloud run jobs execute` |
| `VM_NAME` | the instance the jobs manage |
| `LAUNCH_JOB_SPOT`, `LAUNCH_JOB_STANDARD` | the jobs themselves |

You cannot list any of the last three, so ask the admin for them.

Access to Artifact Registry is granted per *repository*, so the one you can push to may sit in
another project or region entirely. `IMAGE_URI` is therefore written out in full rather than
assembled from a project and a region, and both `bootstrap` and the remote run script take the
Docker auth host from it — which is all it takes to make the split work. The VM's service account
then needs `roles/artifactregistry.reader` on *that* repository's project.

Deliberately absent: machine type, accelerator, disk and run duration (owned by the launcher jobs —
changing one is an admin request); `ZONE` (discovered from the running instance); `SERVICE_ACCOUNT`
(the job picks the VM's identity, and `bootstrap` reads it from the metadata server when it needs to
name one); and anything about the NVE key (`run` takes it from your environment).

Every variable is also overridable per invocation, so a one-off needs no edit at all:

```bash
BUCKET=gs://other-bucket gcp/fjordsim-gcp run --config oslofjorden --steps run_simulation --gpu
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
gcp/fjordsim-gcp up                                  # launcher job: create or start the VM
gcp/fjordsim-gcp bootstrap                           # driver, Docker, toolkit, image pull
gcp/fjordsim-gcp run --config drammensfjorden --steps run_simulation --gpu
gcp/fjordsim-gcp logs <run-id>
gcp/fjordsim-gcp pull-results drammensfjorden
gcp/fjordsim-gcp down                                # stop billing; the disk survives
```

`up` is idempotent and does the right thing in every state: it creates the VM if it is gone, starts
it if it is stopped, and exits cleanly if it is already running. It can also block for a long
while — the job retries every zone for up to ten rounds with backoff capped at ten minutes when
there is no L4 capacity.

`bootstrap` is idempotent too. On the plain Ubuntu image the launcher currently uses it installs the
NVIDIA driver (~5 min, then a reboot it drives itself), Docker and the container toolkit; on a
second call, or after `down`/`up`, every stage is a no-op. After a *preemption* it starts over,
because the disk is gone.

### Prepare locally, simulate in the cloud

Preparing on a laptop is fine and often simplest; only the simulation really needs the cloud.

```bash
julia --project -m FjordSim prepare_bathymetry --config drammensfjorden
julia --project -m FjordSim prepare_forcing    --config drammensfjorden
julia --project -m FjordSim add_rivers         --config drammensfjorden

gcp/fjordsim-gcp push-data drammensfjorden          # only the top-level *.nc, ~970 MB
gcp/fjordsim-gcp run --config drammensfjorden --steps run_simulation --gpu
gcp/fjordsim-gcp pull-results drammensfjorden
```

`push-data` uploads only the prepared NetCDFs, because that is all a simulation opens. `--all` adds
the raw sources (`norkyst/`, `nora3/`, the FileGDB), which only a `prepare_*` step needs.

### Prepare in the cloud

Better for the slow steps, and it keeps the raw sources off your disk entirely.

```bash
gcp/fjordsim-gcp run --config oslofjorden \
    --steps prepare_bathymetry,download_forcing,prepare_forcing,add_rivers
gcp/fjordsim-gcp run --config oslofjorden \
    --steps download_atmosphere,prepare_atmosphere
gcp/fjordsim-gcp run --config oslofjorden --steps run_simulation --gpu
```

Omitting `--gpu` hides the GPU from the container, and `architecture = :auto` resolves to `CPU()` —
the same config runs either way without an edit. (The VM itself always has an L4; that is the
launcher's choice, not ours.) Anything a prep job writes under the fjord's data directory is synced
back to the bucket when it finishes, so the next job finds it.

`add_rivers` needs `NVE_API_KEY` in your own environment — the same variable a local run uses. `run`
writes it to a mode-600 file on the VM rather than passing it as an argument, which would be
world-readable in `ps`. Free key: https://hydapi.nve.no/Users

`drammensfjorden` reads Oslofjord's FileGDB and NORA3 files by absolute path, so preparing it needs
that fjord staged too:

```bash
gcp/fjordsim-gcp run --config drammensfjorden --steps prepare_bathymetry \
    --stage oslofjorden --stage-all
```

### A one-off experiment

A config file ending in `.jl` is uploaded to the bucket and fetched by the VM, so trying a variant
never means touching the repository. `examples/oslofjorden_npzd.jl` is the worked example. Name the
file after the fjord — the basename is what roots the staged directories.

```bash
gcp/fjordsim-gcp run --config ./oslofjorden_npzd.jl --steps run_simulation --gpu
```

### Resuming after a preemption

The VM is Spot, so a long run will be interrupted. Recovery is `up`, `bootstrap`, then the same run
id with `--resume`:

```bash
gcp/fjordsim-gcp run --config oslofjorden --steps run_simulation --gpu --run-id oslo-2020
# after a preemption:
gcp/fjordsim-gcp up && gcp/fjordsim-gcp bootstrap
gcp/fjordsim-gcp run --config oslofjorden --steps run_simulation --gpu --run-id oslo-2020 --resume
```

`--resume` stages `results/<fjord>/` back down and sets `pickup = true` (`gcp/resume.jl`), so the
run continues from the newest checkpoint rather than restarting. This only works if the setup names
a `CheckpointWriter` — both registered setups do.

Checkpoint filenames carry no run tag and are scoped to `results_root`, so **two live runs must not
share a fjord**. With a single shared VM that is easy to do by accident. `fjordsim-gcp` refuses to
start a run whose id is already marked running; give a different `--run-id` for a genuinely separate
experiment, and a `results_root` to match.

### Watching a run

```bash
gcp/fjordsim-gcp status            # every run marker in the bucket
gcp/fjordsim-gcp logs oslo-2020    # tail on the VM, bucket copy once it is gone
gcp/fjordsim-gcp ssh               # the VM is a normal box; docker ps, nvidia-smi, etc.
```

Results sync up every 10 minutes while the run is in flight, so `pull-results` mid-run gives you
something to plot.

## Cost

- **Nothing stops on its own.** A finished run leaves the VM up and billing; `down` is a manual step
  and the only one available. `--max-run-duration=24h` on the instance is the backstop, and it
  *deletes* rather than stops.
- Spot is roughly a 60-70% discount and is survivable here because of checkpoint resume.
- Keep the bucket near the VM (`europe-west4`) — that traffic repeats all run long, so staging is free.
  The image is the exception; it is pulled cross-region once per bootstrap. Pulling results to your laptop is egress — a few
  hundred MB per run.
- The stage-in filter is the other lever: a simulation pulls ~900 MB instead of ~5.5 GB because it
  never touches the raw sources.

## Asking the admin

Ranked by value per unit of admin effort. Nothing here blocks you today — the tooling works as it
stands — but each one removes a sharp edge.

1. **Use a Deep Learning VM image in both launcher jobs.** One flag pair, no new IAM, and it
   deletes the slowest and most fragile part of the loop: `bootstrap` currently spends ~10 minutes
   installing a driver, Docker and the container toolkit on every VM *creation*, and the VM is
   re-created on every preemption.
   ```diff
   -  --image-family=ubuntu-2204-lts  --image-project=ubuntu-os-cloud
   +  --image-family=common-cu129-ubuntu-2204-nvidia-580  --image-project=deeplearning-platform-release
   ```
   `gcp/bootstrap.sh` stays correct either way — it just becomes a no-op.

2. **`roles/compute.instanceAdmin.v1` on the instance `fjordsim-gpu`** (resource-level, not
   project-level). Combined with the `actAs` already granted, this alone restores `stop`, `delete`
   and `setMetadata` — real cost control instead of a guest `poweroff`. An instance-level binding
   dies with the instance, so the jobs' `grant_shamil_admin()` should add it alongside
   `osAdminLogin`. Narrower alternatives: a custom role with just
   `compute.instances.{start,stop,delete}`, or a third Cloud Run job `fjordsim-stop`.

3. **`roles/artifactregistry.reader` for the VM's service account on the image repository.**
   Needed for `docker pull` on the VM and currently unverified: `iam.serviceAccounts.getAccessToken`
   is not granted, so the account cannot be impersonated to test it. `bootstrap` is what finds out,
   and it prints the exact identity to name — read from the VM's metadata server, so it is the one
   actually in use rather than a guess. Workaround without admin: `docker save | gzip` the image to
   the bucket and `docker load` on the VM.

4. **Give the standard launcher its own VM name** (e.g. `fjordsim-gpu-std`), or delete the current
   Spot `fjordsim-gpu`. The two jobs share one instance name and each refuses to touch an instance
   created with the other provisioning model, so while a Spot `fjordsim-gpu` exists —
   and it cannot be deleted from this account — `fjordsim-launch-standard` can never succeed.

5. **A persistent data disk** attached with `auto-delete=no`, holding `/var/lib/docker` and
   `/mnt/stage`. Turns a preemption into a restart rather than a full re-bootstrap plus a 2.7 GB
   image pull. More admin work than §1; only worth it if preemptions turn out to be frequent.

6. **`run.jobs.runWithOverrides` on the two jobs**, so `MAX_ROUNDS` and the backoff can be set per
   launch instead of being baked in at ten rounds. Minor.

7. **`roles/secretmanager.secretAccessor` on the NVE key secret**, only if passing `NVE_API_KEY`
   from your own environment is unacceptable. The env-var route works today.

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
  with `envsubst` and run on the VM), `Dockerfile` / `build_image.sh` (the image), `preflight.jl`
  (fails a GPU run that silently fell back to CPU) and `resume.jl` (sets `pickup = true`).
