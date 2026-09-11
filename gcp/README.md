# Running FjordSim on GCP

`run_simulation` needs a GPU. Everything else — the downloads and the regrids — runs anywhere, but
`download_atmosphere` is close to 10000 OPeNDAP reads per simulated year and the Geonorge
bathymetry source is a 4.5 GB FileGDB, so both halves are often better off in the cloud.

This directory runs any pipeline step of any setup on a GCE VM, and moves data between your machine
and a Cloud Storage bucket. Nothing about a fjord or a step is baked into the image: which setup,
which steps and what hardware are all arguments to one launch.

## How it fits together

```
your machine                    gs://BUCKET                     a GCE VM (deleted when done)
────────────                    ───────────                     ────────────────────────────
~/FjordSim_data/<fjord>  ──push──►  data/<fjord>/  ──stage in──►  /mnt/stage/data/<fjord>
                                                                        │  bind mount
                                    configs/<name>.jl ──────────►  container: FJORDSIM_DATA_ROOT=/data
                                                                        │
~/FjordSim_results/<fjord> ◄─pull─  results/<fjord>/ ◄─stage out─  /mnt/stage/results/<fjord>
```

The package never learns about Cloud Storage. `FJORDSIM_DATA_ROOT` and `FJORDSIM_RESULTS_ROOT` name
the parent of the per-fjord directories (see `fjord_data_root` in `src/Configs.jl`), the VM
bind-mounts staged directories onto them, and every path in every setup follows. Staging runs on the
*host*, which is a Deep Learning VM with `gcloud` already installed — so the image carries no cloud
SDK and no credentials.

Object storage is not used as a filesystem. `add_rivers` copies a 360 MB file, `Checkpointer` lists
a directory, and the NetCDF writers append — all of which a FUSE mount does badly. Copy in, run,
copy out.

## One-time setup

### 1. Fill in `gcp/config.env`

Everything deployment-specific — which project, which bucket, which registry, what hardware — lives
in that one file. It is git-ignored, because it describes your deployment rather than the package,
and it is the *only* file you edit: the scripts read every value from it and hardcode none of them.

```bash
cp gcp/config.env.example gcp/config.env
```

At minimum, set these three:

| Variable | What it is | Where to get it |
|---|---|---|
| `PROJECT_ID` | the GCP project everything is billed to | `gcloud config get-value project`, or whoever gave you access |
| `BUCKET` | `gs://…` — staged inputs, results and uploaded configs | an existing bucket, or create one in step 3 |
| `SERVICE_ACCOUNT` | the identity the job VMs run as | created in step 3; `fjordsim-runner@PROJECT_ID.iam.gserviceaccount.com` with the naming used here |

Set `REGION` and `ZONE` to wherever the bucket already is, rather than keeping the `europe-west4`
default:

```bash
gcloud storage buckets describe gs://YOUR-BUCKET --format='value(location)'
```

The VMs belong in the bucket's region — that is the traffic that repeats all run long, and
cross-region reads are slow and billed as egress. `AR_REPO` is the Artifact Registry repository to
push the image to, `images` unless your project already uses another name.

The registry is the looser of the two. Access to it is often granted per *repository* rather than
per project, so the one you can push to may sit in another project or region entirely. Set
`IMAGE_URI` outright in that case, instead of letting it be derived from `PROJECT_ID`/`REGION`:

```bash
IMAGE_URI=europe-west1-docker.pkg.dev/other-project/other-repo/fjordsim:latest
```

Only the image is affected, and it is pulled once per VM. The runtime service account then needs
`roles/artifactregistry.reader` on that repository's project rather than on your own.

Every variable is also overridable per invocation, so a one-off needs no edit at all:

```bash
BUCKET=gs://other-bucket gcp/fjordsim-gcp run --config oslofjorden --steps run_simulation --gpu
```

### 2. Roles

**On your own account**

| Role | For |
|---|---|
| `roles/artifactregistry.writer` | pushing the image |
| `roles/storage.objectAdmin` on the bucket | syncing data in and out (`roles/storage.admin` only if you also create the bucket) |
| `roles/compute.instanceAdmin.v1` | creating and deleting job VMs |
| `roles/iam.serviceAccountUser` on the runtime SA | launching a VM that runs as it |

**A runtime service account** the VMs run as:

| Role | For |
|---|---|
| `roles/storage.objectAdmin` on the bucket | staging in and out |
| `roles/artifactregistry.reader` | pulling the image |
| `roles/compute.instanceAdmin.v1` | letting a finished job delete its own VM |
| `roles/secretmanager.secretAccessor` | the NVE key, only if you run `add_rivers` on GCP |

**Also**: GPU quota in your region (`NVIDIA_L4_GPUS` for the default machine shape), and the
`compute`, `artifactregistry`, `storage` and `secretmanager` APIs enabled on the project.

> **A 403 from `list` proves nothing.** `gcloud storage buckets list` and
> `gcloud artifacts repositories list` are *project-level* calls, and they are commonly denied to
> accounts that nonetheless have full access to a specific bucket and can push to a specific
> repository. Test what you will actually use: `gcloud storage ls "$BUCKET/"` for the bucket, and a
> real push (`gcp/fjordsim-gcp build`) for the registry. If you cannot list repositories, ask for
> the repository name instead of guessing it.

### 3. Create whatever does not exist yet

A shared project usually has the bucket and the registry already — run only the pieces you are
missing. These read the file you just filled in, so nothing is retyped:

```bash
source gcp/config.env

gcloud storage buckets create "$BUCKET" --project="$PROJECT_ID" --location="$REGION" \
    --uniform-bucket-level-access

gcloud artifacts repositories create "$AR_REPO" --project="$PROJECT_ID" \
    --repository-format=docker --location="$REGION"

gcloud iam service-accounts create fjordsim-runner --project="$PROJECT_ID"

gcloud storage buckets add-iam-policy-binding "$BUCKET" \
    --member="serviceAccount:$SERVICE_ACCOUNT" --role=roles/storage.objectAdmin
for role in roles/artifactregistry.reader roles/compute.instanceAdmin.v1 roles/secretmanager.secretAccessor; do
    gcloud projects add-iam-policy-binding "$PROJECT_ID" \
        --member="serviceAccount:$SERVICE_ACCOUNT" --role="$role"
done

# Only needed if you run `add_rivers` on GCP. Free key: https://hydapi.nve.no/Users
printf '%s' "$NVE_API_KEY" | gcloud secrets create "$NVE_API_KEY_SECRET" \
    --project="$PROJECT_ID" --data-file=- --replication-policy=automatic
```

### 4. Build and push the image

```bash
source gcp/config.env
gcloud auth configure-docker "${REGION}-docker.pkg.dev"
gcp/fjordsim-gcp build
```

The image bakes the precompile cache for the whole Oceananigans / NumericalEarth / CairoMakie
stack. That build is slow once so that no VM pays 10-20 minutes of precompilation on boot, while
being billed for a GPU that is doing nothing.

## Workflows

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

No `--gpu` means no accelerator, and `architecture = :auto` resolves to `CPU()` — the same config
runs on both without an edit. Anything a prep job writes under the fjord's data directory is synced
back to the bucket when it finishes, so the next job finds it.

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

### Long runs on Spot

```bash
gcp/fjordsim-gcp run --config oslofjorden --steps run_simulation --gpu --spot --run-id oslo-2020
# after a preemption, same run id:
gcp/fjordsim-gcp run --config oslofjorden --steps run_simulation --gpu --spot --run-id oslo-2020 --resume
```

`--resume` stages `results/<fjord>/` back down and sets `pickup = true` (`gcp/resume.jl`), so the
run continues from the newest checkpoint rather than restarting. This only works if the setup names
a `CheckpointWriter` — both registered setups do.

Checkpoint filenames carry no run tag and are scoped to `results_root`, so **two live runs must not
share a fjord**. `fjordsim-gcp` refuses to start a run whose id is already marked running; give a
different `--run-id` for a genuinely separate experiment, and a `results_root` to match.

### Watching a run

```bash
gcp/fjordsim-gcp status            # every run marker in the bucket
gcp/fjordsim-gcp logs oslo-2020    # serial console while up, bucket copy afterwards
gcp/fjordsim-gcp ssh oslo-2020     # the VM is a normal box; docker ps, nvidia-smi, etc.
```

Results sync up every 10 minutes while the run is in flight, so `pull-results` mid-run gives you
something to plot.

## Cost

- VMs delete themselves when the last step finishes. `--keep` leaves one up for debugging; you then
  delete it yourself. `MAX_RUN_DURATION` (24h by default) is the backstop if the self-delete never
  runs.
- `--spot` is roughly a 60-70% discount and is safe here because of checkpoint resume.
- Bucket, registry and VMs share a region, so staging is free. Pulling results to your laptop is
  egress — a few hundred MB per run.
- The stage-in filter is the other lever: a GPU VM pulls ~900 MB instead of ~5.5 GB because it never
  touches the raw sources.

## Reproducibility

`Manifest.toml` is git-ignored but *is* copied into the image, so the container resolves to exactly
the package versions on the machine that built it. The trade-off is that an image cannot be rebuilt
from a clean clone. If you want CI builds later, drop `/Manifest*.toml` from `.gitignore` and commit
it — nothing else has to change.

## Notes

- The `Scratch.jl` bathymetry cache (`src/Bathymetry/geonorge.jl`) lives in the Julia depot inside
  the image, so it is rebuilt on each VM. `raw_directory` is an ordinary config field; point it
  under `data_root` in a config file if you would rather it were staged and reused.
- `--dry-run` prints the `gcloud` command and the fully rendered startup script without creating
  anything. It is the fastest way to see what a launch will actually do.
