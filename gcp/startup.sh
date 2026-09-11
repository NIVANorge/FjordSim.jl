#!/usr/bin/env bash
#
# GCE startup script, rendered by gcp/fjordsim-gcp and run once on VM boot as root.
#
# Runs on the host, not in the container: a Deep Learning VM already has gcloud, Docker and the
# NVIDIA runtime, so staging happens here and the image stays free of any cloud SDK or credentials.
# The container sees two bind-mounted directories and nothing else.
#
# Placeholders below are filled by `envsubst` with an explicit allowlist, so ordinary shell
# variables in this file survive rendering untouched.

set -euo pipefail
exec > >(tee -a /var/log/fjordsim-startup.log) 2>&1

readonly BUCKET="${BUCKET}"
readonly FJORD="${FJORD}"
readonly RUN_ID="${RUN_ID}"
readonly STEPS="${STEPS}"
readonly CONFIG_REF="${CONFIG_REF}"
readonly CONFIG_OBJECT="${CONFIG_OBJECT}"
readonly IMAGE_URI="${IMAGE_URI}"
readonly REGION="${REGION}"
readonly ZONE="${ZONE}"
readonly USE_GPU="${USE_GPU}"
readonly RESUME="${RESUME}"
readonly KEEP_VM="${KEEP_VM}"
readonly STAGE_ALL="${STAGE_ALL}"
readonly EXTRA_FJORDS="${EXTRA_FJORDS}"
readonly NVE_API_KEY_SECRET="${NVE_API_KEY_SECRET}"

# The metadata server, not `hostname`: some images return an FQDN there, and the self-delete below
# is the only thing that stops a finished job from billing until MAX_RUN_DURATION expires.
readonly INSTANCE_NAME="$(curl -sf -H 'Metadata-Flavor: Google' \
    http://metadata.google.internal/computeMetadata/v1/instance/name)"

readonly STAGE=/mnt/stage
readonly DATA="$STAGE/data"
readonly RESULTS="$STAGE/results"
readonly RUN_PREFIX="$BUCKET/runs/$RUN_ID"

mark() { echo "$1" | gcloud storage cp - "$RUN_PREFIX/status" --quiet || true; }

finish() {
    local code=$?
    [[ -n "${SYNC_PID:-}" ]] && kill "$SYNC_PID" 2>/dev/null || true

    echo "=== final stage-out ==="
    sync_results || true
    # Newly prepared inputs are a product too: a prep job exists to put them in the bucket.
    gcloud storage rsync -r "$DATA/$FJORD" "$BUCKET/data/$FJORD" --quiet || true
    gcloud storage cp /var/log/fjordsim-startup.log "$RUN_PREFIX/startup.log" --quiet || true

    if [[ $code -eq 0 ]]; then mark succeeded; else mark "failed(exit=$code)"; fi

    if [[ "$KEEP_VM" == "true" ]]; then
        echo "--keep was given; leaving $INSTANCE_NAME running. Delete it when done:"
        echo "  gcloud compute instances delete $INSTANCE_NAME --zone=$ZONE"
    else
        echo "=== deleting $INSTANCE_NAME ==="
        gcloud compute instances delete "$INSTANCE_NAME" --zone="$ZONE" --quiet
    fi
}
trap finish EXIT

sync_results() {
    [[ -d "$RESULTS/$FJORD" ]] || return 0
    gcloud storage rsync -r "$RESULTS/$FJORD" "$BUCKET/results/$FJORD" --quiet
}

# A bucket prefix that does not exist yet is not an error: the first prep job for a fjord starts
# with nothing staged, and `rsync` from a missing source would otherwise abort the run.
stage_in() {
    local remote="$1" local_dir="$2"
    shift 2
    mkdir -p "$local_dir"
    if gcloud storage ls "$remote" >/dev/null 2>&1; then
        echo "staging $remote -> $local_dir"
        gcloud storage rsync -r "$@" "$remote" "$local_dir" --quiet
    else
        echo "nothing at $remote yet; starting empty"
    fi
}

mark running
mkdir -p "$DATA" "$RESULTS" "$STAGE/configs"

echo "=== image ==="
gcloud auth configure-docker "${REGION}-docker.pkg.dev" --quiet
docker pull "$IMAGE_URI"

echo "=== stage in ==="
# A simulation reads only the prepared NetCDFs, which all sit at the top level of the fjord's data
# directory; everything in a subdirectory is raw source that only a prepare_* step opens. Excluding
# subdirectories keeps a GPU VM from pulling the 4.5 GB Geonorge FileGDB it will never read.
if [[ "$STAGE_ALL" == "true" ]]; then
    stage_in "$BUCKET/data/$FJORD" "$DATA/$FJORD"
else
    stage_in "$BUCKET/data/$FJORD" "$DATA/$FJORD" -x '.*/.*'
fi

# A setup may read another fjord's raw downloads by absolute path — drammensfjorden shares
# Oslofjord's FileGDB and NORA3 files — so a prep job for it must stage that fjord too.
# Iterating over a lowercase alias rather than the rendered placeholder: envsubst substitutes in
# place, and an empty value would otherwise leave the syntax error `for extra in ; do`.
extra_fjords="$EXTRA_FJORDS"
for extra in $extra_fjords; do
    stage_in "$BUCKET/data/$extra" "$DATA/$extra"
done

# Results come back down whenever there is a run to continue, so `pickup` finds its checkpoint.
if [[ "$RESUME" == "true" ]]; then
    stage_in "$BUCKET/results/$FJORD" "$RESULTS/$FJORD"
fi

# An out-of-tree config travels through the bucket rather than being baked into the image, which
# is what makes a one-off experiment a local file instead of a repo change.
config_arg="$CONFIG_REF"
if [[ -n "$CONFIG_OBJECT" ]]; then
    gcloud storage cp "$CONFIG_OBJECT" "$STAGE/configs/$(basename "$CONFIG_OBJECT")" --quiet
    config_arg="/configs/$(basename "$CONFIG_OBJECT")"
fi

docker_args=(
    --rm
    -v "$DATA:/data"
    -v "$RESULTS:/results"
    -v "$STAGE/configs:/configs"
    -e FJORDSIM_DATA_ROOT=/data
    -e FJORDSIM_RESULTS_ROOT=/results
)
[[ "$USE_GPU" == "true" ]] && docker_args+=(--gpus all)

# Only `add_rivers` reads it, so a run that never calls that step needs no secret access at all.
if [[ "$STEPS" == *add_rivers* ]]; then
    echo "=== fetching NVE_API_KEY ==="
    nve_key="$(gcloud secrets versions access latest --secret="$NVE_API_KEY_SECRET")"
    docker_args+=(-e "NVE_API_KEY=$nve_key")
fi

if [[ "$USE_GPU" == "true" ]]; then
    echo "=== gpu preflight ==="
    nvidia-smi
    docker run "${docker_args[@]}" --entrypoint julia "$IMAGE_URI" \
        --project /workspace/gcp/preflight.jl
fi

# Sync while the run is in flight, so progress is visible from a laptop and a preempted Spot VM
# loses at most one interval rather than the whole run.
( while true; do sleep 600; sync_results || true; done ) &
SYNC_PID=$!

for step in ${STEPS//,/ }; do
    echo "=== $step --config $config_arg ==="
    if [[ "$step" == "run_simulation" && "$RESUME" == "true" ]]; then
        docker run "${docker_args[@]}" --entrypoint julia "$IMAGE_URI" \
            --project /workspace/gcp/resume.jl "$config_arg"
    else
        docker run "${docker_args[@]}" "$IMAGE_URI" "$step" --config "$config_arg"
    fi
done

echo "=== all steps finished ==="
