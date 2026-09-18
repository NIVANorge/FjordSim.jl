#!/usr/bin/env bash
#
# Rendered by `fjordsim-gcp run` (envsubst), scp'd to the VM and launched detached with nohup. It
# stages inputs down from the bucket, runs each step in the container, syncs results back, and
# leaves a status marker and a copy of its own log in the bucket.
#
# It does not delete the VM at the end: the account has no permission to, and the instance is a
# long-lived singleton shared by every run. Use `fjordsim-gcp down` when you are finished with it.
set -euo pipefail

readonly BUCKET="${BUCKET}"
readonly FJORD="${FJORD}"
readonly RUN_ID="${RUN_ID}"
readonly STEPS="${STEPS}"
readonly CONFIG_REF="${CONFIG_REF}"
readonly CONFIG_OBJECT="${CONFIG_OBJECT}"
readonly IMAGE_URI="${IMAGE_URI}"
readonly USE_GPU="${USE_GPU}"
readonly RESUME="${RESUME}"
readonly STAGE_ALL="${STAGE_ALL}"
readonly EXTRA_FJORDS="${EXTRA_FJORDS}"
readonly DEBUG="${DEBUG}"
readonly CPU="${CPU}"

readonly STAGE=/mnt/stage
readonly DATA="$STAGE/data"
readonly RESULTS="$STAGE/results"
readonly RUN_PREFIX="$BUCKET/runs/$RUN_ID"
readonly RUN_LOG="$HOME/fjordsim-$RUN_ID.log"
readonly NVE_KEY_FILE="$HOME/.fjordsim-nve-key"

# Docker auth is against the host in IMAGE_URI, which need not be in the VM's own region — access
# to Artifact Registry is granted per repository, so the image often lives in another project.
readonly REGISTRY_HOST="${IMAGE_URI%%/*}"

mark() { echo "$1" | gcloud storage cp - "$RUN_PREFIX/status" --quiet || true; }

sync_results() {
    [[ -d "$RESULTS/$FJORD" ]] || return 0
    gcloud storage rsync -r "$RESULTS/$FJORD" "$BUCKET/results/$FJORD" --quiet
}

finish() {
    local code=$?
    [[ -n "${SYNC_PID:-}" ]] && kill "$SYNC_PID" 2>/dev/null || true

    echo "=== final stage-out ==="
    sync_results || true
    # Newly prepared inputs are a product too: a prep job exists to put them in the bucket.
    gcloud storage rsync -r "$DATA/$FJORD" "$BUCKET/data/$FJORD" --quiet || true
    # The VM's disk does not survive a preemption, so this copy is what `logs` falls back to.
    gcloud storage cp "$RUN_LOG" "$RUN_PREFIX/run.log" --quiet || true

    if [[ $code -eq 0 ]]; then mark succeeded; else mark "failed(exit=$code)"; fi

    echo "=== $RUN_ID finished (exit=$code) on $(hostname) ==="
    echo "the VM is left running; stop it with: gcp/fjordsim-gcp down"
}
trap finish EXIT

stage_in() {
    local remote="$1" local_dir="$2"
    shift 2
    mkdir -p "$local_dir"
    # `ls` fails both when the prefix does not exist yet (fine, a first prep run) and when the VM's
    # service account cannot read the bucket (fatal). Only the first message means "empty", so match
    # it and let everything else abort: swallowing a 403 here hands the container an empty /data and
    # the run dies ~12 minutes later complaining about a missing bathymetry file instead.
    local listing
    if listing="$(gcloud storage ls "$remote" 2>&1)"; then
        echo "staging $remote -> $local_dir"
        gcloud storage rsync -r "$@" "$remote" "$local_dir" --quiet
    elif [[ "$listing" == *"matched no objects"* ]]; then
        echo "nothing at $remote yet; starting empty"
    else
        echo "$listing" >&2
        echo "stage_in: cannot list $remote" >&2
        exit 1
    fi
}

mark running
mkdir -p "$DATA" "$RESULTS" "$STAGE/configs"

echo "=== image ==="
gcloud auth configure-docker "$REGISTRY_HOST" --quiet
docker pull "$IMAGE_URI"

echo "=== stage in ==="
if [[ "$STAGE_ALL" == "true" ]]; then
    stage_in "$BUCKET/data/$FJORD" "$DATA/$FJORD"
else
    stage_in "$BUCKET/data/$FJORD" "$DATA/$FJORD" -x '.*/.*'
fi

extra_fjords="$EXTRA_FJORDS"
for extra in $extra_fjords; do
    stage_in "$BUCKET/data/$extra" "$DATA/$extra"
done

if [[ "$RESUME" == "true" ]]; then
    stage_in "$BUCKET/results/$FJORD" "$RESULTS/$FJORD"
fi

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
if [[ "$DEBUG" == "true" ]]; then
    docker_args+=(-e CUDA_LAUNCH_BLOCKING=1)

    # `--check-bounds=yes` and `-g2` (below) both invalidate the image's precompiled pkgimages —
    # neither matches the flags the image was built with — so a debug run otherwise recompiles the
    # whole dependency tree from scratch every time. A persistent depot directory ahead of the
    # image's own in JULIA_DEPOT_PATH lets Julia cache that recompile once and reuse it across
    # debug runs, instead of paying it again on every `--debug` invocation.
    mkdir -p "$STAGE/depot-cache"
    docker_args+=(
        -v "$STAGE/depot-cache:/mnt/depot-cache"
        -e JULIA_DEPOT_PATH=/mnt/depot-cache:/opt/julia
    )
fi
[[ "$CPU" == "true" ]] && docker_args+=(-e FJORDSIM_CPU=true)

# The key is written to a mode-600 file by `fjordsim-gcp run` rather than passed on a command line,
# which would be world-readable in `ps`.
if [[ "$STEPS" == *add_rivers* ]]; then
    echo "=== reading NVE_API_KEY ==="
    [[ -f "$NVE_KEY_FILE" ]] || { echo "$NVE_KEY_FILE is missing" >&2; exit 1; }
    docker_args+=(-e "NVE_API_KEY=$(cat "$NVE_KEY_FILE")")
fi

if [[ "$USE_GPU" == "true" ]]; then
    echo "=== gpu preflight ==="
    nvidia-smi
    docker run "${docker_args[@]}" --entrypoint julia "$IMAGE_URI" \
        --project /workspace/gcp/preflight.jl
fi

( while true; do sleep 600; sync_results || true; done ) &
SYNC_PID=$!

for step in ${STEPS//,/ }; do
    echo "=== $step --config $config_arg ==="
    if [[ "$step" == "run_simulation" && "$RESUME" == "true" ]]; then
        julia_args=(--project)
        if [[ "$DEBUG" == "true" ]]; then
            julia_args+=(--check-bounds=yes)
            # -g2 only helps symbolicate a GPU on-device backtrace; on the CPU, stacktraces are
            # always fully resolved without it, and it has been observed to blow up compile-time
            # memory use without bound on this heavily generic, deeply inlined kernel code — an
            # OOM on a 31 GB VM, scaling with whatever RAM is available rather than the model's
            # actual size.
            [[ "$CPU" == "true" ]] || julia_args+=(-g2)
        fi
        docker run "${docker_args[@]}" --entrypoint julia "$IMAGE_URI" \
            "${julia_args[@]}" /workspace/gcp/resume.jl "$config_arg"
    else
        docker run "${docker_args[@]}" "$IMAGE_URI" "$step" --config "$config_arg"
    fi
done

echo "=== all steps finished ==="
