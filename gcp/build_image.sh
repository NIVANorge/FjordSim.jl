#!/usr/bin/env bash
#
# Build the FjordSim image and push it to Artifact Registry.
#
# Run from the repository root (the build context is the repo, so the Dockerfile can COPY
# Project.toml, Manifest.toml and src/). Values come from gcp/config.env, overridable per call.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(dirname "$here")"

# shellcheck source=/dev/null
source "$here/load_config.sh"
_fjordsim_load_config "$here/config.env"

PROJECT_ID="${PROJECT_ID:?set PROJECT_ID in gcp/config.env}"
REGION="${REGION:-europe-west4}"
AR_REPO="${AR_REPO:-images}"
IMAGE="${IMAGE:-fjordsim}"
IMAGE_TAG="${IMAGE_TAG:-latest}"
IMAGE_URI="${IMAGE_URI:-${REGION}-docker.pkg.dev/${PROJECT_ID}/${AR_REPO}/${IMAGE}:${IMAGE_TAG}}"

# The manifest is git-ignored but must reach the image, so the container resolves to the same
# package versions as this machine. Without it Pkg re-resolves and the image can drift.
[[ -f "$root/Manifest.toml" ]] || {
    echo "build_image: $root/Manifest.toml is missing." >&2
    echo "Run 'julia --project -e \"using Pkg; Pkg.instantiate()\"' first — the image pins what it resolves to." >&2
    exit 1
}

echo "Building ${IMAGE_URI}"
docker build -f "$here/Dockerfile" -t "${IMAGE_URI}" "$root"

echo "Pushing ${IMAGE_URI}"
docker push "${IMAGE_URI}"

echo "${IMAGE_URI}"
