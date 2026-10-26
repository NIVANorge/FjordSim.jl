#!/usr/bin/env bash
#
# Build the FjordSim image and push it to Artifact Registry.
#
# Run from the repository root (the build context is the repo, so the Dockerfile can COPY
# Project.toml, Manifest.toml and src/). Values come from gcp/config.env, overridable per call.
#
# `--environment DIR` builds a second image on top of that one, whose Julia environment is the
# setup environment in DIR (gcp/Dockerfile.environment), and pushes that instead.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(dirname "$here")"

# shellcheck source=/dev/null
source "$here/load_config.sh"
_fjordsim_load_config "$here/config.env"

IMAGE_URI="${IMAGE_URI:?set IMAGE_URI in gcp/config.env}"

environment=""
case "${1:-}" in
    "") ;;
    --environment) environment="${2:?--environment takes a directory}" ;;
    *) echo "build_image: unknown option $1" >&2; exit 1 ;;
esac

# The manifest is git-ignored but must reach the image, so the container resolves to the same
# package versions as this machine. Without it Pkg re-resolves and the image can drift.
[[ -f "$root/Manifest.toml" ]] || {
    echo "build_image: $root/Manifest.toml is missing." >&2
    echo "Run 'julia --project -e \"using Pkg; Pkg.instantiate()\"' first — the image pins what it resolves to." >&2
    exit 1
}

if [[ -z "$environment" ]]; then
    echo "Building ${IMAGE_URI}"
    docker build -f "$here/Dockerfile" -t "${IMAGE_URI}" "$root"
else
    environment="$(cd "$environment" && pwd)"
    [[ -f "$environment/Project.toml" && -f "$environment/Manifest.toml" ]] || {
        echo "build_image: $environment needs a Project.toml and a Manifest.toml." >&2
        echo "Run 'julia --project=$environment -e \"using Pkg; Pkg.instantiate()\"' first." >&2
        exit 1
    }
    repository="$(git -C "$environment" rev-parse --show-toplevel 2>/dev/null || echo "$environment")"

    echo "Building the FjordSim image fjordsim:base"
    docker build -f "$here/Dockerfile" -t fjordsim:base "$root"

    echo "Building ${IMAGE_URI} on it, with the environment in $environment"
    docker build -f "$here/Dockerfile.environment" -t "${IMAGE_URI}" \
        --build-arg BASE_IMAGE=fjordsim:base \
        --build-arg ENVIRONMENT="$(realpath --relative-to="$repository" "$environment")" \
        --build-context environment="$repository" \
        "$here"
fi

echo "Pushing ${IMAGE_URI}"
docker push "${IMAGE_URI}"

echo "${IMAGE_URI}"
