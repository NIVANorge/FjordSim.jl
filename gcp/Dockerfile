# FjordSim on a GPU VM. Built and pushed by gcp/build_image.sh, run by gcp/remote-run.sh.
#
# The image holds Julia, the instantiated package environment and a warm precompile cache — and
# nothing else. Staging to and from Cloud Storage happens on the *host*, which has gcloud and the
# VM's own service account, so this image needs no cloud SDK and no credentials: it sees only the
# two directories the host bind-mounts onto FJORDSIM_DATA_ROOT and FJORDSIM_RESULTS_ROOT.
FROM nvidia/cuda:12.4.1-runtime-ubuntu22.04

ARG JULIA_VERSION=1.12.7

# JULIA_CPU_TARGET is what makes the precompile caches below portable, and it has to be set at
# run time too — Julia checks a cache file's targets against the *current* setting — which is why
# it lives here rather than on the precompiling RUNs. Julia's default is `native`, so an image
# built on a laptop bakes in that laptop's microarchitecture and every pkgimage is then rejected on
# the VM with "Unable to find compatible target in cached code image", silently recompiling the
# whole stack on billed GPU time. x86-64-v3 (haswell baseline) loads on any modern x86 host; the
# AVX-512 it gives up costs little here, since the host only interpolates forcing and does IO.
ENV DEBIAN_FRONTEND=noninteractive \
    JULIA_DEPOT_PATH=/opt/julia \
    JULIA_NUM_THREADS=auto \
    JULIA_CPU_TARGET=x86-64-v3 \
    PATH=/opt/julia-${JULIA_VERSION}/bin:$PATH

# `git` is here because Pkg shells out to it; the rest is what the JLL-provided native libraries
# (GDAL, NetCDF, Cairo) link against once unpacked.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl git libgomp1 \
 && rm -rf /var/lib/apt/lists/*

RUN curl -fsSL "https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_VERSION%.*}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz" \
      | tar -xz -C /opt \
 && julia --version

WORKDIR /workspace

# Dependencies first, so editing `src/` does not re-resolve, re-download or re-precompile 500-odd
# packages. Manifest.toml is copied deliberately — see "Reproducibility" in gcp/README.md.
#
# Only the *dependencies* are precompiled here. `--project` activates FjordSim's own package
# environment, so a bare `Pkg.precompile()` would try to precompile FjordSim too and fail with
# "Missing source file" — `src/` does not exist until the next layer, which is the whole point of
# splitting them. `Pkg.project().dependencies` is the `[deps]` table, which never lists the
# package itself.
COPY Project.toml Manifest.toml ./

# CUDA_Runtime_jll tags its platform with a "cuda" version that Pkg matches artifacts against, and
# it works that version out by inspecting the NVIDIA driver. There is no driver on the build
# machine, so the tag becomes the literal "none" (augment_platform! in the JLL's
# .pkg/platform_augmentation.jl) and no runtime artifact is ever selected. On the VM CUDA.jl then
# dies with "CUDA runtime not found" in the GPU preflight, before the simulation starts. Pinning the
# version short-circuits the inspection — the JLL's own words: "with a version override, we don't
# need to (and shouldn't) inspect the driver: this supports precompiling with a fixed version
# without a driver available."
#
# This has to sit between instantiate and precompile. The tag drives artifact *selection*, so the
# preference must exist before the second instantiate resolves artifacts, and it is a compile-time
# preference, so it must be in place before anything is precompiled against it. Setting it costs one
# throwaway precompile of the CUDA chain, which is why it happens here rather than after the full
# dependency precompile below.
#
# 12.4 matches this image's CUDA_VERSION and is one of the toolkits the JLL ships; the L4's driver
# (615.71.09, CUDA UMD 13.4) is newer, which is the only compatibility requirement. Written through
# `set_runtime_version!` rather than by hand: a hand-rolled LocalPreferences.toml with the same
# `[CUDA_Runtime_jll] version` entry was silently ignored by `Base.get_preferences`.
RUN julia --project -e 'using Pkg; Pkg.instantiate()' \
 && julia --project -e 'using CUDA; CUDA.set_runtime_version!(v"12.4")' \
 && julia --project -e 'using Pkg; Pkg.instantiate()' \
 && { find /opt/julia/artifacts -name 'libcudart.so*' -print -quit | grep -q . \
      || { echo "no CUDA runtime artifact in the depot: the version preference did not take, and" \
                "this image would fail its GPU preflight on the VM" >&2; exit 1; }; }

RUN julia --project -e 'using Pkg; \
        Pkg.precompile(collect(keys(Pkg.project().dependencies)))'

COPY src ./src
COPY examples ./examples
COPY gcp/preflight.jl gcp/resume.jl ./gcp/

# Precompile FjordSim itself into the image. A cold precompile of the Oceananigans /
# NumericalEarth / CairoMakie stack is 10-20 minutes, and without this every VM would pay it
# again on boot, before the GPU it is being billed for does any work.
RUN julia --project -e 'using FjordSim'

# `-m FjordSim` is the package entry point, so `docker run <image> run_simulation --config X`
# reads exactly like the documented local command.
ENTRYPOINT ["julia", "--project", "-m", "FjordSim"]
CMD ["--help"]
