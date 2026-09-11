# FjordSim on a GPU VM. Built and pushed by gcp/build_image.sh, run by gcp/startup.sh.
#
# The image holds Julia, the instantiated package environment and a warm precompile cache — and
# nothing else. Staging to and from Cloud Storage happens on the *host*, which is a Deep Learning
# VM with gcloud already installed, so this image needs no cloud SDK and no credentials: it sees
# only the two directories the host bind-mounts onto FJORDSIM_DATA_ROOT and FJORDSIM_RESULTS_ROOT.
FROM nvidia/cuda:12.4.1-runtime-ubuntu22.04

ARG JULIA_VERSION=1.12.7

ENV DEBIAN_FRONTEND=noninteractive \
    JULIA_DEPOT_PATH=/opt/julia \
    JULIA_NUM_THREADS=auto \
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
RUN julia --project -e 'using Pkg; \
        Pkg.instantiate(); \
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
