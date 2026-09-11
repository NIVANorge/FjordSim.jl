# Fail fast when the GPU a job asked for is not actually usable inside the container.
#
# Without this, `architecture = :auto` silently resolves to `CPU()` (see
# `interpolation_architecture(::Val{:auto})` in src/Forcing/Forcing.jl) and the run proceeds —
# billing a GPU VM for a simulation that is hundreds of times too slow to finish, and looking
# healthy the whole time. The usual cause is the host NVIDIA driver not being injected into the
# container, which `CUDA.functional()` reports and `nvidia-smi` on the host does not.
using CUDA

CUDA.versioninfo()

if !CUDA.functional()
    println(stderr, """
        preflight: no usable GPU inside the container.

        `CUDA.functional()` is false, so a config with `architecture = :auto` would fall back to
        the CPU and a config with `:gpu` would error later. Check that the VM was created with an
        accelerator, that the NVIDIA driver is installed on the host, and that `docker run` was
        given `--gpus all`.
        """)
    exit(1)
end

println("preflight: GPU OK — ", CUDA.name(CUDA.device()))
