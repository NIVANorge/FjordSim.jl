# Resume a simulation from its newest checkpoint: `julia --project gcp/resume.jl <setup-or-path>`.
#
# `pickup` is a `SimulationConfig` field rather than a command-line option, so resuming after a
# Spot preemption would otherwise mean editing the setup. `FjordConfig` and `SimulationConfig` are
# both mutable, so flipping it here needs no change to the package and leaves the setup honest
# about how it starts by default.
#
# `run_simulation` resumes into the highest checkpointed loop it finds under `results_root`
# (`resume_loop`, src/Simulations.jl), so the caller's only job is to have staged that directory
# back down before this runs.
using FjordSim

length(ARGS) == 1 || error("usage: julia --project gcp/resume.jl <setup-name-or-config-path>")

config = fjord_config(ARGS[1])
simulation_config = config.simulation_config

isnothing(simulation_config) && error(
    "$(ARGS[1]) names no simulation config, so there is no run to resume.",
)

simulation_config.pickup = true
run_simulation(config)
