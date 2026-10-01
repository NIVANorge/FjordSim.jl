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
using FjordSim.CLI: log_path, show_compact_error, tee_output

length(ARGS) == 1 || error("usage: julia --project gcp/resume.jl <setup-name-or-config-path>")

config = fjord_config(ARGS[1])
simulation_config = config.simulation_config

isnothing(simulation_config) && error(
    "$(ARGS[1]) names no simulation config, so there is no run to resume.",
)

simulation_config.pickup = true

# Lets a resume be forced onto the CPU on a machine whose GPU would otherwise be picked (`:auto`
# in every setup so far), to get a fully-resolved Julia stacktrace out of a bug that only throws a
# GPU-side exception whose on-device backtrace can't be symbolicated.
get(ENV, "FJORDSIM_CPU", "false") == "true" && (simulation_config.architecture = :cpu)

# Bare `run_simulation(config)` here would leave an uncaught exception to Julia's own top-level
# handler, which prints the full, unabbreviated stacktrace (every type parameter spelled out) and
# skips `show_compact_error`'s timeout guard against a GPU exception hanging while it prints — see
# `run_step` (src/CLI.jl), which this mirrors so a resumed run fails exactly as readably as a fresh one.
#
# Also mirrors `main`'s `tee_output` wrapping (src/CLI.jl), so a resume leaves the same
# `fjordsim_<run_tag>.log` transcript beside its output that a fresh run does. `exit` runs after
# `tee_output` returns, not inside its `do` block, since `exit` skips `tee_output`'s `finally` and
# would leave the log file unflushed and unclosed.
log_file = log_path(simulation_config)
@info "Logging to $(abspath(log_file))"

exit_code = tee_output(log_file) do
    try
        run_simulation(config)
        0
    catch exception
        println(stderr, "fjordsim: resume failed on $(ARGS[1]):")
        show_compact_error(stderr, exception, catch_backtrace())
        1
    end
end

exit(exit_code)
