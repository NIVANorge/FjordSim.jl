# FjordSim Copilot Instructions

These instructions apply to the whole repository unless a prompt explicitly says otherwise.

## Core Behavior

- Make minimal, surgical edits that map directly to the user request.
- Do not refactor unrelated code, rename symbols, or reformat untouched sections.
- If requirements are ambiguous, ask a clarification question before implementing.
- Prefer the simplest working solution; do not add speculative features or abstractions.

## Project Structure and Style

- Keep one Julia module per file in `src/` and maintain explicit `export` lists.
- Follow existing naming style: functions in `snake_case`; modules/types in `PascalCase`.
- Name module files to match their module/type names in `PascalCase` (for example, `Grid.jl` for `module Grid`).
- Prefer keyword arguments for user-facing configurability; keep required core inputs explicit.
- Match local file style for imports (`using` vs `import`) and only `import` functions you extend.

## FjordSim Domain Conventions

- Preserve and extend existing adapter patterns for Oceananigans/NumericalEarth compatibility.
- When adding backend-like types, implement architecture adaptation patterns already used in the repo.
- Use `NamedTuple`-based composition for model configuration and boundary-condition assembly.
- Keep dataset metadata handling dispatch-based and consistent with existing dataset adapters.

## Documentation and Comments

- Add docstrings to new public functions using Julia triple-quote style.
- Keep comments concise and focused on non-obvious domain logic.
- Avoid noisy comments that restate obvious code behavior.

## Validation

- Validate behavior with tests when practical (best-effort based on runtime cost).
- Primary test command: `julia --project -e "using Pkg; Pkg.test()"`.
- For simulation changes, mention whether `examples/oslofjord.jl` should be run and why.

## Safety for Changes

- Hard rule: do not remove or weaken existing ClimaOcean/NumericalEarth compatibility shims unless explicitly requested.
- Do not delete existing code only because it appears unused unless asked.
- If proposing a broader cleanup, present it separately from the requested fix.
