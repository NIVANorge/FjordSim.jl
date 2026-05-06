---
description: "Use when modifying, refactoring, or extending any exported symbol in FjordSim — functions, types, or overloaded methods. Covers public API stability: argument lists, return types, and export list invariants."
---
# Backward Compatibility — FjordSim Public API

FjordSim users interact with the library through its exported symbols. Changing those interfaces breaks user code silently. Treat every exported symbol as a public contract.

## The Contract

For every exported function, type, or overloaded method:
- **Positional arguments** must not be removed, reordered, or have their types narrowed.
- **Keyword arguments** must not be removed or renamed. New optional keyword arguments may be added with a default value.
- **Return type** (or return structure) must remain the same. If a function returns a named tuple, a struct, or a specific Oceananigans type, it must keep returning that same shape.
- **Overloaded methods** (e.g., a FjordSim method on `ImmersedBoundaryGrid`) must keep their full method signature and return the same object type.

You may freely change internals — algorithms, helper calls, data structures — as long as inputs and outputs are preserved.

## API Catalog — The Source of Truth

**See [api-catalog.md](./api-catalog.md) for the definitive record of all exported symbols, their signatures, and return types.**

The catalog is organized by module and updated whenever:
- A new function or type is added to an `export` list
- An exported function's signature changes
- A struct's fields are modified
- An export is removed (with deprecation details)

**Before making changes to exported code, consult the catalog.** After making changes, update it immediately. This prevents silent API drift.

## Exported Symbols Snapshot

For quick reference, here are the main public entry points. For full details, **see the catalog**:

```
# FjordSim main module
ImmersedBoundaryGrid              # overloaded Oceananigans constructor
forcing_from_file                 # forcing from NetCDF time series
top_bottom_boundary_conditions    # BC tuple factory
coupled_hydrostatic_simulation    # main simulation setup
recursive_merge                   # utility for dict merging
progress                          # logging callback
cell_advection_timescale_coupled_model
NORA3PrescribedAtmosphere         # atmospheric boundary layer
NORA3PrescribedRadiation          # downwelling SW/LW radiation from NORA3
MultiYearNORA3                    # NORA3 metadata container

# Submodules (often re-exported):
# Utils: compute_faces, safe_execute, extract_z_faces, netcdf_to_jld2, save_fts
# FDatasets: DSForcing, DSResults, last_date
# BoundaryConditions: top_bottom_boundary_conditions
# Grids: ImmersedBoundaryGrid
# Atmospheres: NORA3PrescribedAtmosphere, NORA3PrescribedRadiation, MultiYearNORA3
```

## Hard Rules

- **Do not remove** any symbol from an `export` list without explicit instruction.
- **Do not rename** exported symbols. If a better name is needed, keep the old name as an alias and export both.
- **Do not add required positional arguments** to an existing exported function. New required inputs must come from keyword arguments with sensible defaults, or a new function.
- **Do not change the return type** of an exported function. If the internal representation changes, add a conversion to preserve the output shape.
- **Do not narrow dispatch** on an exported overload. If an exported method currently accepts `AbstractGrid`, do not tighten it to `RectilinearGrid` unless explicitly asked.

## When You Must Change a Signature

If a change that breaks the above rules is genuinely necessary:
1. **Say so explicitly** before implementing — describe what breaks and why it is unavoidable.
2. Propose a deprecation path: keep the old signature working (e.g., via a wrapper or default argument) until the user confirms removal is safe.
3. Never silently drop or reorder arguments.

## Adding New Functionality

- New exported functions or types are always safe to add.
- New **optional** keyword arguments (with defaults) on existing functions are safe.
- New methods for existing generic functions are safe as long as they do not shadow the existing exported method.

## Workflow: Adding a New Exported Function/Type

1. **Implement the new function/type** with full docstring and examples.
2. **Add to `export` list** in the module file.
3. **Update the API catalog** ([api-catalog.md](./api-catalog.md)):
   - Add a row to the appropriate module table
   - Include: symbol name, type (Function/Struct/Type), full signature, and purpose
   - For struct returns, document the field names and types
4. **Write tests** covering the new export's behavior.
5. **Verify no accidental breaking changes** by checking that existing exported symbols still have the same signatures.

## Workflow: Modifying an Existing Exported Function

1. **Check the catalog** ([api-catalog.md](./api-catalog.md)) for the current signature.
2. **If your change is safe** (only adding optional kwargs or refining internals):
   - Implement the change
   - Update docstring if needed
   - Update the catalog entry if the signature changed
3. **If your change breaks the contract**:
   - **Stop and consult the team** — breaking changes need explicit approval
   - Propose a deprecation path (e.g., keep old signature, add deprecated keyword, wrap old behavior)
   - Document the change in the catalog's changelog section
   - Add deprecation warnings in the code
   - Plan removal timeline (e.g., "remove in v2.0")

## Example: safe vs. unsafe change

```julia
# Original
function top_bottom_boundary_conditions(grid, tracers; top_flux=0.0)
    ...
    return (top=top_bc, bottom=bottom_bc)   # NamedTuple with :top and :bottom
end

# SAFE — adds optional keyword, return shape unchanged
function top_bottom_boundary_conditions(grid, tracers; top_flux=0.0, bottom_flux=0.0)
    ...
    return (top=top_bc, bottom=bottom_bc)
end

# UNSAFE — drops positional arg, breaks all callers
function top_bottom_boundary_conditions(tracers; top_flux=0.0)  # 'grid' removed
    ...
end

# UNSAFE — changes return shape, breaks callers that index .top or .bottom
function top_bottom_boundary_conditions(grid, tracers; top_flux=0.0)
    ...
    return [top_bc, bottom_bc]   # now an Array instead of NamedTuple
end
```

## For AI Agents / Copilot

When an agent (or you) is making changes to FjordSim:

1. **Before implementing changes to exported code:**
   - Consult [api-catalog.md](./api-catalog.md) to understand current contracts
   - Verify proposed changes will not break listed signatures

2. **When adding new exports:**
   - Follow the "Workflow: Adding a New Exported Function/Type" section above
   - Update the catalog as part of the implementation

3. **If uncertain about a change's safety:**
   - Compare the proposed change against the catalog
   - If the change alters arguments, return types, or removes features, flag it and ask for clarification
   - Never silently drop arguments or change return shapes

4. **Integration with project guidelines:**
   - Follow [AGENTS.md](../../AGENTS.md) for surgical, minimal changes
   - This backward-compatibility catalog is the "verifiable success criteria" for API stability
   - Use the catalog to loop and verify that exported interfaces remain unchanged (unless explicitly requested)

---

## References

- **[api-catalog.md](./api-catalog.md)** — Full record of exported functions, types, and signatures (the source of truth)
- **[../copilot-instructions.md](../copilot-instructions.md)** — Project style, naming, and module conventions
- **[../../AGENTS.md](../../AGENTS.md)** — Behavioral guidelines for minimal, non-breaking changes
