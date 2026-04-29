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

## Exported Symbols (as of current `src/FjordSim.jl`)

These are the current public names. Do not remove any of them from the `export` block:

```
ImmersedBoundaryGrid        # overloaded Oceananigans type/constructor
forcing_from_file
top_bottom_boundary_conditions
coupled_hydrostatic_simulation
recursive_merge
progress
cell_advection_timescale_coupled_model
NORA3PrescribedAtmosphere
MultiYearNORA3
```

From submodules (re-exported or used directly by users):
```
# Utils
compute_faces, safe_execute, extract_z_faces, netcdf_to_jld2, save_fts

# FDatasets
DSForcing, DSResults, last_date

# Forcing
forcing_from_file

# NORA3
NORA3PrescribedAtmosphere
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
