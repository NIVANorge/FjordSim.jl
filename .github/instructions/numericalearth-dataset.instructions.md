---
description: "Use when adding or modifying a NumericalEarth-style dataset, metadata wrapper, or DataWrangling adapter. Covers the minimum dataset interface, filename dispatch, axis conventions, and verification with Metadata or Metadatum construction."
---
# Creating A New NumericalEarth-Style Dataset

When adding a new dataset or dataset wrapper, follow the existing DataWrangling adapter pattern instead of inventing a new metadata shape. Start from the closest concrete example:

- Use ERA5-style modules for time-varying gridded datasets with variable maps and NetCDF retrieval.
- Use static bathymetry datasets for datasets without a time dimension.
- Use simple wrapper structs like FjordSim's local dataset metadata wrappers when the dataset is already materialized on disk and only needs to satisfy the `Metadata` / `Metadatum` interface.

## Required Interface

Implement the smallest set of methods needed for `Metadata(...)` and `Metadatum(...)` to work without special cases.

- Define a concrete dataset type.
- Implement `default_download_directory(dataset)`.
- Implement `metadata_filename(dataset, name, date, region)`.
- Implement spatial geometry: `longitude_interfaces`, `latitude_interfaces`, and `z_interfaces`.
- Implement `Base.size(dataset, variable)`.
- Implement date coverage for time-varying datasets with `all_dates(dataset, variable)` and, when needed, `first_date` / `last_date`.

If the dataset participates in direct data loading, also implement only the pieces actually used downstream:

- `available_variables(dataset)`.
- `dataset_variable_name(metadata)`.
- `retrieve_data(metadata)`.
- Axis orientation flags such as `reversed_latitude_axis(dataset)` or `reversed_vertical_axis(dataset)` when the source storage order differs from the logical grid order.

## Filename Dispatch Rule

`Metadata(...)` may build filenames from dataset-level dispatch, while `Metadatum(...)` calls `metadata_filename(dataset, variable_name, date, region)` directly. Do not assume one constructor path.

If your dataset wrapper already stores a ready-made filename, add a permissive dataset-level method such as:

```julia
metadata_filename(dataset, args...) = dataset.metadata_filename
```

This avoids `MethodError` when the same dataset is used through different metadata construction paths.

## Shape And Axis Conventions

- Keep returned data aligned with NumericalEarth's logical grid ordering, not the raw file ordering.
- If the source latitude axis is north-to-south, flip it in `retrieve_data` and declare the reversal method consistently.
- If the dataset is logically 2-D but consumed by 3-D machinery, pad or reshape it in the same way existing single-level datasets do.
- Reuse `DatewiseFilename` only when each date maps to a distinct on-disk file.

## Module And API Conventions

- Keep one Julia module per file in `src/` and preserve explicit `export` lists.
- Match existing naming: `PascalCase` for dataset types and `snake_case` for methods.
- Add only the methods the new dataset actually needs. Do not add speculative hooks or abstractions.
- If you are extending FjordSim with a local wrapper rather than editing NumericalEarth itself, keep the wrapper compatible with NumericalEarth's `Metadata` dispatch instead of forking the interface.

## Verification

Before considering the dataset done, verify the real construction path:

1. Construct both `Metadata(variable; dataset=..., dates=..., region=...)` and `Metadatum(variable; dataset=..., date=..., region=...)` when the dataset supports both.
2. Confirm filename generation works for the intended date and region shapes.
3. If the dataset is used for `FieldTimeSeries` or native-grid reads, run the narrowest example or test that exercises one actual variable.
4. Prefer a small example matching the intended usage pattern, such as the bounded ERA5 pressure-level workflow, over broad package-wide validation.

## Non-Goals

- Do not redesign the DataWrangling interface.
- Do not add new required constructor arguments to existing dataset types unless explicitly requested.
- Do not hide axis-order assumptions; encode them in the dataset methods.