module Grids

export EvenGrid

using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using NCDatasets

import Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
import Oceananigans: LatitudeLongitudeGrid
using ..Utils: compute_faces
using ..Configs: AbstractGridConfig, AbstractBoundaryDataConfig

"""
    domain_grid(config::AbstractGridConfig, architecture)

The bare grid a config describes, with no bathymetry immersed into it.

The geometry every `prepare_*` pipeline regrids onto, and the underlying grid `simulation_grid`
wraps. Declared with no method, so a subtype that does not implement it fails as a `MethodError`
naming the hook rather than silently inheriting some other config's idea of a grid.

Named for the role rather than for a result type — it replaced a method on
`Oceananigans.LatitudeLongitudeGrid`, under which a config describing a rectilinear, single-column
or curvilinear grid could not have implemented the hook without lying about what it returns.
"""
function domain_grid end

"""
    simulation_grid(config::AbstractGridConfig, bathymetry_file, architecture)

The grid the simulation runs on: `domain_grid`'s geometry with the processed bathymetry immersed
into it.

Separate from `domain_grid` because the two are needed at different moments and from different
inputs — the prepare steps have no bathymetry file yet, while `build_simulation` must read back the
very grid the bathymetry was written on.
"""
function simulation_grid end

"""
The four lateral boundaries a regional domain can be open on, in the order every error message
lists them.

`Bathymetry` clears land along the open ones, `Regridding` restores their velocity face rows,
`BoundaryData` regrids along them, and every `Val{edge}` dispatch of `OpenLateralBoundaryFromData`
falls back to an `ArgumentError` naming this tuple.
"""
const LATERAL_EDGES = (:south, :north, :west, :east)

"""
    validate_open_edge(edge)

Return `edge` if it names a lateral boundary, else throw. Checked up front so a typo fails
before any regridding rather than deep inside the variable loop.

Checked here rather than in each pipeline so one typo produces one error message wherever it is
written — in a boundary data config's `open_edges`, or in the `edges` a prepare step is handed.
"""
function validate_open_edge(edge)
    edge in LATERAL_EDGES ||
        throw(ArgumentError("open_edge must be one of $LATERAL_EDGES, got :$edge"))
    return edge
end

"""
    lateral_edges(edges)

Normalize whatever a caller names its open edges as into a validated `Vector{Symbol}`: one `Symbol`,
an iterable of them, or `nothing` for none.

Every consumer iterates, so a domain open on one edge, on all four, or on none is the same code path
rather than a scalar case and a plural one. Duplicates are rejected — a repeated edge would prepare
and write the same variables twice, and mean nothing the second time.
"""
lateral_edges(::Nothing) = Symbol[]
lateral_edges(edge::Symbol) = [validate_open_edge(edge)]

function lateral_edges(edges)
    normalized = Symbol[validate_open_edge(edge) for edge in edges]
    allunique(normalized) ||
        throw(ArgumentError("open_edges names an edge more than once: $normalized"))
    return normalized
end

"""
    open_edges(config)

Which lateral boundaries the domain is open on, as named by an `AbstractBoundaryDataConfig` — a
`Vector{Symbol}`, empty for a setup naming no open-boundary data, whose lateral boundaries are
therefore all closed walls.

Always a collection, never a bare `Symbol`, so a domain open on one edge, on all four, or on none is
one code path in every consumer rather than a scalar case and a plural one. A built-in config's
keyword constructor accepts either shape and normalizes with `lateral_edges`.

Every consumer reads it through this accessor, so the edges are stated in exactly one place — the
config that describes the boundary.
"""
open_edges(::Nothing) = Symbol[]
open_edges(config::AbstractBoundaryDataConfig) = config.open_edges

"""
    EvenGrid

Configuration for a `LatitudeLongitudeGrid` that is regularly spaced in longitude and
latitude, with explicit vertical faces.

# Fields
- `size`: `(Nx, Ny, Nz)` number of cells.
- `halo`: `(Hx, Hy, Hz)` number of halo cells.
- `longitude`: `(west, east)` longitude bounds in degrees.
- `latitude`: `(south, north)` latitude bounds in degrees.
- `z_faces`: Vertical face coordinates in increasing order (bottom to top).
- `float_type`: The float type the simulation grid, and so the model, runs in. `Float32` by
  default, because it is much faster on a GPU. Only `simulation_grid` reads it: `domain_grid`, which
  the prepare steps regrid onto, stays in Oceananigans' own default so it matches the source grids
  they interpolate from.

The float type is stated here, per setup, rather than by setting `Oceananigans.defaults.FloatType`:
that is a process-wide global, so a setup that set it would change the precision of every config
built after it. A setup passes the same `FT` to every model component it constructs — `WENO(FT)`,
`CATKEVerticalDiffusivity(FT; ...)` and so on — since those read that global too when not given one.
"""
Base.@kwdef mutable struct EvenGrid <: AbstractGridConfig
    size::NTuple{3,Int}
    halo::NTuple{3,Int}
    longitude::NTuple{2,Float64}
    latitude::NTuple{2,Float64}
    z_faces::Vector{Float64}
    float_type::DataType = Float32
end

"""
    domain_grid(config::EvenGrid, architecture)

The `LatitudeLongitudeGrid` an `EvenGrid` describes.
"""
domain_grid(config::EvenGrid, architecture) = LatitudeLongitudeGrid(architecture, config)

"""
    simulation_grid(config::EvenGrid, bathymetry_file, architecture)

The `ImmersedBoundaryGrid` an `EvenGrid` runs on.
"""
simulation_grid(config::EvenGrid, bathymetry_file, architecture) =
    ImmersedBoundaryGrid(bathymetry_file, architecture, config.halo, config.float_type)

function LatitudeLongitudeGrid(architecture, config::EvenGrid)
    return LatitudeLongitudeGrid(
        architecture;
        size      = config.size,
        halo      = config.halo,
        longitude = config.longitude,
        latitude  = config.latitude,
        z         = config.z_faces,
    )
end

"""
        ImmersedBoundaryGrid(filepath::String, architecture, halo, FT = Oceananigans.defaults.FloatType)

Construct an immersed-boundary `LatitudeLongitudeGrid` with float type `FT` from a bathymetry NetCDF
file.

The preferred input file layout is:

- dimensions:
    - `lon`: horizontal x/longitude centers, length `Nx`
    - `lat`: horizontal y/latitude centers, length `Ny`
    - `zf`: vertical faces, length `Nz + 1`
- variables:
    - `lon(lon)`: longitude center coordinates
    - `lat(lat)`: latitude center coordinates
    - `z_faces(zf)`: vertical face coordinates
    - `h(lon, lat)`: 2D bathymetry stored on horizontal cell centers

Preferred bathymetry convention:

- `h` should store bottom height, meaning negative values below sea level and
    positive values over land, matching Oceananigans `PartialCellBottom`.

Legacy compatibility retained by this loader:

- If `h` contains only non-negative values, it is interpreted as positive depth
    and converted internally to negative bottom height.
- If the file uses the older axis association where `lat` has length `Nx` and
    `lon` has length `Ny`, that layout is still accepted.

Notes on file construction:

- `lon` and `lat` should contain cell-center coordinates, not faces.
- `z_faces` should contain vertical face coordinates in increasing order from
    bottom to top, for example `[-450.0, ..., -1.0, 0.0]`.
- `h` must have shape `(Nx, Ny)` on tracer centers.
- Missing values in `h` are treated as `0.0` during grid construction.

In short, new files should be written as `lon`, `lat`, `z_faces`, and
`h(lon, lat)` using bottom height, while older files with swapped horizontal
axis vectors or positive depth values are still supported.
"""
function ImmersedBoundaryGrid(filepath::String, architecture, halo, FT = Oceananigans.defaults.FloatType)
    ds = NCDataset(filepath)
    z_faces = ds["z_faces"][:]
    bottom_height = ds["h"][:, :]
    lat_centers = ds["lat"][:]
    lon_centers = ds["lon"][:]

    finite_bottom = bottom_height[isfinite.(bottom_height)]
    if !isempty(finite_bottom) && minimum(finite_bottom) >= 0
        bottom_height = -bottom_height
    end

    Nx, Ny = size(bottom_height)

    if length(lon_centers) == Nx && length(lat_centers) == Ny
        longitude = compute_faces(lon_centers)
        latitude = compute_faces(lat_centers)
    elseif length(lat_centers) == Nx && length(lon_centers) == Ny
        latitude = compute_faces(lat_centers)
        longitude = compute_faces(lon_centers)
    else
        close(ds)
        error("Bathymetry axes do not match h dimensions in $filepath.")
    end

    Nz = length(z_faces)
    # Size should be for grid centers,
    # but z, latitude and langitude should be for faces
    underlying_grid =
        LatitudeLongitudeGrid(architecture, FT; size=(Nx, Ny, Nz - 1), halo=halo, z=z_faces, latitude, longitude)
    bathymetry = Field{Center,Center,Nothing}(underlying_grid)
    set!(bathymetry, coalesce.(bottom_height, 0.0))
    fill_halo_regions!(bathymetry)
    close(ds)
    return ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bathymetry); active_cells_map=true)
end

end  # module Grids
