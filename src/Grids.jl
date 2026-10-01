module Grids

export ImmersedBoundaryGrid, EvenGrid, domain_grid, simulation_grid

using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using NCDatasets

import Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
import Oceananigans: LatitudeLongitudeGrid
using ..Utils: compute_faces
using ..Configs
using ..Configs: AbstractGridConfig, domain_grid, simulation_grid

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

The `LatitudeLongitudeGrid` an `EvenGrid` describes. See `Configs.domain_grid`.
"""
Configs.domain_grid(config::EvenGrid, architecture) = LatitudeLongitudeGrid(architecture, config)

"""
    simulation_grid(config::EvenGrid, bathymetry_file, architecture)

The `ImmersedBoundaryGrid` an `EvenGrid` runs on. See `Configs.simulation_grid`.
"""
Configs.simulation_grid(config::EvenGrid, bathymetry_file, architecture) =
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
