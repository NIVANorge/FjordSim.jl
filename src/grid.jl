import Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid

"""
Return a grid from bathymetry from a netcdf file.
A netcdf file should have 4 variables:
- z_faces - 1d
- h - 2d
- lat and lon - 1d
"""
function ImmersedBoundaryGrid(filepath::String, arch, halo)
    ds = NCDataset(filepath)
    z_faces = ds["z_faces"][:]
    # depths from an nc file are for grid centers
    depth = ds["h"][:, :]
    latitude = compute_faces(ds["lat"][:])
    longitude = compute_faces(ds["lon"][:])

    Nx, Ny = size(depth)
    Nz = length(z_faces)
    # Size should be for grid centers,
    # but z, latitude and langitude should be for faces
    underlying_grid =
        LatitudeLongitudeGrid(arch; size=(Nx, Ny, Nz - 1), halo=halo, z=z_faces, latitude, longitude)
    bathymetry = Field{Center,Center,Nothing}(underlying_grid)
    set!(bathymetry, coalesce.(depth, 0.0))
    fill_halo_regions!(bathymetry)
    return ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bathymetry); active_cells_map=true)
end
