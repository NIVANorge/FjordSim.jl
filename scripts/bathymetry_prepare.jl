##
using Oceananigans
using Oceananigans.Units
using CUDA
using CairoMakie
using Printf
using Statistics
using FjordSim
using FjordSim.Bathymetry
using Oceananigans.Architectures: on_architecture
using Oceananigans.Fields: interior
using Oceananigans.Grids: x_domain, y_domain

function plot_bathymetry(grid, bathymetry; plot_path, title = "Bathymetry", figure_size = (1000, 700))
    isdir(dirname(plot_path)) || mkpath(dirname(plot_path))

    cpu_bathymetry = on_architecture(CPU(), bathymetry)
    bathymetry_data = Array(interior(cpu_bathymetry, :, :, 1))

    Nx, Ny, _ = size(grid)
    longitude = collect(range(x_domain(grid)[1], x_domain(grid)[2], length = Nx))
    latitude = collect(range(y_domain(grid)[1], y_domain(grid)[2], length = Ny))

    figure = Figure(size = figure_size)
    axis = Axis(figure[1, 1]; xlabel = "Longitude", ylabel = "Latitude", title)

    plot = heatmap!(axis, longitude, latitude, bathymetry_data; colormap = :deep, colorrange = extrema(bathymetry_data))
    land_mask = ifelse.(bathymetry_data .>= 0, 1.0f0, NaN32)
    heatmap!(
        axis,
        longitude,
        latitude,
        land_mask;
        colormap = [:ivory, :ivory],
        colorrange = (0, 1),
        nan_color = RGBAf(0, 0, 0, 0),
    )
    Colorbar(figure[1, 2], plot; label = "Bottom height (m)")
    contour!(
        axis,
        longitude,
        latitude,
        bathymetry_data;
        levels = -collect(25:25:300),
        color = (:white, 0.35),
        linewidth = 1,
    )
    contour!(
        axis,
        longitude,
        latitude,
        bathymetry_data;
        levels = [0.0],
        color = :black,
        linewidth = 4,
    )
    save(plot_path, figure)

    return plot_path
end
##

##
arch = CPU()

z_faces = [
    -100.0,
    -75.0,
    -50.0,
    -25.0,
    -15.0,
    -10.0,
    -7.5,
    -5.0,
    -3.0,
    -2.0,
    -1.0,
    0.0,
]

grid = LatitudeLongitudeGrid(
    arch,
    size = (150, 200, 11),
    halo = (7, 7, 7),
    longitude = (10.20, 10.45),
    latitude = (59.58, 59.75),
    z = z_faces,
)
##

##
dx = xspacings(grid)
dy = yspacings(grid)
dz = zspacings(grid)

println("Grid cell side extents (m):")
@printf("  N = (%d, %d, %d)\n", grid.Nx, grid.Ny, grid.Nz)
@printf("  Δx min/max/mean = %.3f / %.3f / %.3f\n", minimum(dx), maximum(dx), mean(dx))
@printf("  Δy min/max/mean = %.3f / %.3f / %.3f\n", minimum(dy), maximum(dy), mean(dy))
@printf("  Δz min/max/mean = %.3f / %.3f / %.3f\n", minimum(dz), maximum(dz), mean(dz))
##

##
name = "drammensfjorden"
bathymetry_name = "bathymetry_$(name)"
output_path = joinpath(homedir(), "FjordSim_data", name, "$(bathymetry_name).nc")
geodatabase_path = joinpath(homedir(), "FjordSim_data", "Basisdata_0000_Norge_25833_Dybdedata_FGDB.gdb")

result = prepare_geonorge_bathymetry(
    grid;
    output_path,
    geodatabase_path,
    raw_resolution_factor = 2,
    padding_cells = 2,
    include_contours = false,
    interpolation_passes = 1,
    major_basins = 1,
    cache = false,
)
##

##
bathymetry = result.bottom_height

plot_path = joinpath(homedir(), "FjordSim_data", name, "$(bathymetry_name).png")
plot_bathymetry(grid, bathymetry; plot_path)

@info "Raw Geonorge bathymetry saved to $(result.raw_path)"
@info "Processed FjordSim bathymetry saved to $(result.output_path)"
@info "Bathymetry plot saved to $plot_path"
##
