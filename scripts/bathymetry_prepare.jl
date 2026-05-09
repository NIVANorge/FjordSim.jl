using Oceananigans
using Oceananigans.Units
using CUDA
using CairoMakie
using FjordSim
using FjordSim.Bathymetry
using Oceananigans.Architectures: on_architecture
using Oceananigans.Fields: interior
using Oceananigans.Grids: x_domain, y_domain

arch = CPU()

z_faces = [
    -450.0,
    -400.0,
    -350.0,
    -300.0,
    -250.0,
    -200.0,
    -150.0,
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
    size = (105, 232, 18),
    halo = (7, 7, 7),
    longitude = (10.2, 11.02),
    latitude = (59.0, 59.93),
    z = z_faces,
)

output_path = joinpath(homedir(), "FjordSim_data", "oslofjord", "bathymetry_105to232+.nc")

result = prepare_geonorge_bathymetry(
    grid;
    output_path,
    raw_resolution_factor = 4,
    padding_cells = 2,
    interpolation_passes = 8,
    major_basins = 1,
    cache = true,
)

bathymetry = result.bottom_height

plot_path = joinpath(homedir(), "FjordSim_data", "oslofjord", "bathymetry_105to232+.png")
isdir(dirname(plot_path)) || mkpath(dirname(plot_path))

cpu_bathymetry = on_architecture(CPU(), bathymetry)
bathymetry_data = Array(interior(cpu_bathymetry, :, :, 1))
Nx, Ny, _ = size(grid)
longitude = collect(range(x_domain(grid)[1], x_domain(grid)[2], length = Nx))
latitude = collect(range(y_domain(grid)[1], y_domain(grid)[2], length = Ny))

figure = Figure(size = (1000, 700))
axis = Axis(figure[1, 1]; xlabel = "Longitude", ylabel = "Latitude", title = "Oslofjord Bathymetry")

plot = heatmap!(axis, longitude, latitude, bathymetry_data; colormap = :deep, colorrange = extrema(bathymetry_data))
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
save(plot_path, figure)

@info "Raw Geonorge bathymetry saved to $(result.raw_path)"
@info "Processed FjordSim bathymetry saved to $(result.output_path)"
@info "Bathymetry plot saved to $plot_path"
