using Oceananigans
using Oceananigans.Units
using CUDA
using FjordSim
using FjordSim.Bathymetry

arch = GPU()

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

output_path = joinpath(homedir(), "FjordSim_data", "oslofjord", "bathymetry_105to232.nc")

result = prepare_geonorge_bathymetry(
    grid;
    output_path,
    raw_resolution_factor = 4,
    padding_cells = 2,
    interpolation_passes = 8,
    major_basins = 1,
)

bathymetry = result.bottom_height

@info "Raw Geonorge bathymetry saved to $(result.raw_path)"
@info "Processed FjordSim bathymetry saved to $(result.output_path)"
