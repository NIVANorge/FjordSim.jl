using Oceananigans
using Oceananigans.Units
using Dates
using CUDA
import NumericalEarth

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

bathymetry = NumericalEarth.regrid_bathymetry(grid)
