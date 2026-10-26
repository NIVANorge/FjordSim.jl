using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using FjordSim

config = fjord_config(joinpath(@__DIR__, "..", "examples", "oslofjorden.jl"))
prepare_bathymetry(config)
# add_rivers(config)
# download_boundary_data(config)
# prepare_boundary_data(config)
# download_atmosphere(config)
# prepare_atmosphere(config)

simulation = build_simulation(config)
