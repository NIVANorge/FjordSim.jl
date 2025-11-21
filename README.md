# FjordSim.jl

A framework for ocean simulations built on top of [Oceananigans](https://github.com/CliMA/Oceananigans.jl) and [ClimaOcean](https://github.com/CliMA/ClimaOcean.jl).

One of the main contributions to ClimaOcean is providing a streamlined way to set up simulations. A FjordSim simulation is organized around three core components:

1. **Bathymetry file**
   Contains the domain coordinates (longitude, latitude, depth) along with the bathymetric data.

2. **Forcing file**
   Includes information about sinks and sources (e.g., rivers), boundary conditions, and custom forcings.
   It can also be used to load initial conditions.

3. **Atmospheric data**
   Currently relies on JRA55 from [ClimaOcean](https://github.com/CliMA/ClimaOcean.jl).

Preparation scripts for the Oslofjord are available in the following repository:
[https://github.com/NIVANorge/oslofjord-sim](https://github.com/NIVANorge/oslofjord-sim)

## Installation (from GitHub)

1. Clone the git repository `git clone https://github.com/AquaBadgerForge/FjordSim.jl.git`.
2. Move to the downloaded folder `cd FjordSim.jl`.
3. Run Julia REPL and activate the FjordSim environment `julia --project`.
4. Enter the Pkg REPL by pressing `]` from Julia REPL.
5. Type `instantiate` to 'resolve' a `Manifest.toml` from a `Project.toml` to install and precompile dependency packages.

## Run an example Oslofjord simulation

1. Download the [grid, forcing, atmospheric forcing](https://www.dropbox.com/scl/fo/gc3yc155b5eohi7998wgh/AGN2Yt3HyQ0LlZGImpcca6o?rlkey=x6okc3uxe2avud6sbxgd00l14&st=093llyqp&dl=0).
2. In `FjordSim.jl/app/oslofjord.jl` it is possible to specify the location of the input data files.
By default, the files should be in `$HOME/FjordSim_data/oslofjord/` and `$HOME/FjordSim_data/JRA55/`.
Also, it is possible to specify the results folder destination.
By default, the result will go to `$HOME/FjordSim_results/oslofjord/`.
3. Run `julia --project app/oslofjord.jl`.
This will generate a netcdf result file.

![example_result](./artifacts/phytoplankton_multi.png)
