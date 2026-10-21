module FjordSim

export
    # the setup container and the helpers a setup builds its paths with
    FjordConfig,
    fjord_config,
    fjord_data_root,
    fjord_results_root,
    # pipeline steps
    prepare_bathymetry,
    download_forcing,
    prepare_forcing,
    download_initial_state,
    prepare_initial_state,
    add_rivers,
    download_boundaries,
    prepare_boundaries,
    download_atmosphere,
    prepare_atmosphere,
    build_simulation,
    run_simulation,
    validate_simulation,
    # built-in sources
    EvenGrid,
    DybdedataConfig,
    NorKystConfig,
    NorKystHindcastConfig,
    NorKystInitialStateConfig,
    NVERiversConfig,
    NVERiver,
    NorKystBoundariesConfig,
    NorKystHindcastBoundariesConfig,
    NORA3Config,
    # simulation
    SimulationConfig,
    CoupledHydrostaticSimulation,
    SplitExplicitFreeSurfaceConfig,
    BoundarySponge,
    RiverMouthMixing,
    SnapshotWriter,
    FieldSnapshotWriter,
    Station,
    StationWriter,
    FieldStationWriter,
    CheckpointWriter,
    ProgressCallback,
    progress,
    AdaptiveTimeStep,
    FromForcing,
    FromInitialState,
    FromResults,
    # boundary conditions
    MergedBoundaryConditions,
    AirSeaFluxes,
    QuadraticBottomDrag,
    OpenLateralBoundaryFromData

using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Units
using Oceananigans.Utils
using NumericalEarth
using NCDatasets
using Adapt

include("Configs.jl")
include("Utils.jl")
# Grids comes before every pipeline: each builds its target grid with `domain_grid`, and most read
# the edge vocabulary.
include("Grids.jl")
# Plotting comes before the pipelines that call it: their setup-level drivers plot as their last
# step, and Plotting itself only needs Configs and Grids.
include("Plotting.jl")
include("Bathymetry/Bathymetry.jl")
include("Regridding.jl")
include("Atmospheres/Atmospheres.jl")
include("DataSources/NorKyst.jl")
include("Forcing/Forcing.jl")
include("InitialState/InitialState.jl")
include("Rivers/Rivers.jl")
include("Boundaries/BoundaryConditions/BoundaryConditions.jl")
include("Boundaries/BoundaryData/BoundaryData.jl")
include("Simulations.jl")
include("Validation/Validation.jl")
# CLI names every driver, so it comes last.
include("CLI.jl")

using .Configs
using .Utils
using .Grids
using .Plotting
using .Bathymetry
using .Regridding
using .Atmospheres
using .Forcing
using .InitialState
using .Rivers
using .BoundaryConditions
using .BoundaryData
using .Simulations
using .Validation
using .CLI

#####
##### Upstream workaround: net ocean fluxes for a setup with no sea ice
#####
#
# NumericalEarth's own flux assembler, with one value given a concrete float type, because without
# it `run_simulation` cannot compile on a GPU at all.
#
# Upstream reads the sea-ice contribution as
# `computed_fluxes(coupled_model.interfaces.sea_ice_ocean_interface)`. With no sea-ice component
# that interface is `nothing`, so the result is a `ZeroFluxes` built from `ZeroField()`, which
# defaults to `ZeroField{Int64}`. Indexing one yields an `Int64` zero, and
# `NumericalEarth/src/Oceans/assemble_net_ocean_fluxes.jl`'s
#
#     Jˢ[i, j, 1] = ifelse(inactive, zero(grid), Jˢio)
#
# then infers `Union{FT, Int64}` rather than `FT`: `zero(grid)` is a float, `Jˢio` is that integer,
# and `ifelse` does not promote its branches. Every other use of these zeros is arithmetic
# (`Jᵀao + Jᵀio`), which *does* promote, so that one bare `ifelse` is the only unstable spot.
#
# On the CPU that is a small performance wart. On the GPU it is fatal: Julia emits the dynamic
# branch as a read of a `const` global, which lowers to a load that is both `atomic` and tagged
# `!invariant.load`, and LLVM's NVPTX back-end segfaults selecting exactly that pair in
# `NVPTXDAGToDAGISel::tryLDG`. It does so in 22.1.5, 22.1.7 and 23.1.1 — every
# `NVPTX_LLVM_Backend_jll` in the General registry — so no version bump escapes it, and neither
# does any `sm_*`, PTX ISA or stack-size setting. The run dies in codegen before the first time
# step. The offending load is dead code, kept only because `volatile` forbids deleting it.
#
# A setup with a real sea-ice component never reaches this: its fluxes are ordinary float fields,
# so the line is already stable. That is why NumericalEarth's own ocean-sea-ice tests never see it.
#
# `computed_fluxes(::Nothing)` is where the fix belongs, but Julia forbids overwriting another
# module's method during precompilation, and the interface field cannot be swapped either: it is
# typed `Nothing` on the built model, and it is *that* `Nothing` which dispatches
# `compute_sea_ice_ocean_fluxes!` to a no-op. So the narrowest hook left is this caller, dispatched
# on `NoSeaIceInterfaceModel` — precisely the configuration that breaks — which leaves every other
# setup on upstream's code path. The body is upstream's apart from `sea_ice_ocean_fluxes`, and the
# kernel it launches is still upstream's.
#
# Delete this and `typed_zero_fluxes` once NumericalEarth types those zeros itself.

# The body below is a copy, so an upgrade could silently leave it running stale flux assembly.
# Warn while precompiling — the moment the manifest changes — rather than fail a run much later.
if pkgversion(NumericalEarth) != v"0.8.1"
    @warn """FjordSim carries a copy of NumericalEarth's `update_net_ocean_fluxes!` for setups \
             without sea ice (see the comment above it in src/FjordSim.jl). It was written against \
             NumericalEarth 0.8.1; this is $(pkgversion(NumericalEarth)). Check whether the \
             upstream `ZeroField()` zeros are typed yet — if they are, delete the copy."""
end

function NumericalEarth.Oceans.update_net_ocean_fluxes!(
    coupled_model::NumericalEarth.EarthSystemModels.NoSeaIceInterfaceModel,
    ocean_model,
    grid,
)
    sea_ice = coupled_model.sea_ice
    arch = Oceananigans.Architectures.architecture(grid)
    clock = coupled_model.clock

    net_ocean_fluxes = coupled_model.interfaces.net_fluxes.ocean
    # The one line that differs from upstream: zeros typed to the grid, not `ZeroField{Int64}`.
    sea_ice_ocean_fluxes = typed_zero_fluxes(eltype(grid))

    atmos_ocean_fluxes = NumericalEarth.Oceans.atmos_ocean_flux(coupled_model)
    rainfall = NumericalEarth.Oceans.rainfall_flux(coupled_model)
    snowfall = NumericalEarth.Oceans.snowfall_flux(coupled_model)

    land_exchanger = coupled_model.interfaces.exchanger.land
    runoff_freshwater_flux = NumericalEarth.Oceans.land_freshwater_flux(land_exchanger, :runoff_freshwater_flux)
    iceberg_freshwater_flux = NumericalEarth.Oceans.land_freshwater_flux(land_exchanger, :iceberg_freshwater_flux)

    ice_concentration = NumericalEarth.EarthSystemModels.sea_ice_concentration(sea_ice)
    intercepted_snowfall_flux = NumericalEarth.EarthSystemModels.intercepted_snowfall(sea_ice)
    ocean_surface_temperature = NumericalEarth.EarthSystemModels.ocean_surface_temperature(ocean_model)
    ocean_properties = coupled_model.interfaces.ocean_properties

    Oceananigans.Utils.launch!(arch, grid, :xy,
                               NumericalEarth.Oceans._assemble_net_ocean_fluxes!,
                               net_ocean_fluxes,
                               grid,
                               clock,
                               atmos_ocean_fluxes,
                               sea_ice_ocean_fluxes,
                               ocean_surface_temperature,
                               ice_concentration,
                               rainfall,
                               snowfall,
                               intercepted_snowfall_flux,
                               runoff_freshwater_flux,
                               iceberg_freshwater_flux,
                               ocean_properties)

    if grid isa Oceananigans.ImmersedBoundaries.MutableGridOfSomeKind
        Oceananigans.BoundaryConditions.fill_halo_regions!(net_ocean_fluxes.η)
    end

    return nothing
end

"""
    typed_zero_fluxes(FT)

`NumericalEarth`'s `ZeroFluxes` with every field a `ZeroField{FT}` instead of the `ZeroField{Int64}`
its own zero-argument constructor produces. The field count is read from the struct rather than
written out, because upstream grows it as interfaces gain flux names.
"""
function typed_zero_fluxes(::Type{FT}) where FT
    ZeroFluxes = NumericalEarth.EarthSystemModels.InterfaceComputations.ZeroFluxes
    return ZeroFluxes(ntuple(_ -> Oceananigans.Fields.ZeroField(FT), fieldcount(ZeroFluxes))...)
end

"""
    main(args)

Entry point for `julia --project -m FjordSim SUBCOMMAND --config SETUP`. Returns a process exit
code; see `FjordSim.CLI.USAGE` for the subcommands.

Deliberately *not* exported. Julia's startup runs `Main.main` after a script's body whenever that
binding resolves to an entry point, so exporting this would make every `using FjordSim` in a
script — `test/runtests.jl`, or a config file run directly rather than through `--config` — run the
CLI on the way out.
"""
function main(args)
    return CLI.main(args)
end

# Bare `@main`, applied after the definition: `@main function main(args) ... end` expands to a
# *call*, which would run the CLI while the package precompiles.
@main

end  # module
