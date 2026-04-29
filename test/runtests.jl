using FjordSim
using Test
using NCDatasets
using Oceananigans
using Oceananigans.BoundaryConditions: FluxBoundaryCondition

@testset "Backward Compatibility — API Exports" begin
    # Verify all exported symbols are present in the public interface
    exported_symbols = [
        :ImmersedBoundaryGrid,
        :forcing_from_file,
        :top_bottom_boundary_conditions,
        :coupled_hydrostatic_simulation,
        :recursive_merge,
        :progress,
        :cell_advection_timescale_coupled_model,
        :NORA3PrescribedAtmosphere,
        :MultiYearNORA3,
    ]

    for sym in exported_symbols
        @test isdefined(FjordSim, sym)  # All core exports must be defined
    end

    # Verify submodule exports
    submodule_exports = [
        (:Utils, :compute_faces),
        (:Utils, :safe_execute),
        (:Utils, :extract_z_faces),
        (:Utils, :netcdf_to_jld2),
        (:Utils, :save_fts),
        (:FDatasets, :DSForcing),
        (:FDatasets, :DSResults),
        (:FDatasets, :last_date),
    ]

    for (module_name, sym) in submodule_exports
        mod = getfield(FjordSim, module_name)
        @test isdefined(mod, sym)  # All submodule exports must be defined
    end
end

@testset "Backward Compatibility — Function Signatures" begin
    mktempdir() do tmp
        # Create minimal test files
        bathymetry_path = joinpath(tmp, "bathymetry.nc")
        ds = NCDataset(bathymetry_path, "c")
        defDim(ds, "x", 2)
        defDim(ds, "y", 2)
        defDim(ds, "zf", 3)
        z_faces = defVar(ds, "z_faces", Float64, ("zf",))
        h = defVar(ds, "h", Float64, ("x", "y"))
        lat = defVar(ds, "lat", Float64, ("x",))
        lon = defVar(ds, "lon", Float64, ("y",))
        z_faces[:] = [-20.0, -10.0, 0.0]
        h[:, :] = [10.0 10.0; 10.0 10.0]
        lat[:] = [59.0, 60.0]
        lon[:] = [10.0, 11.0]
        close(ds)

        arch = CPU()
        grid = ImmersedBoundaryGrid(bathymetry_path, arch, (1, 1, 1))

        # Test top_bottom_boundary_conditions signature
        @test_nowarn top_bottom_boundary_conditions(; grid, bottom_drag_coefficient = 0.003)
        bcs = top_bottom_boundary_conditions(; grid, bottom_drag_coefficient = 0.003)
        
        # Verify return type: should be a NamedTuple with u, v, T, S
        @test haskey(bcs, :u)  # Return must include u
        @test haskey(bcs, :v)  # Return must include v
        @test haskey(bcs, :T)  # Return must include T
        @test haskey(bcs, :S)  # Return must include S
        
        # Verify nested structure for u and v: (top, bottom)
        @test haskey(bcs.u, :top)  # u must have top
        @test haskey(bcs.u, :bottom)  # u must have bottom
        @test haskey(bcs.v, :top)  # v must have top
        @test haskey(bcs.v, :bottom)  # v must have bottom
        
        # Verify return types are callable boundary conditions (not checking exact type, just that they exist and are not nothing)
        @test !isnothing(bcs.u.top)  # u.top must exist
        @test !isnothing(bcs.u.bottom)  # u.bottom must exist
        @test !isnothing(bcs.v.top)  # v.top must exist
        @test !isnothing(bcs.v.bottom)  # v.bottom must exist
        @test !isnothing(bcs.T.top)  # T.top must exist
        @test !isnothing(bcs.S.top)  # S.top must exist
    end
end

@testset "Backward Compatibility — Struct Fields" begin
    mktempdir() do tmp
        # Create minimal NORA3 file
        nora3_filename = "NORA3_test.nc"
        nora3_path = joinpath(tmp, nora3_filename)
        ds = NCDataset(nora3_path, "c")
        defDim(ds, "x", 2)
        defDim(ds, "y", 2)
        defDim(ds, "time", 2)
        air_temperature_2m = defVar(ds, "air_temperature_2m", Float64, ("x", "y", "time"))
        time_var = defVar(ds, "time", Float64, ("time",))
        air_temperature_2m[:, :, :] .= 273.15  # Use broadcasting assignment
        time_var[:] = [0.0, 3600.0]
        close(ds)

        # Test MultiYearNORA3 struct fields
        nora3 = MultiYearNORA3(nora3_filename, tmp)
        @test hasfield(typeof(nora3), :metadata_filename)  # Field required
        @test hasfield(typeof(nora3), :default_download_directory)  # Field required
        @test hasfield(typeof(nora3), :size)  # Field required
        @test hasfield(typeof(nora3), :all_dates)  # Field required
        
        # Verify field types
        @test isa(nora3.metadata_filename, String)  # metadata_filename is String
        @test isa(nora3.default_download_directory, String)  # default_download_directory is String
        @test isa(nora3.size, Tuple)  # size is Tuple
        @test isa(nora3.all_dates, Vector)  # all_dates is Vector
    end
end

@testset "FjordSim.jl" begin
    mktempdir() do tmp
        bathymetry_path = joinpath(tmp, "bathymetry.nc")
        forcing_path = joinpath(tmp, "forcing.nc")
        nora3_filename = "NORA3_test.nc"
        nora3_path = joinpath(tmp, nora3_filename)

        # Minimal bathymetry file for Grids.ImmersedBoundaryGrid.
        ds = NCDataset(bathymetry_path, "c")
        defDim(ds, "x", 2)
        defDim(ds, "y", 2)
        defDim(ds, "zf", 3)
        z_faces = defVar(ds, "z_faces", Float64, ("zf",))
        h = defVar(ds, "h", Float64, ("x", "y"))
        lat = defVar(ds, "lat", Float64, ("x",))
        lon = defVar(ds, "lon", Float64, ("y",))
        z_faces[:] = [-20.0, -10.0, 0.0]
        h[:, :] = fill(10.0, (2, 2))  # Create and assign a 2x2 matrix of 10.0
        lat[:] = [59.0, 60.0]
        lon[:] = [10.0, 11.0]
        close(ds)

        # Minimal forcing file for Forcing.forcing_from_file.
        ds = NCDataset(forcing_path, "c")
        defDim(ds, "Nx", 2)
        defDim(ds, "Ny", 2)
        defDim(ds, "Nz", 2)
        defDim(ds, "time", 2)
        time = defVar(ds, "time", Float64, ("time",))
        u = defVar(ds, "u", Float64, ("Nx", "Ny", "Nz", "time"))
        v = defVar(ds, "v", Float64, ("Nx", "Ny", "Nz", "time"))
        T = defVar(ds, "T", Float64, ("Nx", "Ny", "Nz", "time"))
        u_lambda = defVar(ds, "u_lambda", Float64, ("Nx", "Ny", "Nz", "time"))
        v_lambda = defVar(ds, "v_lambda", Float64, ("Nx", "Ny", "Nz", "time"))
        T_lambda = defVar(ds, "T_lambda", Float64, ("Nx", "Ny", "Nz", "time"))
        time[:] = [0.0, 3600.0]
        u[:, :, :, :] .= 0.0
        v[:, :, :, :] .= 0.0
        T[:, :, :, :] .= 4.0
        u_lambda[:, :, :, :] .= 2.0
        v_lambda[:, :, :, :] .= -2.0
        T_lambda[:, :, :, :] .= 0.1
        close(ds)

        # Minimal NORA3 file for MultiYearNORA3 dataset construction.
        ds = NCDataset(nora3_path, "c")
        defDim(ds, "x", 2)
        defDim(ds, "y", 2)
        defDim(ds, "time", 2)
        air_temperature_2m = defVar(ds, "air_temperature_2m", Float64, ("x", "y", "time"))
        time = defVar(ds, "time", Float64, ("time",))
        air_temperature_2m[:, :, :] .= 273.15
        time[:] = [0.0, 3600.0]
        close(ds)

        arch = CPU()
        grid = @test_nowarn ImmersedBoundaryGrid(bathymetry_path, arch, (0, 0, 0))
        @test_nowarn top_bottom_boundary_conditions(; grid, bottom_drag_coefficient = 0.003)
        @test_nowarn forcing_from_file(; grid, filepath = forcing_path, tracers = (:T,))
        @test_nowarn MultiYearNORA3(nora3_filename, tmp)
    end
end
