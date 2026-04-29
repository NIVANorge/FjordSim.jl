using FjordSim
using Test

@testset "FjordSim.jl" begin
    @test isdefined(FjordSim, :coupled_hydrostatic_simulation)
    @test isdefined(FjordSim, :forcing_from_file)
    @test isdefined(FjordSim, :top_bottom_boundary_conditions)
    @test isdefined(FjordSim, :NORA3PrescribedAtmosphere)
    @test isdefined(FjordSim, :ImmersedBoundaryGrid)
    @test isdefined(FjordSim, :Atmospheres)
    @test isdefined(FjordSim.Atmospheres, :NORA3)

    example_path = joinpath(dirname(@__DIR__), "examples", "oslofjord.jl")
    @test isfile(example_path)
end
