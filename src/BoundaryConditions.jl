module BoundaryConditions

export top_bottom_boundary_conditions

using Oceananigans
using Oceananigans.BoundaryConditions: FluxBoundaryCondition
using NumericalEarth.Oceans: u_quadratic_bottom_drag, v_quadratic_bottom_drag

""" Return a named tuple with boundary conditions """
function top_bottom_boundary_conditions(; grid, bottom_drag_coefficient)
    top_zonal_momentum_flux = τx = Field{Face,Center,Nothing}(grid)
    top_meridional_momentum_flux = τy = Field{Center,Face,Nothing}(grid)
    top_ocean_heat_flux = Jᵀ = Field{Center,Center,Nothing}(grid)
    top_salt_flux = Jˢ = Field{Center,Center,Nothing}(grid)

    u_bot_bc =
        FluxBoundaryCondition(u_quadratic_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)
    v_bot_bc =
        FluxBoundaryCondition(v_quadratic_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)

    return (
        u = (top = FluxBoundaryCondition(τx), bottom = u_bot_bc),
        v = (top = FluxBoundaryCondition(τy), bottom = v_bot_bc),
        T = (top = FluxBoundaryCondition(Jᵀ),),
        S = (top = FluxBoundaryCondition(Jˢ),),
    )
end

end  # module BoundaryConditions
