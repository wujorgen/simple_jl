module BoundaryConditions

export set_velocity_boundaries!
export set_pressure_boundaries!

using ..Types

"""
    set_velocity_boundaries!
Sets the velocity at the boundaries of the domain.
"""
function set_velocity_boundaries!(velocity::VectorFieldStaggered2D)

end

"""
    set_pressure_boundaries!
Sets the pressure at the boundaries of the domain.
"""
function set_pressure_boundaries!(pressure::Array{Float64,2})

end

end