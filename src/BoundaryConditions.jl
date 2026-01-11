module BoundaryConditions

export set_velocity_boundaries!
export set_pressure_boundaries!

using ..Types

"""
    set_velocity_boundaries!
Sets the velocity at the boundaries of the domain.
"""
function set_velocity_boundaries!(velocity::VectorFieldStaggered2D;
    top::Tuple{Float64,Float64}=(0.0, 0.0), bottom::Tuple{Float64,Float64}=(0.0, 0.0), left::Tuple{Float64,Float64}=(0.0, 0.0), right::Tuple{Float64,Float64}=(0.0, 0.0))
    # u-velocity
    # left
    velocity.u[1, :] .= 2 * left[1] .- velocity.u[2, :]
    # right
    velocity.u[end, :] .= 2 * right[1] .- velocity.u[end-1, :]
    # top
    velocity.u[:, end] .= top[1]
    # bottom
    velocity.u[:, 1] .= bottom[1]
    # v-velocity
    # top
    velocity.v[:, end] .= 2 * top[2] .- velocity.v[:, end-1]
    # bottom
    velocity.v[:, 1] .= 2 * bottom[2] .- velocity.v[:, 2]
    # left
    velocity.v[1, :] .= left[2]
    # right
    velocity.v[end, :] .= right[2]
end

"""
    set_pressure_boundaries!
Sets the pressure at the boundaries of the domain. Hardcoded for no gradients.
"""
function set_pressure_boundaries!(pressure::Array{Float64,2})
    # bottom
    pressure[:, 1] = pressure[:, 2]
    # top
    pressure[:, end] = pressure[:, end-1]
    # left
    pressure[1, :] = pressure[2, :]
    # right
    pressure[end, :] = pressure[end-1, :]
end

end