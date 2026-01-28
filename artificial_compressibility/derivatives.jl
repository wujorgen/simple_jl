module Derivatives

export calc_dudtau!, calc_dvdtau!, calc_dpdtau!

"""
    calc_dudtau!

du/dtau + u du/dx + v du/dy = -1/rho * dp/dx + mu/rho * (ddu/ddx + ddu/ddy)

"""
function calc_dudtau!(dudtau::Array{Float64,2}, u::Array{Float64,2}, v::Array{Float64,2}, p::Array{Float64,2}, dx::Float64, dy::Float64, rho::Float64, mu::Float64)
    dudx = (u[3:end, 2:end-1] - u[1:end-2, 2:end-1]) / (2 * dx)
    dudy = (u[2:end-1, 3:end] - u[2:end-1, 1:end-2]) / (2 * dy)
    dpdx = (p[3:end, 2:end-1] - p[1:end-2, 2:end-1]) / (2 * dx)
    dduddx = (u[3:end, 2:end-1] - 2 * u[2:end-1, 2:end-1] + u[1:end-2, 2:end-1]) / (dx^2)
    dduddy = (u[2:end-1, 3:end] - 2 * u[2:end-1, 2:end-1] + u[2:end-1, 1:end-2]) / (dy^2)

    dudtau[2:end-1, 2:end-1] = -(u[2:end-1, 2:end-1] .* dudx .+ v[2:end-1, 2:end-1] .* dudy) - 1 / rho * dpdx + (mu / rho) * (dduddx + dduddy)
end

"""
    calc_dvdtau!

dv/dtau + u dv/dx + v dv/dy = -1/rho * dp/dy + mu/rho * (ddv/ddx + ddv/ddy)

"""
function calc_dvdtau!(dvdtau::Array{Float64,2}, u::Array{Float64,2}, v::Array{Float64,2}, p::Array{Float64,2}, dx::Float64, dy::Float64, rho::Float64, mu::Float64)
    dvdx = (v[3:end, 2:end-1] - v[1:end-2, 2:end-1]) / (2 * dx)
    dvdy = (v[2:end-1, 3:end] - v[2:end-1, 1:end-2]) / (2 * dy)
    dpdy = (p[2:end-1, 3:end] - p[2:end-1, 1:end-2]) / (2 * dy)
    ddvddx = (v[3:end, 2:end-1] - 2 * v[2:end-1, 2:end-1] + v[1:end-2, 2:end-1]) / (dx^2)
    ddvddy = (v[2:end-1, 3:end] - 2 * v[2:end-1, 2:end-1] + v[2:end-1, 1:end-2]) / (dy^2)

    dvdtau[2:end-1, 2:end-1] = -(u[2:end-1, 2:end-1] .* dvdx .+ v[2:end-1, 2:end-1] .* dvdy) - 1 / rho * dpdy + (mu / rho) * (ddvddx + ddvddy)
end

"""
    calc_dpdtau!

1/beta * dP/dtau + div (u,v) = 0

"""
function calc_dpdtau!(dpdtau::Array{Float64,2}, u::Array{Float64,2}, v::Array{Float64,2}, p::Array{Float64,2}, dx::Float64, dy::Float64, beta::Float64)
    dudx = (u[3:end, 2:end-1] - u[1:end-2, 2:end-1]) / (2 * dx)
    dvdy = (v[2:end-1, 3:end] - v[2:end-1, 1:end-2]) / (2 * dy)
    # println("Size of dudx: ", size(dudx))
    # println("Size of dvdy: ", size(dvdy))
    # println("Max |u|: ", maximum(abs.(u)))
    # println("Max |v|: ", maximum(abs.(v)))
    # println("Max |dudx|: ", maximum(abs.(dudx)))
    # println("Max |dvdy|: ", maximum(abs.(dvdy)))
    div = dudx + dvdy
    dpdtau[2:end-1, 2:end-1] = -div * beta

    # add artificial dissipation
    epsilon = 0.01  # Small dissipation coefficient

    # 4th order dissipation in x-direction
    dissipation_x = epsilon * (p[5:end, 3:end-2] - 4 * p[4:end-1, 3:end-2] + 6 * p[3:end-2, 3:end-2] - 4 * p[2:end-3, 3:end-2] + p[1:end-4, 3:end-2])

    # 4th order dissipation in y-direction  
    dissipation_y = epsilon * (p[3:end-2, 5:end] - 4 * p[3:end-2, 4:end-1] + 6 * p[3:end-2, 3:end-2] - 4 * p[3:end-2, 2:end-3] + p[3:end-2, 1:end-4])

    # Apply to interior points
    # dpdtau[3:end-2, 3:end-2] .-= (dissipation_x + dissipation_y)
end

end