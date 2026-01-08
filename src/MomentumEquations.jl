module MomentumEquations

using ..Types

"""
    solve_u_momentum!

# Arguments

# Returns

"""
function solve_u_momentum!(tentative_velocity::VectorField2D, d_e::Matrix{Float64},
    velocity::VectorField2D, pressure::Matrix{Float64},
    props::FluidProperties, grid::RegularGrid2D;
    max_iter::Int=1000, tol::Float64=1e-6)
    for itr in 1:max_iter
        residual::Float64 = 0.0
        for j in 2:grid.NY-1
            for i in 2:grid.NX-1
                # u-cell CV face velocities
                u_E::Float64 = (velocity.u[i+1, j] + velocity.u[i, j]) / 2
                u_W::Float64 = (velocity.u[i-1, j] + velocity.u[i, j]) / 2
                v_N::Float64 = (velocity.v[i+1, j] + velocity.v[i+1, j+1]) / 2
                v_S::Float64 = (velocity.v[i, j] + velocity.v[i-1, j]) / 2
                # dimensions of current u-cell
                dx_cell::Float64 = grid.dx
                dy_cell::Float64 = grid.dy
                # distances to next u-nodes
                dx_e::Float64 = grid.dx
                dx_w::Float64 = grid.dx
                dy_n::Float64 = grid.dy
                dy_s::Float64 = grid.dy
                #
                a_ij::Float64 = dy_cell * (u_E - u_W) / 2 + dx_cell * (v_N - v_S) / 2 + dy_cell * (props.mu / props.rho) * (1 / dx_e + 1 / dx_w) + dx_cell * (props.mu / props.rho) * (1 / dy_n + 1 / dy_s)
                #if dt != Inf64
                #    # TODO: move transients to their own function call
                #    a_ij += dx_cell * dy_cell / dt
                #end
                a_ip1::Float64 = -dy_cell * u_E / 2 + dy_cell * (props.mu / props.rho) / dx_e
                a_im1::Float64 = dy_cell * u_W / 2 + dy_cell * (props.mu / props.rho) / dx_w
                a_jp1::Float64 = -dx_cell * v_N / 2 + dx_cell * (props.mu / props.rho) / dy_n
                a_jm1::Float64 = dx_cell * v_S / 2 + dx_cell * (props.mu / props.rho) / dy_s
                b::Float64 = A_e * (pressure[i, j] - pressure[i-1, j])
                #if dt != Inf64
                #    # TODO: move transients to their own function call
                #    # b += velocity_old(i, j) * dx_cell * dy_cell / dt
                #end
                A_e::Foat64 = -dy_cell / props.rho
                # TODO: add body forces to b
                #
                u_old = tentative_velocity.u[i, j]
                tentative_velocity.u[i, j]
                d_e[i, j] = A_e / a_ij
                #
                residual += (tentative_velocity.u[i, j] - u_old)^2
            end
        end
        # TODO: check convergence here
        residual /= (grid.nx * grid.ny)

        if residual < tol
            return iter
        end
    end
    @warn "Gauss-Seidel did not converge in $max_iter iterations"
    return max_iter
end

function solve_v_momentum!()
end

function solve_u_momentum!(grid::Grid3D; dt::Float64=Inf64)
end

function solve_v_momentum!(grid::Grid3D; dt::Float64=Inf64)
end

function solve_w_momentum!(grid::Grid3D; dt::Float64=Inf64)
end

end