module MomentumEquations

export solve_u_momentum!
export solve_v_momentum!

using ..Types

"""
    solve_u_momentum!
SIMPLE u-velocity predictor for steady-state simulations
# Arguments

# Returns

"""
function solve_u_momentum!(tentative_velocity::VectorField2D, ccfs::CorrectionCoefficients2D,
    velocity::VectorField2D, pressure::Array{Float64,2},
    props::FluidProperties, grid::GridRegular2D;
    max_itr::Int=10, tol::Float64=1e-9)::Int64
    for itr in 1:max_itr
        residual::Float64 = 0.0
        for j in 2:grid.ny-1
            for i in 2:grid.nx-1
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
                # coefficients
                a_ij::Float64 = dy_cell * (u_E - u_W) / 2 + dx_cell * (v_N - v_S) / 2 + dy_cell * (props.mu / props.rho) * (1 / dx_e + 1 / dx_w) + dx_cell * (props.mu / props.rho) * (1 / dy_n + 1 / dy_s)
                #if dt != Inf64
                #    # TODO: move transient terms to their own function call - use multiple dispatch!
                #    a_ij += dx_cell * dy_cell / dt
                #end
                a_ip1::Float64 = -dy_cell * u_E / 2 + dy_cell * (props.mu / props.rho) / dx_e
                a_im1::Float64 = dy_cell * u_W / 2 + dy_cell * (props.mu / props.rho) / dx_w
                a_jp1::Float64 = -dx_cell * v_N / 2 + dx_cell * (props.mu / props.rho) / dy_n
                a_jm1::Float64 = dx_cell * v_S / 2 + dx_cell * (props.mu / props.rho) / dy_s
                A_e::Float64 = -dy_cell / props.rho
                b::Float64 = A_e * (pressure[i, j] - pressure[i-1, j])
                #if dt != Inf64
                #    # TODO: move transients to their own function call
                #    # b += velocity_old(i, j) * dx_cell * dy_cell / dt
                #end
                # TODO: add body forces to b
                #
                u_old = tentative_velocity.u[i, j]
                tentative_velocity.u[i, j] = (a_ip1 * tentative_velocity.u[i+1, j] + a_im1 * tentative_velocity.u[i-1, j] + a_jp1 * tentative_velocity.u[i, j+1] + a_jm1 * tentative_velocity.u[i, j-1] + b) / a_ij
                ccfs.d_u[i, j] = A_e / a_ij
                #
                residual += (tentative_velocity.u[i, j] - u_old)^2
            end
        end
        # TODO: check convergence here
        residual /= (grid.nx * grid.ny)

        if residual < tol
            return itr
        end
    end
    # @warn "Gauss-Seidel did not converge in $max_iter iterations"  # because this is an inner loop within SIMPLE, it doesn't really need to fully converge
    return max_itr
end

function solve_v_momentum!(tentative_velocity::VectorField2D, ccfs::CorrectionCoefficients2D,
    velocity::VectorField2D, pressure::Array{Float64,2},
    props::FluidProperties, grid::GridRegular2D;
    max_itr::Int=10, tol::Float64=1e-9)::Int64
    for itr in 1:max_itr
        residual::Float64 = 0.0
        for j in 2:grid.ny-1
            for i in 2:grid.nx-1
                # v-cell CV face velocities
                # dimensions of current v-cell
                # distances to next v-nodes
                # coefficients
            end
        end
        residual /= (grid.nx * grid.ny)

        if residual < tol
            return itr
        end
    end
    return max_itr
end

function solve_u_momentum!(tentative_velocity::VectorField3D, ccfs::CorrectionCoefficients3D,
    velocity::VectorField3D, pressure::Array{Float64,3},
    props::FluidProperties, grid::GridRegular3D;
    max_itr::Int=10, tol::Float64=1e-9)::Int64
end

function solve_v_momentum!(tentative_velocity::VectorField3D, ccfs::CorrectionCoefficients3D,
    velocity::VectorField3D, pressure::Array{Float64,3},
    props::FluidProperties, grid::GridRegular3D;
    max_itr::Int=10, tol::Float64=1e-9)::Int64
end

function solve_w_momentum!(tentative_velocity::VectorField3D, ccfs::CorrectionCoefficients3D,
    velocity::VectorField3D, pressure::Array{Float64,3},
    props::FluidProperties, grid::GridRegular3D;
    max_itr::Int=10, tol::Float64=1e-9)::Int64
end

end