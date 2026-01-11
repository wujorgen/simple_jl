module Corrections

export solve_pressure_correction!
export correct_pressure!
export correct_velocity!

using ..Types

"""
    solve_pressure_correction!

# Arguments
- `p_corr::Array{Float64, 2}`: modified in place
"""
function solve_pressure_correction!(p_corr::Array{Float64,2}, velocity::VectorFieldStaggered2D, ccfs::CorrectionCoefficients2D, props::FluidProperties, grid::GridRegular2D; max_itr::Int64=1000, tol::Float64=1e-9)
    for itr in 1:max_itr
        residual::Float64 = 0.0
        for j in 2:grid.ny-1
            for i in 2:grid.nx-1
                #
                u_E::Float64 = velocity.u[i+1, j]
                u_W::Float64 = velocity.u[i, j]
                v_N::Float64 = velocity.v[i, j+1]
                v_S::Float64 = velocity.v[i, j]
                #
                a_ip1::Float64 = -props.rho * ccfs.d_u[i+1, j] * grid.dy
                a_im1::Float64 = -props.rho * ccfs.d_u[i, j] * grid.dy
                a_jp1::Float64 = -props.rho * ccfs.d_v[i, j+1] * grid.dx
                a_jm1::Float64 = -props.rho * ccfs.d_v[i, j] * grid.dx
                a_ij::Float64 = a_ip1 + a_im1 + a_jp1 + a_jm1
                b::Float64 = grid.dy * props.rho * (u_W - u_E) + grid.dx * props.rho * (v_S - v_N)
                #
                p_corr_old = p_corr[i, j]
                p_corr[i, j] = (a_ip1 * p_corr[i+1, j] + a_im1 * p_corr[i-1, j] + a_jp1 * p_corr[i, j+1] + a_jm1 * p_corr[i, j-1] + b) / a_ij
                residual += (p_corr[i, j] - p_corr_old)^2
            end
        end
        residual /= (grid.nx * grid.ny)
        if residual < tol
            return itr
        end
    end
    @warn "Pressure correction not converged."
    return max_itr
end

function correct_pressure!(pressure::Array{Float64,2}, p_corr::Array{Float64,2}, relax::RelaxationFactors, grid::GridRegular2D)
    for j in 2:grid.ny-1
        for i in 2:grid.nx-1
            pressure[i, j] += relax.pressure * p_corr[i, j]
        end
    end
end

"""
    correct_velocity!
Members
- `velocity::VectorFieldStaggered2D`: modified in place
- `tentative_velocity::VectorFieldStaggered2D`
- `p_corr::Array{Float64, 2}`
- `ccfs::CorrectionCoefficients2D`
- `relax::RelaxationFactors`
"""
function correct_velocity!(velocity::VectorFieldStaggered2D, tentative_velocity::VectorFieldStaggered2D, p_corr::Array{Float64,2}, ccfs::CorrectionCoefficients2D, grid::GridRegular2D, relax::RelaxationFactors)
    for j in 2:grid.ny-1
        for i in 2:grid.nx
            velocity.u[i, j] = tentative_velocity.u[i, j] + relax.velocity * ccfs.d_u[i, j] * (p_corr[i, j] - p_corr[i-1, j])
        end
    end

    for j in 2:grid.ny
        for i in 2:grid.nx-1
            velocity.v[i, j] = tentative_velocity.v[i, j] + relax.velocity * ccfs.d_v[i, j] * (p_corr[i, j] - p_corr[i, j-1])
        end
    end
end

end