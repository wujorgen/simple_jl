module Utilities

export compute_continuity_error

using ..Types

function compute_continuity_error(velocity::VectorFieldStaggered2D, grid::GridRegular2D)
    max_error = 0.0
    for j in 2:grid.ny-1
        for i in 2:grid.nx-1
            div = (velocity.u[i+1, j] - velocity.u[i, j]) / grid.dx +
                  (velocity.v[i, j+1] - velocity.v[i, j]) / grid.dy
            max_error = max(max_error, abs(div))
        end
    end
    return max_error
end

end