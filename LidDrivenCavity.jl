include("src/Types.jl")
include("src/MomentumEquations.jl")
include("src/Corrections.jl")
include("src/BoundaryConditions.jl")
using .Types
using .MomentumEquations
using .Corrections
using .BoundaryConditions
using Plots


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


function main()
    NX::Int64 = 51
    NY::Int64 = 51

    LX::Float64 = 1
    LY::Float64 = 1

    dx::Float64 = LX / (NX - 1)
    dy::Float64 = LY / (NY - 1)

    grid = GridRegular2D(
        nx=NX, ny=NY,
        dx=dx, dy=dy,
        x=collect(0:dx:LX),
        y=collect(0:dy:LY),
    )

    properties = FluidProperties(rho=1, mu=0.01)

    rlx = RelaxationFactors(velocity=0.3, pressure=0.3)

    velocity_star = VectorFieldStaggered2D(NX, NY)

    velocity = VectorFieldStaggered2D(NX, NY)

    ccfs = CorrectionCoefficients2D(
        d_u=zeros(Float64, NX + 1, NY),
        d_v=zeros(Float64, NX, NY + 1)
    )

    pressure = zeros(Float64, NX, NY)
    p_corr = zeros(Float64, NX, NY)

    LID_VELOCITY = (1.0, 0.0)
    set_velocity_boundaries!(velocity, top=LID_VELOCITY)
    set_velocity_boundaries!(velocity_star, top=LID_VELOCITY)

    residual_history = []
    continuity_history = []
    for otr in 1:500
        velocity_star.u .= velocity.u
        velocity_star.v .= velocity.v

        solve_u_momentum!(velocity_star, ccfs, velocity, pressure, properties, grid)
        solve_v_momentum!(velocity_star, ccfs, velocity, pressure, properties, grid)
        set_velocity_boundaries!(velocity_star, top=LID_VELOCITY)

        p_corr .= 0
        solve_pressure_correction!(p_corr, velocity_star, ccfs, properties, grid)

        continuity_error = compute_continuity_error(velocity_star, grid)
        push!(continuity_history, continuity_error)
        println("Continuity error: $continuity_error")

        correct_pressure!(pressure, p_corr, rlx, grid)

        push!(residual_history, sum(p_corr .* p_corr)^(0.5))

        correct_velocity!(velocity, velocity_star, p_corr, ccfs, grid, rlx)
        set_velocity_boundaries!(velocity, top=LID_VELOCITY)
        set_pressure_boundaries!(pressure)
    end

    plot(residual_history,
        yscale=:log10,
        xlabel="Iteration",
        ylabel="Residual",
        title="SIMPLE Convergence",
        marker=:circle,
        legend=false)
    savefig("convergence.png")

    Plots.CURRENT_PLOT.nullableplot = nothing

    X, Y = [grid.x[i] for i in 1:grid.nx, j in 1:grid.ny], [grid.y[j] for i in 1:grid.nx, j in 1:grid.ny]

    u = (velocity.u[1:grid.nx, :] .+ velocity.u[2:grid.nx+1, :]) ./ 2

    v = (velocity.v[:, 1:grid.ny] .+ velocity.v[:, 2:grid.ny+1]) ./ 2

    u = u ./ 100
    v = v ./ 100

    contourf(grid.x, grid.y, pressure', color=:viridis)
    quiver!(X', Y', quiver=(u', v'), color=:black)
    savefig("quiverplot.png")
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end