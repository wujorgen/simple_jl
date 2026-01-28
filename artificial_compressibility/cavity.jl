include("derivatives.jl")
using .Derivatives

using Plots


"""
    main

Solves the Navier Stokes equations for lid driven cavity, using the artifical compressiblity method. Tau denotes pseudotime.

"""
function main()
    NX::Int64 = 41
    NY::Int64 = 41

    LX::Float64 = 1
    LY::Float64 = 1

    dx::Float64 = LX / (NX - 1)
    dy::Float64 = LY / (NY - 1)

    tau::Float64 = 0.0001
    beta::Float64 = 100

    lid_velocity::Float64 = 1.0
    rho::Float64 = 1
    mu::Float64 = 1
    Re::Float64 = rho * lid_velocity * LY / mu
    println("Reynolds Number: $Re")

    u::Array{Float64,2} = zeros((NX, NY))
    v::Array{Float64,2} = zeros((NX, NY))
    p::Array{Float64,2} = zeros((NX, NY))

    dudtau::Array{Float64,2} = zeros((NX, NY))
    dvdtau::Array{Float64,2} = zeros((NX, NY))
    dpdtau::Array{Float64,2} = zeros((NX, NY))

    function apply_bcs!(u, v, p)
        # No-slip on all walls
        u[:, 1] .= 0.0      # bottom
        u[:, end] .= lid_velocity  # top (lid)
        u[1, :] .= 0.0      # left
        u[end, :] .= 0.0    # right

        v[:, 1] .= 0.0      # bottom
        v[:, end] .= 0.0    # top
        v[1, :] .= 0.0      # left
        v[end, :] .= 0.0    # right

        # Neumann BCs for pressure (zero gradient)
        p[1, :] .= p[2, :]
        p[end, :] .= p[end-1, :]
        p[:, 1] .= p[:, 2]
        p[:, end] .= p[:, end-1]
    end

    calc_dudtau!(dudtau, u, v, p, dx, dy, rho, mu)
    calc_dvdtau!(dvdtau, u, v, p, dx, dy, rho, mu)
    calc_dpdtau!(dpdtau, u, v, dx, dy, beta)

    # println(dudtau * tau)

    for taudx in 1:100
        # Gradual lid startup
        lid_vel_current = lid_velocity * min(1.0, taudx * tau / 0.1)  # Ramp over 0.1 time units
        println("lid velocity setting is $lid_vel_current")

        # apply boundary conditions with ramped velocity
        u[:, 1] .= 0.0
        u[:, end] .= lid_vel_current  # Use ramped velocity
        u[1, :] .= 0.0
        u[end, :] .= 0.0

        v[:, 1] .= 0.0
        v[:, end] .= 0.0
        v[1, :] .= 0.0
        v[end, :] .= 0.0

        p[1, :] .= p[2, :]
        p[end, :] .= p[end-1, :]
        p[:, 1] .= p[:, 2]
        p[:, end] .= p[:, end-1]

        # calculate derivatives
        calc_dudtau!(dudtau, u, v, p, dx, dy, rho, mu)
        calc_dvdtau!(dvdtau, u, v, p, dx, dy, rho, mu)
        calc_dpdtau!(dpdtau, u, v, dx, dy, beta)
        if taudx == 1 & false
            i, j = 10, 10
            println("At point ($i, $j):")
            println("  u = ", u[i, j])
            println("  v = ", v[i, j])
            println("  p = ", p[i, j])
            println("  dudtau = ", dudtau[i, j])
            println("  dvdtau = ", dvdtau[i, j])
            println("  dpdtau = ", dpdtau[i, j])

            # Manually compute what dudtau should be at this point
            u_local = u[i-1:i+1, j-1:j+1]
            println("  Local u neighborhood:")
            println(u_local)

            manual_dudx = (u[i+1, j] - u[i-1, j]) / (2 * dx)
            manual_dudy = (u[i, j+1] - u[i, j-1]) / (2 * dy)
            manual_dpdx = (p[i+1, j] - p[i-1, j]) / (2 * dx)
            manual_laplacian = ((u[i+1, j] - 2 * u[i, j] + u[i-1, j]) / dx^2 +
                                (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) / dy^2)

            manual_dudtau = -(u[i, j] * manual_dudx + v[i, j] * manual_dudy) - (1 / rho) * manual_dpdx + (mu / rho) * manual_laplacian

            println("  Manual dudtau calculation = ", manual_dudtau)
            println("  Function dudtau = ", dudtau[i, j])
        end
        # track correction magnitude
        dpdtau_mag = sum(dpdtau .* dpdtau)
        println("Iteration $taudx, |dpdtau| = $dpdtau_mag")
        #println("Sum(u): ", sum(u))
        #println("Sum(v): ", sum(v))
        # integrate
        u .+= dudtau * tau
        v .+= dvdtau * tau
        p .+= dpdtau * tau
    end

end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end