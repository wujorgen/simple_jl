include("derivatives.jl")
using .Derivatives

using Plots


"""
    main

Solves the Navier Stokes equations for pipe flow, using the artifical compressiblity method. Tau denotes pseudotime.

"""
function main()
    NX::Int64 = 21
    NY::Int64 = 21

    LX::Float64 = 1
    LY::Float64 = 1

    dx::Float64 = LX / (NX - 1)
    dy::Float64 = LY / (NY - 1)

    tau::Float64 = 0.001
    beta::Float64 = 1000

    rho::Float64 = 1
    mu::Float64 = 0.1
    # Re::Float64 = rho * lid_velocity * LY / mu
    # println("Reynolds Number: $Re")

    u::Array{Float64,2} = zeros((NX, NY))
    v::Array{Float64,2} = zeros((NX, NY))
    p::Array{Float64,2} = zeros((NX, NY))

    dudtau::Array{Float64,2} = zeros((NX, NY))
    dvdtau::Array{Float64,2} = zeros((NX, NY))
    dpdtau::Array{Float64,2} = zeros((NX, NY))

    function apply_bcs!(u, v, p)
        # u BCs
        u[:, 1] .= 0.0
        u[:, end] .= 1.0
        u[1, :] .= u[2, :]
        u[end, :] .= u[end-1, :]

        v[:, 1] .= 0.0
        v[:, end] .= 0.0
        v[1, :] .= v[2, :]  # or 0.0
        v[end, :] .= v[end-1, :]  # or 0.0

        # p BCs
        p[1, :] .= p[2, :]
        p[end, :] .= p[end-1, :]
        p[:, 1] .= p[:, 2]
        p[:, end] .= p[:, end-1]
    end

    apply_bcs!(u, v, p)

    for taudx in 1:5000
        calc_dudtau!(dudtau, u, v, p, dx, dy, rho, mu)
        calc_dvdtau!(dvdtau, u, v, p, dx, dy, rho, mu)
        calc_dpdtau!(dpdtau, u, v, p, dx, dy, beta)

        if taudx % 25 == 0
            dudtau_mag = sum(dudtau .* dudtau)
            dvdtau_mag = sum(dvdtau .* dvdtau)
            dpdtau_mag = sum(dpdtau .* dpdtau)
            println("Iteration $taudx:")
            println("\t|dudtau| = $dudtau_mag")
            println("\t|dvdtau| = $dvdtau_mag")
            println("\t|dpdtau| = $dpdtau_mag")
        end

        u .+= dudtau * tau
        v .+= dvdtau * tau
        p .+= dpdtau * tau

        apply_bcs!(u, v, p)
    end

    println("u-velocity at y=0 (bottom): ", u[:, 1])
    println("u-velocity at y=mid: ", u[:, div(NY, 2)])
    println("u-velocity at y=end (top): ", u[:, end])

    xlist = collect(0:dx:LX)
    ylist = collect(0:dy:LY)
    contourf(xlist, ylist, p, color=:viridis)
    savefig("contour.png")
    Plots.CURRENT_PLOT.nullableplot = nothing

    # TODO
    X, Y = [xlist[i] for i in 1:NX, j in 1:NY], [ylist[j] for i in 1:NX, j in 1:NY]
    quiver!(X, Y, quiver=(u, v), color=:black)
    savefig("quiverplot.png")
    Plots.CURRENT_PLOT.nullableplot = nothing
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end