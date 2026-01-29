include("derivatives.jl")
using .Derivatives

using Plots


"""
    main

Solves the Navier Stokes equations for pipe flow, using the artifical compressiblity method. Tau denotes pseudotime.

"""
function main()
    NX::Int64 = 61
    NY::Int64 = 41

    LX::Float64 = 1
    LY::Float64 = 0.4

    dx::Float64 = LX / (NX - 1)
    dy::Float64 = LY / (NY - 1)

    tau::Float64 = 0.00001
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
        u[:, 1] .= 0.0  # bottom
        u[:, end] .= 0.0  # top
        ystart::Int64 = 2
        yend::Int64 = div(NY, 2) - 1
        ymid::Int64 = div(ystart + yend, 2)
        for ydx in ystart:yend
            u[1, ydx] = 0.1 * (1 - (abs(ydx - ymid) / ymid)^2)
        end
        # println(u[1, :])
        # u[1, 5:end-4] .= 1.0  # left
        u[end, :] .= u[end-1, :]  # right

        v[:, 1] .= 0.0  # bottom
        v[:, end] .= 0.0  # top
        v[1, :] .= 0.0  # left
        v[end, :] .= v[end-1, :]  # right

        # p BCs
        p[1, :] .= p[2, :]  # left
        p[end, :] .= p[end-1, :]  # right
        p[:, 1] .= p[:, 2]  # bottom
        p[:, end] .= p[:, end-1]  # top

        # block!
        for idx in 1:div(NX, 3)
            for jdx in yend:NY
                u[idx, jdx] = 0.0
                v[idx, jdx] = 0.0
            end
        end
        p[1:div(NX, 3), yend] = p[1:div(NX, 3), yend-1]
        p[div(NX, 3), yend:NY] = p[div(NX, 3)+1, yend:NY]
    end

    apply_bcs!(u, v, p)

    for taudx in 1:20000
        calc_dudtau!(dudtau, u, v, p, dx, dy, rho, mu)
        calc_dvdtau!(dvdtau, u, v, p, dx, dy, rho, mu)
        calc_dpdtau!(dpdtau, u, v, p, dx, dy, beta)

        if taudx % 500 == 0
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

    #println("u-velocity at y=0 (bottom): ", u[:, 1])
    #println("u-velocity at y=mid: ", u[:, div(NY, 2)])
    #println("u-velocity at y=end (top): ", u[:, end])

    xlist = collect(0:dx:LX)
    ylist = collect(0:dy:LY)
    X, Y = [xlist[i] for i in 1:NX, j in 1:NY], [ylist[j] for i in 1:NX, j in 1:NY]

    contourf(xlist, ylist, p', color=:viridis)
    savefig("pipe_contour.png")
    Plots.CURRENT_PLOT.nullableplot = nothing

    # TODO
    quiver!(X[1:2:end, 1:2:end], Y[1:2:end, 1:2:end], quiver=(u[1:2:end, 1:2:end], v[1:2:end, 1:2:end]), arrow=arrow(:auto, :head, 0.1, 0.1), color=:black)  # arrow(headlength=10, headwidth=5)  # arrow(:auto, :head, 0.1, 0.1)
    savefig("pipe_quiverplot.png")
    Plots.CURRENT_PLOT.nullableplot = nothing
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end