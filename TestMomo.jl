include("src/Types.jl")
include("src/MomentumEquations.jl")
include("src/Corrections.jl")
using .Types
using .MomentumEquations
using .Corrections

NX::Int64 = 5
NY::Int64 = 5

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

properties = FluidProperties(10, 1)

velocity_star = VectorFieldStaggered2D(NX, NY)

velocity = VectorFieldStaggered2D(NX, NY)

ccfs = CorrectionCoefficients2D(
    d_u=zeros(Float64, NX + 1, NY),
    d_v=zeros(Float64, NX, NY + 1)
)

velocity.u[:, end] .= 2.0
velocity_star.u[:, end] .= 2.0

println(velocity.u)

pressure = zeros(Float64, NX, NY)
p_corr = zeros(Float64, NX, NY)

@time u_itr::Int64 = solve_u_momentum!(velocity_star, ccfs, velocity, pressure, properties, grid)
@time v_itr::Int64 = solve_v_momentum!(velocity_star, ccfs, velocity, pressure, properties, grid)

println(u_itr)

println(velocity_star.u)

println("u_star")
for row in eachrow(velocity_star.u)
    println(row)
end

println("v_star")
for row in eachrow(velocity_star.v)
    println(row)
end

@time p_itr::Int64 = solve_pressure_correction!(p_corr, velocity_star, ccfs, properties, grid)
for row in eachrow(p_corr)
    println(row)
end
