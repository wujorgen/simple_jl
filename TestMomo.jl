include("src/Types.jl")
include("src/MomentumEquations.jl")
using .Types
using .MomentumEquations

NX::Int64 = 11
NY::Int64 = 11

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

velocity_star = VectorField2D(
    u=zeros(Float64, NX + 1, NY),
    v=zeros(Float64, NX, NY + 1)
)

velocity = VectorField2D(
    u=zeros(Float64, NX + 1, NY),
    v=zeros(Float64, NX, NY + 1)
)

ccfs = CorrectionCoefficients2D(
    d_u=zeros(Float64, NX + 1, NY),
    d_v=zeros(Float64, NX, NY + 1)
)

velocity.u[:,end] .= 2.0
velocity_star.u[:,end] .= 2.0

println(velocity.u)

pressure = zeros(Float64, NX, NY)

@time u_itr::Int64 = solve_u_momentum!(velocity_star, ccfs, velocity, pressure, properties, grid)

println(u_itr)

println(velocity_star.u[:, end-1])