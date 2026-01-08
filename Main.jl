include("src/Types.jl")
using .Types

NX::Int64 = 11
NY::Int64 = 11

VelocityField = VectorField2D(
    zeros(NX + 1, NY),
    zeros(NX, NY + 1)
)