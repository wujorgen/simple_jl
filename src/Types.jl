module Types

export VectorField
export Grid

export VectorField2D, VectorField3D
export Mesh2D, Mesh3D
export FluidProperties

abstract type AbstractVectorField end
abstract type AbstractGrid end

"""
Members:
- u
- v
"""
@kwdef mutable struct VectorField2D <: AbstractVectorField
    u::Matrix{Float64}
    v::Matrix{Float64}
end

@kwdef mutable struct VectorField3D <: AbstractVectorField
    u::Matrix{Float64}
    v::Matrix{Float64}
    w::Matrix{Float64}
end

@kwdef mutable struct RegularGrid2D <: AbstractGrid
    nx::Float64
    ny::Float64
    x::Array{Float64}
    y::Array{Float64}
    dx::Float64
    dy::Float64
end

@kwdef mutable struct Grid2D <: AbstractGrid
    nx::Float64
    ny::Float64
    x::Array{Float64}
    y::Array{Float64}
    dx::Array{Float64}
    dy::Array{Float64}
end

@kwdef mutable struct Grid3D <: AbstractGrid
    nx::Float64
    ny::Float64
    nz::Float64
    x::Array{Float64}
    y::Array{Float64}
    z::Array{Float64}
    dx::Array{Float64}
    dy::Array{Float64}
    dz::Array{Float64}
end

@kwdef mutable struct FluidProperties
    rho::Float64
    mu::Float64
end

end