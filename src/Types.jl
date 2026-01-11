module Types

export AbstractVectorField, AbstractGrid
export AbstractVectorField2D, AbstractGrid2D
export AbstractVectorField3D, AbstractGrid3D
export VectorField2D, VectorField3D
export VectorFieldStaggered2D, VectorFieldStaggered3D
export GridRegular2D, Grid2D
export GridRegular3D, Grid3D

export FluidProperties
export RelaxationFactors

export CorrectionCoefficients2D, CorrectionCoefficients3D

################## Abstract Types ##################

abstract type AbstractVectorField end
abstract type AbstractGrid end

abstract type AbstractVectorField2D <: AbstractVectorField end
abstract type AbstractGrid2D <: AbstractGrid end

abstract type AbstractVectorField3D <: AbstractVectorField end
abstract type AbstractGrid3D <: AbstractGrid end

################## Vector Fields ##################
"""
    VectorFieldStaggered2D
Members:
- u
- v
"""
@kwdef mutable struct VectorFieldStaggered2D <: AbstractVectorField2D
    u::Array{Float64,2}
    v::Array{Float64,2}
    function VectorFieldStaggered2D(u::Array{Float64,2}, v::Array{Float64,2})
        new(u, v)
    end
    function VectorFieldStaggered2D(nx::Int64, ny::Int64)
        new(zeros(nx + 1, ny), zeros(nx, ny + 1))
    end
end

"""
    VectorField2D
Members:
- u
- v
"""
@kwdef mutable struct VectorField2D <: AbstractVectorField2D
    u::Array{Float64,2}
    v::Array{Float64,2}
end

"""
    VectorFieldStaggered3D
Members:
- u
- v
- w
"""
@kwdef mutable struct VectorFieldStaggered3D <: AbstractVectorField3D
    u::Array{Float64,3}
    v::Array{Float64,3}
    w::Array{Float64,3}
    function VectorFieldStaggered3D(u::Array{Float64,3}, v::Array{Float64,3}, w::Array{Float64,3})
        new(u, v, w)
    end
    function VectorFieldStaggered3D(nx::Int64, ny::Int64, nz::Int64)
        new(zeros(nx + 1, ny, nz), zeros(nx, ny + 1, nz), zeros(nx, ny, nz + 1))
    end
end

"""
    VectorField3D
Members:
- u
- v
- w
"""
@kwdef mutable struct VectorField3D <: AbstractVectorField3D
    u::Array{Float64,3}
    v::Array{Float64,3}
    w::Array{Float64,3}
end

################## Grid ##################

@kwdef mutable struct GridRegular2D <: AbstractGrid2D
    nx::Int64
    ny::Int64
    x::Array{Float64}
    y::Array{Float64}
    dx::Float64
    dy::Float64
end

@kwdef mutable struct Grid2D <: AbstractGrid2D
    nx::Int64
    ny::Int64
    x::Array{Float64}
    y::Array{Float64}
    dx::Array{Float64}
    dy::Array{Float64}
end

@kwdef mutable struct GridRegular3D <: AbstractGrid3D
    nx::Int64
    ny::Int64
    nz::Int64
    x::Array{Float64}
    y::Array{Float64}
    z::Array{Float64}
    dx::Float64
    dy::Float64
    dz::Float64
end

@kwdef mutable struct Grid3D <: AbstractGrid3D
    nx::Int64
    ny::Int64
    nz::Int64
    x::Array{Float64}
    y::Array{Float64}
    z::Array{Float64}
    dx::Array{Float64}
    dy::Array{Float64}
    dz::Array{Float64}
end

################## uh ##################
"""
Members:
- rho:
- mu:
"""
@kwdef mutable struct FluidProperties
    rho::Float64
    mu::Float64
end

"""
Members:
- velocity:
- pressure:
"""
@kwdef mutable struct RelaxationFactors
    velocity::Float64 = 0.8
    pressure::Float64 = 0.8
end

################## Coefficients ##################
"""
Members:
- d_u: also notated as d_e (for east)
- d_v: also notated as d_n (for north)
"""
@kwdef mutable struct CorrectionCoefficients2D
    d_u::Array{Float64,2}
    d_v::Array{Float64,2}
end

"""
Members:
- d_u: also notated as d_e (for east)
- d_v: also notated as d_n (for north)
- d_w: also notated as d_t (for top)
"""
@kwdef mutable struct CorrectionCoefficients3D
    d_u::Array{Float64,3}
    d_v::Array{Float64,3}
    d_w::Array{Float64,3}
end

end