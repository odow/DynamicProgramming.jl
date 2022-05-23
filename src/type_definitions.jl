#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

abstract type ModelType end
abstract type HereAndNow <: ModelType end
abstract type WaitAndSee <: ModelType end
abstract type ExpectedValue <: ModelType end

abstract type ModelSense end
abstract type Minimisation <: ModelSense end
abstract type Maximisation <: ModelSense end

abstract type RewardType end
abstract type TerminalReward <: RewardType end
abstract type InterpolatedReward <: RewardType end

abstract type SolveType end
abstract type Parallel <: SolveType end
abstract type Serial <: SolveType end

import Base: product, axes, length, convert, rand, isless

struct NestedCVaRType
    lambda::Float64
    beta::Float64
end

"""
    NestedCVaR(;lambda=0.0, beta=0.0)

# Description

A risk measure that is a convex combination of Expectation and Conditional Value @ Risk
(also called Average Value @ Risk).

    λ * E[x] + (1 - λ) * CV@R(1-β)[x]

# Keyword Arguments

 * `lambda`
Convex weight on the expectation (`(1-lambda)` weight is put on the CV@R component.
Inreasing values of `lambda` are less risk averse (more weight on expecattion)

 * `beta`
 The quantile at which to calculate the Conditional Value @ Risk. Increasing values
 of `beta` are less risk averse. If `beta=0`, then the CV@R component is the
 worst case risk measure.

# Returns

A risk measure to use
"""
function NestedCVaR(; lambda = 1.0, beta = 0.5)
    if beta < 1e-10
        warn("Beta is very small.")
    end
    return NestedCVaRType(lambda, beta)
end
Expectation() = NestedCVaRType(1.0, 1.0)

"""
A storage container for weighted probabilities
"""
struct WeightedProbability{T}
    value::T
    probability::Float64
end
function Base.isless(
    x::WeightedProbability{T},
    y::WeightedProbability{T},
) where {T}
    return isless(x.value, y.value)
end

"""
Sample from a vector of WeightedProbability
"""
function Base.rand(x::Vector{WeightedProbability{T}}) where {T}
    r = rand()
    @inbounds for i in x
        r -= i.probability
        if r <= 0
            return i.value
        end
    end
    return x[end].value
end

"""
Conversion helpers for Vectors of WeightedProbability
"""
function Base.convert(
    ::Type{Vector{WeightedProbability{T}}},
    x::V,
) where {T,V<:AbstractVector}
    y = WeightedProbability{T}[]
    prob = 1.0 / length(x)
    for xi in x
        push!(y, WeightedProbability(xi, prob))
    end
    return y
end
function Base.convert(
    ::Type{Vector{WeightedProbability{T}}},
    x::Vector{WeightedProbability{T}},
) where {T}
    return x
end

"""
Some helper functions to return type
"""
getType(x::AbstractVector{T}) where {T} = T
getType(x::AbstractVector{WeightedProbability{T}}) where {T} = T

"""
A GenericSpace is a space defined by discrete points along N dimensions.
"""
struct GenericSpace{T,T2,N}
    dimensions::T                             # list of dimensions
    indices::Tuple{Vararg{UnitRange{Int},N}} # List of indices
    nameindices::Dict{Symbol,Int}            # a reference naming dict
    minimum::T2 # Minimum values
    maximum::T2 # Minimum values
    bounded::Tuple{Vararg{Bool,N}}
end

GenericSpace(; kwargs...) = GenericSpace(v -> v; kwargs...)
function GenericSpace(T::DataType; kwargs...)
    return GenericSpace(v -> convert(Vector{T}, v); kwargs...)
end
function GenericSpace(::Type{WeightedProbability}; kwargs...)
    return GenericSpace(
        v -> convert(Vector{WeightedProbability{getType(v)}}, v);
        kwargs...,
    )
end
function GenericSpace(conversion::Function; kwargs...)
    nameindices = Dict{Symbol,Int}()
    dims = Any[]
    indices = UnitRange{Int}[]
    mins = Any[]
    maxs = Any[]
    bounded = Bool[]

    for (j, (key, value)) in enumerate(kwargs)
        nameindices[key] = j
        push!(dims, conversion(value))
        push!(indices, 1:length(value))
        push!(mins, minimum(value))
        push!(maxs, maximum(value))
        # always true, we are considering only discretized states
        push!(bounded, true)
    end

    return GenericSpace{
        typeof(tuple(dims...)),
        typeof(tuple(mins...)),
        length(dims),
    }(
        tuple(dims...),
        tuple(indices...),
        nameindices,
        tuple(mins...),
        tuple(maxs...),
        tuple(bounded...),
    )
end

Base.length(gs::GenericSpace) = length(gs.dimensions)
Base.product(gs::GenericSpace) = product(gs.dimensions...)
Base.product() = [()]
Base.axes(gs::GenericSpace) = product(gs.indices...)

# ===================================
mutable struct Stage{T,T2,M,U<:GenericSpace,V<:GenericSpace}
    statespace::GenericSpace{T,T2,M}
    controlspace::U
    noisespace::V
    bellmansurface::Array{Float64,M}
    interpolatedsurface::Interpolations.GriddedInterpolation #{Float64, M, Float64, Interpolations.Gridded{Interpolations.Linear}, T, 0}
    dynamics!::Function                 # dynamics!(x', x, u, w)
    reward::Function                    # r(x, u, w)
    terminalcost::Function              # k(x, u, w)
    isfeasible::Function                # k(x, u, w)
end

function Stage(M, tmpdict)
    if !haskey(tmpdict, :dynamics)
        error("Must specify dynamics")
    end

    if !haskey(tmpdict, :reward)
        error("Must specify reward function")
    end

    return Stage(
        tmpdict[:statespace],
        tmpdict[:controlspace],
        tmpdict[:noisespace],
        Array{Float64}(
            undef,
            tuple(map(length, tmpdict[:statespace].dimensions)...),
        ),
        interpolate(
            tuple(fill(Float64[], M)...),
            Array{Float64}(undef, tuple(zeros(Int, M)...)),
            Gridded(Linear()),
        ),
        tmpdict[:dynamics],
        tmpdict[:reward],
        tmpdict[:terminalcost],
        tmpdict[:isfeasible],
    )
end

mutable struct SDPModel{T<:Stage,N}
    stages::Tuple{Vararg{T,N}}
    sense::Type{<:ModelSense}
end

function SDPModel(buildstage!::Function; stages::Int = 1, sense::Symbol = :Max)
    modsense = (sense == :Max ? Maximisation : Minimisation)
    outstages = Stage[]
    for t in 1:stages
        tmpdict = Dict{Symbol,Any}(
            :noisespace => GenericSpace(),
            :terminalcost => (x) -> 0.0,
            :isfeasible => (x, u, w) -> true,
            :reward => (x, u, w) -> 0.0,
        )

        buildstage!(tmpdict, t)

        # Compilation trigger
        x = tuple([dim[1] for dim in tmpdict[:statespace].dimensions]...)
        y = [dim[1] for dim in tmpdict[:statespace].dimensions]
        u = tuple([dim[1] for dim in tmpdict[:controlspace].dimensions]...)
        w = tuple([dim[1].value for dim in tmpdict[:noisespace].dimensions]...)
        tmpdict[:terminalcost](x)
        tmpdict[:isfeasible](x, u, w)
        tmpdict[:dynamics](y, x, u, w)
        tmpdict[:reward](x, u, w)
        # ====================

        sp = Stage(length(tmpdict[:statespace]), tmpdict)
        push!(outstages, sp)
    end
    return SDPModel(tuple(outstages...), modsense)
end

mutable struct TimingLog
    dynamics::Float64
    constraints::Float64
    reward::Float64
    interpolation::Float64
    paralleloverhead::Float64
    innerloop::Float64
    total::Float64
end
TimingLog() = TimingLog(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
