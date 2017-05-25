#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

macro dispatchhelper(T, args...)
    code = quote end
    ex = Expr(:curly, :Union)
    for arg in args
        push!(code.args, Expr(:typealias, arg, Expr(:curly, :Val, QuoteNode(arg))))
        push!(ex.args, Expr(:curly, :Type, arg))
    end
    push!(code.args, Expr(:typealias, T, ex))
    code
end
@dispatchhelper(ModelType, HereAndNow, WaitAndSee, ExpectedValue)
@dispatchhelper(ModelSense, Minimisation, Maximisation)
@dispatchhelper(RewardType, TerminalReward, InterpolatedReward)
@dispatchhelper(SolveType, Parallel, Serial)

import Base: product, indices, length, convert, rand, isless

immutable NestedCVaRType
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
function NestedCVaR(; lambda=1., beta=0.5)
    if beta < 1e-10
        warn("Beta is very small.")
    end
    NestedCVaRType(lambda, beta)
end
Expectation() = NestedCVaRType(1., 1.)

"""
A storage container for weighted probabilities
"""
immutable WeightedProbability{T}
    value::T
    probability::Float64
end
Base.isless{T}(x::WeightedProbability{T}, y::WeightedProbability{T}) = isless(x.value, y.value)
"""
Sample from a vector of WeightedProbability
"""
function Base.rand{T}(x::Vector{WeightedProbability{T}})
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
function Base.convert{T, V<:AbstractVector}(::Type{Vector{WeightedProbability{T}}}, x::V)
    y = WeightedProbability{T}[]
    prob = 1. / length(x)
    for xi in x
        push!(y, WeightedProbability(xi, prob))
    end
    y
end
Base.convert{T}(::Type{Vector{WeightedProbability{T}}}, x::Vector{WeightedProbability{T}}) = x

"""
Some helper functions to return type
"""
getType{T}(x::AbstractVector{T}) = T
getType{T}(x::AbstractVector{WeightedProbability{T}}) = T

"""
A GenericSpace is a space defined by discrete points along N dimensions.
"""
immutable GenericSpace{T, T2, N}
    dimensions::T                             # list of dimensions
    indices::Tuple{Vararg{UnitRange{Int}, N}} # List of indices
    nameindices::Dict{Symbol, Int}            # a reference naming dict
    minimum::T2 # Minimum values
    maximum::T2 # Minimum values
    bounded::Tuple{Vararg{Bool, N}}
end
GenericSpace(;kwargs...) = GenericSpace(v->v;kwargs...)
GenericSpace(T::DataType; kwargs...) = GenericSpace(v->convert(Vector{T}, v); kwargs...)
GenericSpace(::Type{WeightedProbability}; kwargs...) = GenericSpace(v->convert(Vector{WeightedProbability{getType(v)}}, v); kwargs...)
function GenericSpace(conversion::Function;kwargs...)
    nameindices = Dict{Symbol, Int}()
    dims         = Any[]
    indices      = UnitRange{Int}[]
    mins = Any[]
    maxs = Any[]
    bounded = Bool[]
    j=1
    for (key, value) in kwargs
        nameindices[key] = j
        push!(dims, conversion(value))
        push!(indices, 1:length(value))
        push!(mins, minimum(value))
        push!(maxs, maximum(value))
        push!(bounded, false)
        j+=1
    end
    GenericSpace{typeof(tuple(dims...)), typeof(tuple(mins...)), length(dims)}(tuple(dims...), tuple(indices...), nameindices, tuple(mins...), tuple(maxs...), tuple(bounded...))
end

Base.length(gs::GenericSpace) = length(gs.dimensions)
Base.product(gs::GenericSpace) = product(gs.dimensions...)
Base.product() = [()]
Base.indices(gs::GenericSpace) = product(gs.indices...)

# ===================================
type Stage{T,T2,M,U<:GenericSpace,V<:GenericSpace}
    statespace::GenericSpace{T,T2,M}
    controlspace::U
    noisespace::V
    bellmansurface::Array{Float64, M}
    interpolatedsurface::Interpolations.GriddedInterpolation{Float64, M, Float64, Interpolations.Gridded{Interpolations.Linear}, T,0}
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
    Stage(
    tmpdict[:statespace],
    tmpdict[:controlspace],
    tmpdict[:noisespace],
    Array(Float64, tuple(map(length, tmpdict[:statespace].dimensions)...)),
    interpolate(tuple(fill(Float64[], M)...), Array(Float64, tuple(zeros(Int, M)...)), Gridded(Linear())),
    tmpdict[:dynamics],
    tmpdict[:reward],
    tmpdict[:terminalcost],
    tmpdict[:isfeasible]
    )
end

type SDPModel{T<:Stage, N}
    stages::Tuple{Vararg{T, N}}
    sense::ModelSense
end

function SDPModel(buildstage!::Function;
    stages::Int   = 1,
    sense::Symbol = :Max
    )

    modsense = (sense==:Max?Maximisation:Minimisation)
    outstages = Stage[]
    for t in 1:stages
        tmpdict = Dict{Symbol, Any}(
        :noisespace   => GenericSpace(),
        :terminalcost => (x)->0.,
        :isfeasible => (x,u,w)->true,
        :reward => (x,u,w) -> 0.0
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

type TimingLog
    dynamics::Float64
    constraints::Float64
    reward::Float64
    interpolation::Float64
    paralleloverhead::Float64
    innerloop::Float64
    total::Float64
end
TimingLog() = TimingLog(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
