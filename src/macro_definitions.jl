#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

struct ModelVariable
    idx::Int
end

import Base: getindex

Base.getindex(x::AbstractVector, i::ModelVariable) = x[i.idx]
Base.getindex(x::Tuple, i::ModelVariable) = x[i.idx]
Base.getindex(x::AbstractVector, i::Tuple{ModelVariable}) = x[i...]

Base.getindex(x::AbstractVector{WeightedProbability}, i::ModelVariable) = x[i.idx].value
Base.getindex(x::Tuple{Vararg{WeightedProbability, N}}, i::ModelVariable) where {N} = x[i.idx].value
Base.getindex(x::AbstractVector{WeightedProbability}, i::Tuple{ModelVariable}) = x[i...].value
Base.getindex(x::Tuple, i::Tuple{ModelVariable}) = x[i...]
Base.setindex!(y::AbstractVector, v::T, i::ModelVariable) where {T<:Number} = (y[i.idx] = v)
Base.setindex!(y::AbstractVector, v::T, i::Tuple{ModelVariable}) where {T<:Number} = (y[i] = v)

function macrobody!(ex, blk)
    code = quote end
    i = 1
    filter!(arg -> typeof(arg) != LineNumberNode, blk.args)
    for line in blk.args
        if !Base.Meta.isexpr(line, :line)
            if line.head == :call
                @assert line.args[1] == :in || line.args[1] == :(∈)
                varname = line.args[2]
                points  = line.args[3]
            elseif line.head == :(=)
                varname = line.args[1]
                points  = line.args[2]
            else
                error("Unvalid syntax in $(line)")
            end
            push!(ex.args, Expr(:kw, varname, esc(points)))
            push!(code.args, esc(Expr(:(=), varname, Expr(:call, Expr(:(.), :DynamicProgramming, QuoteNode(:ModelVariable)), i))))
            i += 1
        end
    end
    push!(code.args, ex)
    code
end

macro states(sp, blk)
    sp = esc(sp)
    ex = Expr(:call, :addstates!, sp)
    macrobody!(ex, blk)
end

macro controls(sp, blk)
    sp = esc(sp)
    ex = Expr(:call, :addcontrols!, sp)
    macrobody!(ex, blk)
end

macro noises(sp, blk)
    sp = esc(sp)
    ex = Expr(:call, :addnoises!, sp)
    macrobody!(ex, blk)
end
