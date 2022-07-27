#  Copyright 2021, Emanuele Concas
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#  Dynamic pricing with customer choice model (MNL)
#   
#  An enterprise needs to determine the price for its inventory, composed of 
#  n substitutable products. There are three main factors that affect an user
#  choice: current inventory x, quality a and price r. The utility for a 
#  product is given by Uᵢ = (aᵢ-rᵢ) + eᵢ, where eᵢ is an independent 
#  identically distributed Gumbel random variables having μ as scale 
#  parameter. The utility of no-buy is given by u₀. Customers arrivals follow
#  a non-homogeneous Poisson process with intensity λ. The probability that a 
#  customer buy product i, with vector of price r, is modelled with a MNL
#  model, and it's denoted by P(r, i). Given the initial inventory, the 
#  enterprise needs to find the best policy to maximize the expected profit 
#  for each period t ∈ T.
#
#  πₜ(xₜ) = λ ∑ᵗₛ₌₁(mₛ(xₜ)-μ) where m is the solution of mₜ(xₜ)P⁰ₜ(rₜ*) = μ
#
#  From 10.1287/msom.1080.0221
#
#  To run with N processors, use
#      julia -p N better_dynamic_pricing.jl

using DynamicProgramming, Distributed, Test

@everywhere begin
    using Roots

    const T = 100
    const inventory = (5, 5, 5)
    const a = [11.75, 9, 6.25]
    const μ = 1
    const λ = 0.1
    const u₀ = 0
end

margins = Dict{Tuple{Int,Tuple},Float64}()

m = SDPModel(stages = T, sense = :Max) do sp, t
    S(x) = findall(x -> x > 0, x)

    function Δπ(t, i, x)
        if t == T
            return 0
        else
            e = [x...]
            e[i] = e[i] - 1
            return m.stages[t+1].interpolatedsurface(x...) -
                   m.stages[t+1].interpolatedsurface(e...)
        end
    end

    @states(sp, begin
        x₁ in 0:inventory[1]
        x₂ in 0:inventory[2]
        x₃ in 0:inventory[3]
    end)

    @controls(sp, begin
        margin = 0.0
    end)

    dynamics!(sp) do y, x, u, w
        y[x₁] = x[x₁]
        y[x₂] = x[x₂]
        y[x₃] = x[x₃]

        return λ * (u[margin] - μ)
    end

    constraints!(sp) do x, u, w
        return all([all([x...] .≤ inventory), sum([x...]) > sum(inventory) - t])
    end

    presolvecallback!(sp) do x, s
        if !haskey(margins, (t, x))
            margins[(t, x)] = find_zero(
                margin::Float64 ->
                    (margin / μ - 1) * exp((margin + u₀) / μ) - sum(
                        exp((a[i] - Δπ(t, i, x)) / μ) for i in S(x);
                        init = 0,
                    ),
                μ,
            )
        end

        return s.controlspace = GenericSpace(; margin = margins[(t, x)])
    end
end

DynamicProgramming.solve(m)

@test isapprox(
    sum(m.stages[t].interpolatedsurface[5, 5, 5] for t in 1:T),
    4145.263,
    atol = 1e-2,
)
