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
#  πₜ(x) = max rₜ∈ℝⁿ {λ ∑i∈S(x) P(rₜ, i)[rₜ(i) Δπ(t, i, x)] + πₜ₋₁(x) }
#  where Δπ = πₜ₋₁(x) - πₜ₋₁(x-eⁱ), eⁱ is the unit vector with 1 in
#  position i and S(x) returns the non zero inventory items.
#
#  To run with N processors, use
#      julia -p N dynamic_pricing.jl

using DynamicProgramming
using Test

T = 4

m = SDPModel(stages = T, sense = :Max) do sp, t
    inventory = (2, 1)
    a = (10, 5)
    μ = 1
    λ = 0.1
    u₀ = 0

    S(x) = findall(x -> x > 0, x)
    function P(r, i)
        return exp((a[i] - r[i]) / μ) / (
            sum(exp((a[j] - r[j]) / μ) for j in eachindex(a); init = 0.0) +
            exp(u₀ / μ)
        )
    end

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
    end)

    @controls(sp, begin
        r = collect(Base.product(0:1:5, 0:1:5))
    end)

    dynamics!(sp) do y, x, u, w
        y[x₁] = x[x₁]
        y[x₂] = x[x₂]

        return sum(
            λ * P(u[r], i) * (u[r][i] - Δπ(t, i, x)) for i in S(x);
            init = 0.0,
        )
    end

    constraints!(sp) do x, u, w
        return all([all([x...] .≤ inventory), sum([x...]) > sum(inventory) - t])
    end
end

solve(m)

@test isapprox(m.stages[1].interpolatedsurface[2.0, 1.0], 1.96, atol = 1e-2)
