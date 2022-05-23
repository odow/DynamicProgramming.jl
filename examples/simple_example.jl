#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
#
#   To run with N processors, use
#       julia -p N battery_storage.jl
#
using DynamicProgramming, Distributed

@everywhere begin
    using Distributions

    const sampled_errors = rand(Normal(0, 0.0785220888754297), 20)
    const σ² = range(1, 0, 28) # Decreasing variance in changes in price over time
    const transaction_cost = 0.01
end

m = SDPModel(stages = 28, sense = :Max) do sp, t
    @states(sp, begin
        contracts in range(0, stop=1.5, length=16)
        price in range(3, stop=9, length=10)
        production in range(0, stop=1.5, length=16)
    end)

    @controls(sp, begin
        buy in range(0, stop=1.2, length=13)
    end)

    @noises(sp, begin
        alpha in sampled_errors
        beta in range(0.0, stop=0.05, length=5)
    end)

    dynamics!(sp) do y, x, u, w
        # contracts at end of period is incremented by how many we buy
        y[contracts] = x[contracts] + u[buy]
        # the Futures price follows a log(AR(1)) process with decreasing variance
        y[price] =
            min(9, max(3, 1.01 * exp(log(x[price]) + σ²[t] * w[alpha])))
        # Production
        y[production] = x[production] + w[beta]

        # return stage objective
        return u[buy] * x[price] - transaction_cost * abs(u[buy])
    end

    terminalobjective!(sp) do x
        return (x[production] - x[contracts]) * x[price]
    end

    constraints!(sp) do x, u, w
        return x[contracts] + u[buy] <= 1.2
    end
end

@time solve(
    m,
    realisation = WaitAndSee,
    riskmeasure = NestedCVaR(beta = 0.5, lambda = 0.5),
)

@time results = simulate(
    m,
    500,              # number of simulations
    contracts = 0,   # Initial states
    price = 4.5, #
    production = 0.0,   #
)

@visualise(
    results,
    t,
    i,
    begin
        results[:objective][i], (title = "Objective", ylabel = "\$")
        results[:contracts][t, i],
        (title = "Number of Contracts", ylabel = "Units")
        results[:price][t, i], (title = "Price", ylabel = "\$/Units")
        results[:production][t, i],
        (title = "Production to Date", ylabel = "Units")
        results[:buy][t, i], (title = "Quantity to buy", ylabel = "Units")
        results[:alpha][t, i], (title = "Alpha Noise")
        results[:beta][t, i], (title = "Beta Noise")
    end
)
