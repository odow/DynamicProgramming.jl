#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

# This example is a modified version of that found at
# https://github.com/JuliaOpt/StochDynamicProgramming.jl/blob/v0.2.0/examples/battery_storage_parallel.jl
# which was orignally authored by Vincent Leclere, Francois Pacaud and Henri Gerard
#
#   To run with N processors, use
#       julia -p N battery_storage.jl
#
using DynamicProgramming

@everywhere begin
    srand(1111)
    const discharge_efficiency = 0.97
    const charge_efficiency    = 0.98
    const cost   = rand(50)
    const demand = rand(50)
end

m = SDPModel(
        stages = 50,
        sense  = :Min
            ) do sp, t

    @states(sp, begin
        storage in 0.:0.01:1
    end)

    @controls(sp, begin
        charge in 0:0.05:1
    end)

    @noises(sp, begin
        # charge leakage
        xi in linspace(0, 0.1, 10)
    end)

    dynamics!(sp) do y, x, u, w
        # update state
        y[storage] = x[storage] + charge_efficiency * u[charge] - (demand[t] / discharge_efficiency + w[xi])

        # return stage objective
        return cost[t] * u[charge] + 2 * abs(min(0, y[storage]))
    end

    constraints!(sp) do x, u, w
        x[storage] + charge_efficiency * u[charge] - (demand[t] / discharge_efficiency + w[xi]) <= 1
    end

end

println("Expected Value")
@time solve(m, realisation=ExpectedValue)

@time results = simulate(m,
    200,          # number of simulations
    storage = 0.5 # initial storage
)

@visualise(results, t, i, begin
    results[:objective][i], (title="Objective",                    ylabel="\$")
    results[:storage][t,i], (title="Storage Level of the Battery", ylabel="")
    results[:charge][t,i],  (title="Charging Decision",            ylabel="")
    results[:xi][t,i],      (title="Random Noise",                 ylabel="")
end)

println("Here and Now")
@time solve(m, realisation=HereAndNow)

@time results = simulate(m,
    200,          # number of simulations
    storage = 0.5 # initial storage
)

@visualise(results, t, i, begin
    results[:objective][i], (title="Objective",                    ylabel="\$")
    results[:storage][t,i], (title="Storage Level of the Battery", ylabel="")
    results[:charge][t,i],  (title="Charging Decision",            ylabel="")
    results[:xi][t,i],      (title="Random Noise",                 ylabel="")
end)

println("Wait and See")
@time solve(m, realisation=WaitAndSee)

@time results = simulate(m,
    200,          # number of simulations
    storage = 0.5 # initial storage
)

@visualise(results, t, i, begin
    results[:objective][i], (title="Objective",                    ylabel="\$")
    results[:storage][t,i], (title="Storage Level of the Battery", ylabel="")
    results[:charge][t,i],  (title="Charging Decision",            ylabel="")
    results[:xi][t,i],      (title="Random Noise",                 ylabel="")
end)
