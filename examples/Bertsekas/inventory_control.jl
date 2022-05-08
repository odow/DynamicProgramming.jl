#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=
    Example 1.3.2 from
        Bertsekas, D. (2005). Dynamic Programming and Optimal Control:
        Volume I (3rd ed.). Bellmont, MA: Athena Scientific.
=#

using DynamicProgramming, Test, Statistics, Random

m = SDPModel(stages = 3, sense = :Min) do sp, t
    @states(sp, begin
        xₖ in 0:1:2
    end)

    @controls(sp, begin
        uₖ = 0:1:2
    end)

    @noises(
        sp,
        begin
            wₖ in DiscreteDistribution([0.0, 1.0, 2.0], [0.1, 0.7, 0.2])
        end
    )

    dynamics!(sp) do y, x, u, w
        y[xₖ] = max(0, x[xₖ] + u[uₖ] - w[wₖ])
        return u[uₖ] + (x[xₖ] + u[uₖ] - w[wₖ])^2
    end

    terminalobjective!(sp) do x
        return 0.0
    end

    constraints!(sp) do x, u, w
        return x[xₖ] + u[uₖ] <= 2
    end
end

# ==============================
#   Here-and-Now Solution
# ==============================
solve(m, realisation = HereAndNow, print_level = 0)

here_and_now_solution = [
    3.7 2.5 1.3
    2.7 1.5 0.3
    2.818 1.68 1.1
]
for t in 1:3
    for xk in 0:1:2
        @test isapprox(
            m.stages[t].interpolatedsurface[xk],
            here_and_now_solution[xk+1, t],
            atol = 1e-2,
        )
    end
end

Random.seed!(123)
sims = simulate(m, 1_000_000, xₖ = 0.0)
A = sims[:uₖ] + (sims[:xₖ] + sims[:uₖ] - sims[:wₖ]) .^ 2
@test sum(A, dims = 1)[:] == sims[:objective]
@test isapprox(mean(sims[:objective]), 3.29, atol = 1e-2)

# ==============================
#   Wait-and-See Solution
# ==============================
solve(m, realisation = WaitAndSee, print_level = 0)

wait_and_see_solution = [
    3.300 2.20 1.10
    2.412 1.32 0.30
    2.628 1.64 1.10
]
for t in 1:3
    for xk in 0:1:2
        @test isapprox(
            m.stages[t].interpolatedsurface[xk],
            wait_and_see_solution[xk+1, t],
            atol = 1e-2,
        )
    end
end

Random.seed!(123)
ws_sims = simulate(m, 1_000_000, xₖ = 0.0)
ws_A = ws_sims[:uₖ] + (ws_sims[:xₖ] + ws_sims[:uₖ] - ws_sims[:wₖ]) .^ 2
@test sum(ws_A, dims = 1)[:] == ws_sims[:objective]
@test isapprox(mean(ws_sims[:objective]), 3.29, atol = 1e-2)

# ==============================
#   Expected-Value Solution
# ==============================
solve(m, realisation = ExpectedValue, print_level = 0)

expected_value_solution = [
    3.03 2.02 1.01
    2.03 1.02 0.01
    1.93 0.92 0.81
]
for t in 1:3
    for xk in 0:1:2
        # @show t, xk, m.stages[t].interpolatedsurface[xk]
        @test isapprox(
            m.stages[t].interpolatedsurface[xk],
            expected_value_solution[xk+1, t],
            atol = 1e-2,
        )
    end
end

Random.seed!(123)
ev_sims = simulate(m, 1_000_000, xₖ = 0.0)
ev_A = ev_sims[:uₖ] + (ev_sims[:xₖ] + ev_sims[:uₖ] - ev_sims[:wₖ]) .^ 2
@test sum(ev_A, dims = 1)[:] == ev_sims[:objective]
@test isapprox(mean(ev_sims[:objective]), 3.29, atol = 1e-2)

# ==============================
#   Risk-Averse Here-and-Now
# ==============================
solve(
    m,
    realisation = HereAndNow,
    riskmeasure = NestedCVaR(beta = 1.0, lambda = 0.0),
    print_level = 0,
)

for t in 1:3
    for xk in 0:1:2
        @test isapprox(
            m.stages[t].interpolatedsurface[xk],
            here_and_now_solution[xk+1, t],
            atol = 1e-2,
        )
    end
end

solve(
    m,
    realisation = HereAndNow,
    riskmeasure = NestedCVaR(beta = 0.5, lambda = 0.5),
    print_level = 0,
)

risk_averse_solution = [
    4.05 2.75 1.45
    3.05 1.75 0.45
    3.243 2.035 1.35
]

for t in 1:3
    for xk in 0:1:2
        @test isapprox(
            m.stages[t].interpolatedsurface[xk],
            risk_averse_solution[xk+1, t],
            atol = 1e-2,
        )
    end
end

# ==============================
#   Errors
# ==============================
@test_throws Exception simulate(m, 1, x = 0.0)
