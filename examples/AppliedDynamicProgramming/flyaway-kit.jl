#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=

    This example comes from Chapter II, Section 8.
        Bellman, R. & Dreyfus, S. (1962). Applied Dynamic
        Programming. Princeton, NJ: Princeton University Press.

    Let us assume that we are about to dispatch a transport plane
    to an overseas base with a cargo which is to consist of
    replacement parts for airplanes. Suppose that there are N types
    of replacement parts, and that each has associated with it a cost
    that is incurred if that part is needed at the base, but not
    available. Let us futher assume that the demand for eah part can
    be described by a known Poisson distribution. As constraints by
    the weight capacity W and the space capacity S of the cargo
    vehicle, how many items of each type shall be loaded so as to
    minimize the expected cost due to shortages at the base?
        This problem is typical of many that are met in the study of
    inventory and equipment replacement processses.

    N.B.: somewhat annoyingly, the book fails to give the values
    for the costs and lambda so we just make some up.
=#

using DynamicProgramming, Test

m = SDPModel(sense = :Min, stages = 10) do sp, t
    weights = [3, 5, 2, 5, 6, 2, 4, 7, 5, 3]
    sizes = [5, 4, 7, 4, 6, 2, 5, 7, 3, 4]
    costs = [22, 22, 12, 6, 12, 2, 12, 12, 22, 4]
    lambda = [1, 1, 4, 4, 1, 4, 1, 1, 1, 3]
    @states(sp, begin
        weight in 0.0:1.0:30.0
        size in 0.0:1.0:30.0
    end)
    @controls(sp, begin
        quantity in 0.0:1.0:3.0
    end)
    P(k, i) = exp(-lambda[i]) * lambda[i]^k / factorial(k)

    dynamics!(sp) do y, x, u, w
        y[weight] = x[weight] - weights[t] * u[quantity]
        y[size] = x[size] - sizes[t] * u[quantity]
        return costs[t] * sum(
            (k - u[quantity]) * P(k, t) for k in Int(u[quantity] + 1):10
        )
    end
    constraints!(sp) do x, u, w
        return all([
            x[weight] >= weights[t] * u[quantity],
            x[size] >= sizes[t] * u[quantity],
        ])
    end
end

solve(m, print_level = 0)
s = simulate(m, 1, weight = 30.0, size = 30.0)
@test isapprox(s[:objective][1], 122.37, atol = 1e-2)
