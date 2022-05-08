#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=

    This example comes from Chapter I, Sections 22-26.
        Bellman, R. & Dreyfus, S. (1962). Applied Dynamic
        Programming. Princeton, NJ: Princeton University Press.

    Let us suppose that we are engaged in loading a vessel with
    a cargo comprised of different types of items. Since these
    different items have various weights and values, the problem
    arises as to how to load a ship of limited capacity with the
    most valuable cargo. The reader, by means of a simple
    transformation of the situations, will be able to conceive
    many questions of similar natures.
=#

using DynamicProgramming, Test

m = SDPModel(sense=:Max, stages=8) do sp, t
    weights = [20, 18, 14, 12, 10, 16, 22, 24]
    values  = [72, 60, 40, 27, 20, 50, 85, 96]
    @states(sp, begin
        remaining_space in 0.0:1.0:100.0
    end)
    @controls(sp, begin
        item in 0.0:1.0:5.0
    end)
    dynamics!(sp) do y, x, u, w
        y[remaining_space] = x[remaining_space] - weights[t] * u[item]
        return values[t] * u[item]
    end
    constraints!(sp) do x, u, w
        x[remaining_space] >= weights[t] * u[item]
    end
end

solve(m, print_level=0)

solutions = Dict(
    100 => (384, [(8,4)]),
    95  => (373, [(7,1), (8,3)]),
    91  => (351, [(7,3), (8,1)]),
    89  => (340, [(7,4)]),
    87  => (328, [(3,1), (8,3)]),
    85  => (317, [(3,1), (7,1), (8,2)]),
    83  => (308, [(5,1), (8,3)]),
    81  => (297, [(5,1), (7,1), (8,2)])
)
for (z, sol) in solutions
    s = simulate(m, 1, remaining_space=z)
    @test s[:objective][1] == sol[1]
    for (i, v) in sol[2]
        @test s[:item][i] == v
    end
end
