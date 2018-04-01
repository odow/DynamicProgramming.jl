#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=
	inspired by @flowald's problem
https://discourse.julialang.org/t/parametrize-jump-problem-for-repeat-evaluation

	This is very similar to AppliedDynamicProgramming/cargo-loading.jl.
=#

using DynamicProgramming

TOT_COURSES = 10
MAX_COURSES = 5
BUDGET      = 3
price = collect(linspace(0.1,1,TOT_COURSES))
prefs = collect(linspace(0.1,1,TOT_COURSES))

model = SDPModel(stages=TOT_COURSES, sense=:Max) do sp, t
	@states(sp, begin
		budget in 0:0.1:BUDGET
		classes in 0:1:MAX_COURSES
	end)
	@controls(sp, begin
		accept_course in [0, 1]
	end)
	dynamics!(sp) do state_out, state_in, control, noise
		state_out[budget] = state_in[budget] - control[accept_course] * price[t]
		state_out[classes] = state_in[classes] + control[accept_course]
		return prefs[t] * control[accept_course]
	end
	constraints!(sp) do state_in, control, noise
		all([
			state_in[classes] + control[accept_course] <= MAX_COURSES,
			state_in[budget] - control[accept_course] * price[t] >= 0
		])
	end
end

solve(model)
println("Objective is: ", DynamicProgramming.getbound(model, 1, classes=0,budget=3))
println("Solution is:")
sim = simulate(model, 1, classes=0, budget=3)
println(sim[:accept_course][:])
println("total cost = ", dot(price, sim[:accept_course][:]))
