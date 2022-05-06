#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#==
    An implementation of the Hydro-thermal example from FAST
    https://github.com/leopoldcambier/FAST/tree/daea3d80a5ebb2c52f78670e34db56d53ca2e778/examples/hydro%20thermal
==#

using DynamicProgramming

m = SDPModel(
                sense  = :Min,
                stages = 2
                    ) do sp, t
    RAINFALL = t == 1 ? [6] : [2, 10]
    @states(sp, begin
        x in 0:1:8
    end)
    @controls(sp, begin
        p in 0:1:6
        y in 0:1:18
    end)
    @noises(sp, begin
        rainfall in RAINFALL
    end)
    dynamics!(sp) do state_out, state_in, control, noise
        state_out[x] <= state_in[x] + noise[rainfall] - control[y]
        return 5.0 * control[p]
    end
    constraints!(sp) do state_in, control, noise
        all([
            control[p] + control[y] >= 6,
            state_in[x] + noise[rainfall] - control[y] >= 0
        ])
    end
end

solve(m, realisation=WaitAndSee, print_level=0)
println(DynamicProgramming.getbound(m,1,x=0))
