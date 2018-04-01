#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#==
    An implementation of the Production Management example from FAST
    https://github.com/leopoldcambier/FAST/blob/daea3d80a5ebb2c52f78670e34db56d53ca2e778/examples/production management multiple stages/
==#

using DynamicProgramming, Base.Test

m = SDPModel(
                sense  = :Min,
                stages = 3
                    ) do sp, t

    DEMAND = t==1 ? [0] : [2, 10]
    N = 10
    @states(sp, begin
        x1 in 0:2:N
        x2 in 0:2:N
    end)
    @controls(sp, begin
        s1    in 0:2:N
        s2    in 0:2:N
        newx1 in 0:2:N
        newx2 in 0:2:N
    end)
    @noises(sp, begin
        demand in DEMAND
    end)
    constraints!(sp) do state, control, noise
        all([
            control[s1] <= state[x1],
            control[s2] <= state[x2],
            control[s1] + control[s2] <= noise[demand]
        ])
    end
    dynamics!(sp) do state_out, state_in, control, noise
        state_out[x1] = control[newx1]
        state_out[x2] = control[newx2]
        return (0.2 * state_out[x1] + 0.7 * state_out[x2]) -
            (2.33 * control[s1] + 2.54 * control[s2])
    end
end

solve(m, realisation=WaitAndSee, print_level=0)
@test isapprox(DynamicProgramming.getbound(m, 1, x1=0,x2=0), -23.96)
