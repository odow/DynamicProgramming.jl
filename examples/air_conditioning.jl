#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=
    The air conditioning example from Anthony Papavasiliou
    https://perso.uclouvain.be/anthony.papavasiliou/public_html/SDDP.pdf

    Consider the following problem
        Produce air conditioners for 3 months
        200 units/month at 100 $/unit
        Overtime costs 300 $/unit
        Known demand of 100 units for period 1
        Equally likely demand, 100 or 300 units, for periods 2, 3
        Storage cost is 50 $/unit
        All demand must be met

    Optimal bound $62,500
=#
using DynamicProgramming, Base.Test

function airconditioningmodel()
    m = SDPModel(
                 stages = 3,
                  sense = :Max,
                            ) do sp, stage
        DEMAND = [
            [100],
            [100, 300],
            [100, 300]
        ]

        @states(sp, begin
            # number of units
            storage in [0, 100, 200, 300]
        end)

        @controls(sp, begin
            # number of units produced during normal production
            # circumstances
            production in [0, 100, 200]
            # number of units produced during overtime production
            # circumstances
            overtime   in [0, 100, 200, 300]
        end)

        @noises(sp, begin
            demand in DEMAND[stage]
        end)

        dynamics!(sp) do y, x, u, w
            y[storage] = x[storage] + u[production] + u[overtime] - w[demand]
            return -(100 * u[production] + 300 * u[overtime] + 50 * y[storage])
        end

        constraints!(sp) do x, u, w
            x[storage] + u[production] + u[overtime] - w[demand] >= 0
        end
    end
end

srand(1234)
m = airconditioningmodel()
solve(m, realisation=WaitAndSee)
@test isapprox(DynamicProgramming.getbound(m, 1, storage=0), -62_500.0)
