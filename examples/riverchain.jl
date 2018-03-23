#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using DynamicProgramming

m = SDPModel(
    stages   = 12,
    sense    = :Max
            ) do sp, t

    function generationquantity(v)
        if 0 <= v <= 50
            return 1.1v
        elseif 50 < v <= 60
            return v + 5.0
        elseif 60 < v <= 70
            0.5 * (v - 60) + 65.0
        else
            0.0
        end
    end
    @states!(sp, begin
          upstream in 0:5:100
        downstream in 0:5:100
             price in vcat(collect(50:0.5:61), 61.261, collect(62:0.5:90))
    end)

    @controls!(sp, begin
          uprelease in 0:5:70
        downrelease in 0:5:70
    end)

    @noises!(sp, begin
        noise in vcat(-[0.5,1.5,2.5,3.5,4.5], [0.5,1.5,2.5,3.5,4.5])
    end)

    dynamics!(sp) do y, x, u, w
        y[upstream]    = x[upstream] - u[uprelease]
        y[downstream]  = x[downstream] - u[downrelease] + u[uprelease]
        spill = 0.0
        upstream_spill = max(0.0, y[upstream] - 200.0)
        if upstream_spill > 0.0
            spill += upstream_spill
            y[upstream] = 200.0
        end
        downstream_spill = max(0.0, y[downstream] - 200.0)
        if downstream_spill > 0.0
            spill += downstream_spill
            y[downstream] = 200.0
        end
        b_t = [61.261, 56.716, 59.159, 66.080, 72.131, 76.708, 76.665,
            76.071, 76.832, 69.970, 69.132, 67.176]
        if t > 1
            y[price] = 0.5 * x[price] + 0.5 * b_t[t] + w[noise]
        else
            y[price] = x[price]
        end

        # return stage objective
        total_generation = generationquantity(u[uprelease] - upstream_spill) + generationquantity(u[downrelease] - downstream_spill)

        return total_generation * y[price] - 1000.0 * spill
    end

    terminalobjective!(sp) do x
        0.0
    end

    constraints!(sp) do x, u, w
        (x[upstream] - u[uprelease] >= 0) &&
        (x[downstream] - u[downrelease] >= 0)
    end

end

@time solve(m,
    realisation=WaitAndSee
)

@time results = simulate(m,
    1000,
    upstream   = 100.0,
    downstream = 100.0,
    price      = 61.261
)

@visualise(results, t, i, begin
    results[:upstream][t,i],    (title="Upstream", ylabel="")
    results[:downstream][t,i],  (title="Downstream",            ylabel="")
    results[:price][t,i],       (title="Price",                 ylabel="")
    results[:uprelease][t,i],   (title="Uprelease", ylabel="")
    results[:downrelease][t,i], (title="Downrelease",            ylabel="")
end)

@show m.stages[1].interpolatedsurface[100.0, 100.0, 61.261]
