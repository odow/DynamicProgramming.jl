#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=
    Example 7.2 from
        Bertsekas, D. (2005). Dynamic Programming and Optimal Control:
        Volume I (3rd ed.). Bellmont, MA: Athena Scientific.

    A quarterback can choose between running and passing the ball on any given
    play. The number of yards gained by running is integer and is Poisson
    distributed wiht parameter λᵣ. A pass is incomplete with probability p, is
    intercepted with probabilty q, and is completed with probability 1 - p - q.
    When completed, a pass gains an integer number of yards that is Poisson
    distributed with parameter  λₚ. We assume that the probabilty of scoring a
    touchdown on a single play starting i yards from the goal is equal to the
    probability of gaining a number of yards greater  than or equal to i. We
    assume also that the yardage cannot be lost on any play and  that there are
    no penalties. The ball is tured over to the other team on a fourth down or
    when an interception occurs.
=#

using DynamicProgramming, Test

m = SDPModel(
        stages = 4, # four downs
        sense  = :Max
            ) do sp, t
    p, q, λᵣ, λₚ = 0.4, 0.05, 3, 10
    # Poisson probability mass function
    pmf(k, λ) = convert(Float64, (λ^k) * exp(-λ) / factorial(BigInt(k)))
    # Not very efficient. We sample 51 x 51 = 2601 noise terms each stage
    # we could instead make the noise terms linspace(0, 1, N) and associate
    # each random value with a meter.
    YARDS = 0:50
    run_prob  = pmf.(YARDS, λᵣ)  # compute probability
    run_prob ./= sum(run_prob)   # normalize
    pass_prob = pmf.(YARDS, λₚ)  # compute probability
    pass_prob ./= sum(pass_prob) # normalize

    @states(sp, begin
        yards    in 0:100
        has_ball in [0, 1]
    end)

    @controls(sp, begin
        play in [:none, :run, :pass]
    end)

    @noises(sp, begin
        pass_result in DiscreteDistribution([:complete, :interception, :incomplete], [p, q, 1-p-q])
        pass_yard   in DiscreteDistribution(YARDS, pass_prob)
        run_yard    in DiscreteDistribution(YARDS, run_prob)
    end)

    dynamics!(sp) do y, state, control, noise
        yards_play, possession = 0, true
        if control[play] == :run
            yards_play = noise[run_yard]
        elseif control[play] == :pass && noise[pass_result] == :complete
            yards_play = noise[pass_yard]
        elseif control[play] == :pass && noise[pass_result] == :interception
            possession = false
        end
        y[yards] = min(100, state[yards] + yards_play)
        if y[yards] > 99.5 || !possession
            y[has_ball] = 0.0
        else
            y[has_ball] = 1.0
        end
        return 0.0
    end

    terminalobjective!(sp) do state
        state[yards] > 99.5 ? 1.0 : 0.0
    end
end

@time solve(m, realisation=HereAndNow)
win_probability = [
    DynamicProgramming.getbound(m, d, has_ball=1.0, yards=i)
    for i in 0:99, d in 1:4
]
for i in 1:size(win_probability, 1)-1
    for j in 1:size(win_probability, 2)-1
        # win probability increases with more downs
        @test win_probability[i,j] >= win_probability[i,j+1] - 1e-4
        # win probability increases closer to end-zone
        @test win_probability[i,j] <= win_probability[i+1,j] + 1e-4
    end
end
