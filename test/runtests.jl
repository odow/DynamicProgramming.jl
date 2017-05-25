#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using DynamicProgramming
using Base.Test

const DP = DynamicProgramming

@testset "Simple helper functions" begin
    @test DP.comparison(DP.Minimisation, 1., 2.) == true
    @test DP.comparison(DP.Minimisation, 2., 1.) == false
    @test DP.comparison(DP.Minimisation, 0, 0) == false

    @test DP.comparison(DP.Maximisation, 1., 2.) == false
    @test DP.comparison(DP.Maximisation, 2., 1.) == true
    @test DP.comparison(DP.Maximisation, 0, 0) == false
end

for example in ["simple_example.jl", "battery_storage.jl"]
    include(joinpath(dirname(dirname(@__FILE__)), "examples", example))
end
