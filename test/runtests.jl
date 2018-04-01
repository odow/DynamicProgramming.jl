#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using DynamicProgramming, Base.Test

const DP = DynamicProgramming

@testset "Simple helper functions" begin
    @test DP.comparison(DP.Minimisation, 1., 2.) == true
    @test DP.comparison(DP.Minimisation, 2., 1.) == false
    @test DP.comparison(DP.Minimisation, 0, 0) == false

    @test DP.comparison(DP.Maximisation, 1., 2.) == false
    @test DP.comparison(DP.Maximisation, 2., 1.) == true
    @test DP.comparison(DP.Maximisation, 0, 0) == false
end

const examples_dir = joinpath(dirname(dirname(@__FILE__)), "examples")

const Examples = Dict(
    "AppliedDynamicProgramming" => [
            "cargo-loading.jl",
            "flyaway-kit.jl"
        ],
    "Bertsekas" => [
        "inventory_control.jl"
    ],
    "FAST" => [
        "hydro-thermal.jl",
        "production-management.jl"
    ],
    "" => [
        "air_conditioning.jl"
    ]
)

@testset "Examples" begin
    for (key, examples) in Examples
        @testset "$(key)" begin
            for example in examples
                @testset "$example" begin
                    println("Running $(key):$(example)")
                    include(joinpath(examples_dir, key, example))
                end
            end
        end
    end
end
