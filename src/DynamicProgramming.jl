#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

__precompile__()

module DynamicProgramming

using Interpolations, JSON

export
    SDPModel, solve, simulate,
    # macros
    @states!, @controls!, @noises!, @visualise,
    # model functions
    dynamics!, terminalobjective!, constraints!,
    WeightedDist,
    # risk measures
    Expectation, NestedCVaR,
    # uncertainty realisations
    WaitAndSee, HereAndNow, ExpectedValue

include("MIT_licencedcode.jl")
include("type_definitions.jl")
include("macro_definitions.jl")
include("function_definitions.jl")
include("visualise.jl")

end
