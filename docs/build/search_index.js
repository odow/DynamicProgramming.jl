var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#DynamicProgramming.jl-1",
    "page": "Home",
    "title": "DynamicProgramming.jl",
    "category": "section",
    "text": "A Julia package for discrete stochastic dynamic programming."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "This package is not yet registered. You can install it by typingjulia> Pkg.clone(\"https://github.com/odow/DynamicProgramming.jl.git\")"
},

{
    "location": "model_formulation.html#",
    "page": "Model Formulation",
    "title": "Model Formulation",
    "category": "page",
    "text": ""
},

{
    "location": "model_formulation.html#Initialise-Model-1",
    "page": "Model Formulation",
    "title": "Initialise Model",
    "category": "section",
    "text": "A SDP model is stored in the SDPModel type. A SDPModel can be intialised using the following SDPModel() do ... end block structure:m = SDPModel(\n        stages = 1.   # Int\n        sense  = :Min # Symbol (:Min or :Max)\n    ) do sp, t\n\n        ... problem definition ...\n\nendsp can be given any name but refers to the stage problem. t can also be given any name, but refers to an index that runs from 1 to T where T is the number of stages.Inside the SDPModel definition, we define our subproblems. We first need to add some state variables, some control (or action) variables, and some noise (or stochastic) variables."
},

{
    "location": "model_formulation.html#Initialise-Variables-1",
    "page": "Model Formulation",
    "title": "Initialise Variables",
    "category": "section",
    "text": "States can be added with the following macro:@addstates!(sp, begin\n    x = linspace(0, 1, 10)\nend)This creates a state variable x that is discretised into the set linspace(0, 1, 10). Note that currently, all state dimensions get converted into Float64 representations. The discretisation should be any type that can be converted to a Vector{Float64} type.Controls can be added with the @addcontrols! macro that has similar syntax. However there is less restriction on the type. The discretisation should just be an iterable subtype of AbstractVector.Noise (or stochastic variables) can be added with the @addnoises! macro:@addnoises!(sp, begin\n    u = WeightedDist([1,2,3], [0.5, 0.25, 0.25])\n    v = 1:10\nend)In contrast to the other two macros, there is a slight subtlety. The discretisations can either be subtypes of AbstractVector (in which case their realisations are assumed to be uniformly sampled), or a WeightedDist.The WeighedDist constructor is WeightedDist(values::AbstractVector, probability::AbstractVector). This realisations values are sampled with probability probability.If more than one noise is defined, then the multiple noises are assumed to be independent."
},

{
    "location": "model_formulation.html#Dynamics-1",
    "page": "Model Formulation",
    "title": "Dynamics",
    "category": "section",
    "text": "You must provide a function that takes four inputs.function foo(states_out, states, controls, noises)\nendHowever, we prefer the anonymous function syntax:dynamics!(sp, (y, x, u, w) -> (\n        ... definitions ...\n    )\n)This anonymous function must take the current state x,  a control u and a noise w and update the new state y.You can refer to model variables using the [] indexing operator. For example, if we defined a state variable quantity, we could refer it as x[quantity].By thinking about variable scopes it is possible to encapsulate all the necessary data into this syntax."
},

{
    "location": "model_formulation.html#Stage-Objective-1",
    "page": "Model Formulation",
    "title": "Stage Objective",
    "category": "section",
    "text": "The stage objective function takes as input vectors for the initial state, control and noise. It should return a single Float64 value corresponding to the cost (or profit) accrued in the current stage.stageobjective!(sp, (x, u, w) -> (\n        ... definitions ...\n    )\n)"
},

{
    "location": "model_formulation.html#Terminal-Objective-1",
    "page": "Model Formulation",
    "title": "Terminal Objective",
    "category": "section",
    "text": "The terminal objective function takes as input a vector of the final state at the end of the finite time horizon. It returns a single Float64 value corresponding to the cost (or profit) of ending in that state.terminalobjective!(sp, (x) -> (\n        ... definitions ...\n    )\n)"
},

{
    "location": "model_formulation.html#Constraints-1",
    "page": "Model Formulation",
    "title": "Constraints",
    "category": "section",
    "text": "The constraints function takes as input vectors for the initial state, control and noise. It should return a single Bool value indicating if the state, control, noise combination is feasible. Typically this can be implemented by a chained series of boolean comparisons.constraints!(sp, (x, u, w) -> (\n        ... definitions ...\n    )\n)"
},

{
    "location": "solve.html#",
    "page": "Solve",
    "title": "Solve",
    "category": "page",
    "text": ""
},

{
    "location": "solve.html#Solve-1",
    "page": "Solve",
    "title": "Solve",
    "category": "section",
    "text": "The solve function takes as input the initalised SDPModel object, as well as two keyword arguments.The realisation must be either WaitAndSee or HereAndNow. A wait and see model observes the noise before choosing the optimal control whereas a here and now model chooses the best control before observing the noise.The riskmeasure is a nested λE[x] + (1-λ)CVaRᵦ[x]solve(m::SDPModel,\n    realisation=WaitAndSee,\n    riskmeasure=NestedCVaR(beta=0.5, lambda=0.5)\n)"
},

{
    "location": "simulate.html#",
    "page": "Simulate",
    "title": "Simulate",
    "category": "page",
    "text": ""
},

{
    "location": "simulate.html#Simulate-1",
    "page": "Simulate",
    "title": "Simulate",
    "category": "section",
    "text": "Once a SDPModel has been solved, it is possible to simulate the performance of the policy using the function simulate(m::SDPModel, N::Int; kwargs...). m is the solved SDPModel to be simulated. N is the number of realisations to perform. Initial values for the state variables are given via the keyword arguments.For example:results = simulate(m,\n    500,\n    contracts  = 0,\n    price      = 4.5,\n    production = 0.\n)"
},

{
    "location": "visualisation.html#",
    "page": "Visualise",
    "title": "Visualise",
    "category": "page",
    "text": ""
},

{
    "location": "visualisation.html#Visualise-1",
    "page": "Visualise",
    "title": "Visualise",
    "category": "section",
    "text": "It is possible to create an interactive visualisation of the simulated policy with the @visualise macro. The following keywords should be wrapped with parentheses.\"cumulative\"  = false Plot the cumulation of the variable over stages\n\"title\"       = \"\" Plot title\n\"xlabel\"      = \"Stages\" Label for x axis\n\"ylabel\"      = \"\" Label for y axis\n\"interpolate\" = \"linear\" D3.js interpolation method to use. See the D3 wiki for more.The following example gives an example of possible syntax:@visualise(results, stage, replication, begin\nresults[:Current][stage][replication],  (title=\"Accumulated Profit\", ylabel=\"Accumulated Profit (\\$)\", cumulative=true)\nresults[:x][stage][replication],    (title=\"Value of a State\", ylabel=\"Level\")\nresults[:u][stage][replication],    (title=\"Value of a Control\")\nresults[:w][stage][replication],    (title=\"Value of a Noise\", interpolate=\"step\")\nend)"
},

]}
