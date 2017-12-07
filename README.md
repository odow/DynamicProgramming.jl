# DynamicProgramming.jl

A Julia package for discrete stochastic dynamic programming.

## Installation

This package is not yet registered. **Until it is, things may change. It is perfectly
usable but should not be considered stable**.

You can install it by typing

```julia
julia> Pkg.clone("https://github.com/odow/DynamicProgramming.jl.git")
```

## Initialise Model

A SDP model is stored in the `SDPModel` type. A `SDPModel` can be intialised
using the following `SDPModel() do ... end` block structure:

```julia
m = SDPModel(
        stages = 1.   # Int
        sense  = :Min # Symbol (:Min or :Max)
    ) do sp, t

        ... problem definition ...

end
```

`sp` can be given any name but refers to the *stage problem*. `t` can also be
given any name, but refers to an index that runs from `1` to `T` where `T` is
the number of stages.

Inside the `SDPModel` definition, we define our subproblems. We first need to
add some state variables, some control (or action) variables, and some noise (or
     stochastic) variables.


## Initialise Variables

States can be added with the following macro:
```julia
@states!(sp, begin
    x in linspace(0, 1, 10)
end)
```

This creates a state variable `x` that is discretised into the set
`linspace(0, 1, 10)`. Note that currently, all state dimensions get converted
into `Float64` representations. The discretisation should be any type that can
be converted to a `Vector{Float64}` type.

Controls can be added with the `@controls!` macro that has similar syntax.
However there is less restriction on the type. The discretisation should just be
an iterable subtype of `AbstractVector`.

Noise (or stochastic variables) can be added with the `@noises!` macro:

```julia
@noises!(sp, begin
    u in DiscreteDistribution([1,2,3], [0.5, 0.25, 0.25])
    v in 1:10
end)
```

In contrast to the other two macros, there is a slight subtlety. The
discretisations can either be subtypes of `AbstractVector` (in which case their
    realisations are assumed to be uniformly sampled), or a `DiscreteDistribution`.

The `DiscreteDistribution` constructor is
`DiscreteDistribution(observations::AbstractVector, probability::AbstractVector)`.
This realisations `observations` are sampled with probability `probability`.

If more than one noise is defined, then the multiple noises are assumed to be
independent.

## Dynamics

You must provide a function that takes four inputs.

```julia
function foo(states_out, states, controls, noises)
end
```

However, we prefer the anonymous function syntax:
```julia
dynamics!(sp) do y, x, u, w
        ... definitions ...
end
```

This anonymous function must take the current state `x`,  a control `u` and a
noise `w` and update the new state `y`.

You can refer to model variables using the `[]` indexing operator. For example,
if we defined a state variable `quantity`, we could refer it as `x[quantity]`.

By thinking about variable scopes it is possible to encapsulate all the
necessary data into this syntax.

**Important:** the function should return a single `Float64` value corresponding
 to the cost (or profit) accrued in the current stage.

## Terminal Objective

The terminal objective function takes as input a vector of the final state at
the end of the finite time horizon. It returns a single `Float64` value
corresponding to the cost (or profit) of ending in that state.

```julia
terminalobjective!(sp) do x
        ... definitions ...
end
```

## Constraints

The constraints function takes as input vectors for the initial state, control
and noise. It should return a single `Bool` value indicating if the state,
control, noise combination is feasible. Typically this can be implemented by a
chained series of boolean comparisons.

```julia
constraints!(sp) do x, u, w
        ... definitions ...
end
```

## Solve

The `solve` function takes as input the initalised `SDPModel` object, as well as
 two keyword arguments.

The `realisation` must be one of `WaitAndSee`, `HereAndNow` or `ExpectedValue`.

The `WaitAndSee` model observes the noise before choosing the optimal control.
The `HereAndNow` model chooses the best control before observing the noise.
The `ExpectedValue` model substitutes the noise for the expected value of each
 independent noise.

The `riskmeasure` is a nested `λE[x] + (1-λ)CVaRᵦ[x]`.

```julia
solve(m::SDPModel,
    realisation=WaitAndSee,
    riskmeasure=NestedCVaR(beta=0.5, lambda=0.5)
)
```

## Simulate

Once a `SDPModel` has been solved, it is possible to simulate the performance of
 the policy using the function `simulate(m::SDPModel, N::Int; kwargs...)`.
`m` is the solved `SDPModel` to be simulated. `N` is the number of realisations
to perform. Initial values for the state variables are given via the keyword
arguments.

For example:

```julia
results = simulate(m,
    500,
    contracts  = 0,
    price      = 4.5,
    production = 0.
)
```

## Visualise

It is possible to create an interactive visualisation of the simulated policy
with the `@visualise` macro. The following keywords should be wrapped with
parentheses.
 - `"cumulative"  = false` Plot the cumulation of the variable over stages
 - `"title"       = ""` Plot title
 - `"xlabel"      = "Stages"` Label for x axis
 - `"ylabel"      = ""` Label for y axis
 - `"interpolate" = "linear"` D3.js interpolation method to use. See the [D3 wiki](https://github.com/d3/d3/wiki/SVG-Shapes#line_interpolate) for more.

The following example gives an example of possible syntax:
```julia
@visualise(results, stage, replication, begin
    results[:Current][stage][replication],  (title="Accumulated Profit", ylabel="Accumulated Profit (\$)", cumulative=true)
    results[:x][stage][replication],    (title="Value of a State", ylabel="Level")
    results[:u][stage][replication],    (title="Value of a Control")
    results[:w][stage][replication],    (title="Value of a Noise", interpolate="step")
end)
```
