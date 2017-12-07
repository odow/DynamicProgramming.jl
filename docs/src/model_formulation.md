## Initialise Model

A SDP model is stored in the `SDPModel` type. A `SDPModel` can be intialised using the following `SDPModel() do ... end` block structure:

```julia
m = SDPModel(
        stages = 1.   # Int
        sense  = :Min # Symbol (:Min or :Max)
    ) do sp, t

        ... problem definition ...

end
```

`sp` can be given any name but refers to the *stage problem*. `t` can also be given any name, but refers to an index that runs from `1` to `T` where `T` is the number of stages.

Inside the `SDPModel` definition, we define our subproblems. We first need to add some state variables, some control (or action) variables, and some noise (or stochastic) variables.


## Initialise Variables

States can be added with the following macro:
```julia
@addstates!(sp, begin
    x = linspace(0, 1, 10)
end)
```

This creates a state variable `x` that is discretised into the set `linspace(0, 1, 10)`. Note that currently, all state dimensions get converted into `Float64` representations. The discretisation should be any type that can be converted to a `Vector{Float64}` type.

Controls can be added with the `@addcontrols!` macro that has similar syntax. However there is less restriction on the type. The discretisation should just be an iterable subtype of `AbstractVector`.

Noise (or stochastic variables) can be added with the `@addnoises!` macro:

```julia
@addnoises!(sp, begin
    u = DiscreteDistribution([1,2,3], [0.5, 0.25, 0.25])
    v = 1:10
end)
```

In contrast to the other two macros, there is a slight subtlety. The discretisations can either be subtypes of `AbstractVector` (in which case their realisations are assumed to be uniformly sampled), or a `DiscreteDistribution`.

The `WeighedDist` constructor is `DiscreteDistribution(values::AbstractVector, probability::AbstractVector)`. This realisations `values` are sampled with probability `probability`.

If more than one noise is defined, then the multiple noises are assumed to be independent.

## Dynamics

You must provide a function that takes four inputs.

```julia
function foo(states_out, states, controls, noises)
end
```

However, we prefer the anonymous function syntax:
```julia
dynamics!(sp, (y, x, u, w) -> (
        ... definitions ...
    )
)
```

This anonymous function must take the current state `x`,  a control `u` and a noise `w` and update the new state `y`.

You can refer to model variables using the `[]` indexing operator. For example, if we defined a state variable `quantity`, we could refer it as `x[quantity]`.

By thinking about variable scopes it is possible to encapsulate all the necessary data into this syntax.

## Stage Objective

The stage objective function takes as input vectors for the initial state, control and noise. It should return a single `Float64` value corresponding to the cost (or profit) accrued in the current stage.

```julia
stageobjective!(sp, (x, u, w) -> (
        ... definitions ...
    )
)
```

## Terminal Objective

The terminal objective function takes as input a vector of the final state at the end of the finite time horizon. It returns a single `Float64` value corresponding to the cost (or profit) of ending in that state.

```julia
terminalobjective!(sp, (x) -> (
        ... definitions ...
    )
)
```

## Constraints

The constraints function takes as input vectors for the initial state, control and noise. It should return a single `Bool` value indicating if the state, control, noise combination is feasible. Typically this can be implemented by a chained series of boolean comparisons.

```julia
constraints!(sp, (x, u, w) -> (
        ... definitions ...
    )
)
```
