## Simulate

Once a `SDPModel` has been solved, it is possible to simulate the performance of the policy using the function `simulate(m::SDPModel, N::Int; kwargs...)`.
`m` is the solved `SDPModel` to be simulated. `N` is the number of realisations to perform. Initial values for the state variables are given via the keyword arguments.

For example:

```julia
results = simulate(m,
    500,
    contracts  = 0,
    price      = 4.5,
    production = 0.
)
```
