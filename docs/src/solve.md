## Solve

The `solve` function takes as input the initalised `SDPModel` object, as well as two keyword arguments.

The `realisation` must be one of `WaitAndSee`, `HereAndNow` or `ExpectedValue`.

The `WaitAndSee` model observes the noise before choosing the optimal control.
The `HereAndNow` model chooses the best control before observing the noise.
The `ExpectedValue` model substitutes the noise for the expected value of each independent noise.

The `riskmeasure` is a nested `λE[x] + (1-λ)CVaRᵦ[x]`.

```julia
solve(m::SDPModel,
    realisation=WaitAndSee,
    riskmeasure=NestedCVaR(beta=0.5, lambda=0.5)
)
```
