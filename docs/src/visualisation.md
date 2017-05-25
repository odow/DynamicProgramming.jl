## Visualise

It is possible to create an interactive visualisation of the simulated policy with the `@visualise` macro. The following keywords should be wrapped with parentheses.
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
