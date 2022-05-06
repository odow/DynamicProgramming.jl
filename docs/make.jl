using Documenter, DynamicProgramming

makedocs(
    # modules = [DynamicProgramming],
    clean = false,
    format = Documenter.HTML(
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://odow.github.io/DynamicProgramming.jl",
        assets=String[],
    ),
    sitename = "DynamicProgramming.jl",
    authors = "Oscar Dowson",
    pages = [
        "Home" => "index.md",
        "Model Formulation" => "model_formulation.md",
        "Solve" => "solve.md",
        "Simulate" => "simulate.md",
        "Visualise" => "visualisation.md"
    ]
)
