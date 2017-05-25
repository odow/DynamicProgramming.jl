using Documenter, DynamicProgramming

makedocs(
    # modules = [DynamicProgramming],
    format = :html,
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
