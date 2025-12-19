using Documenter
using ConstantFunctionMarketMakers

makedocs(
    sitename = "ConstantFunctionMarketMakers.jl",
    modules = [ConstantFunctionMarketMakers],
    repo = Remotes.GitHub("sschuldenzucker", "ConstantFunctionMarketMakers.jl"),
    pages = ["Home" => "index.md", "API Reference" => "api.md", "Internal API" => "api_internal.md"],
    format = Documenter.HTML(
        # I don't like these and they make it harder to browse locally.
        prettyurls = false,
    ),
    # checkdocs=:public,  # Check that all exported or public identifiers with a docstring actually occur in the docs
    checkdocs_ignored_modules=[
        ConstantFunctionMarketMakers.Util,
        ConstantFunctionMarketMakers.ECLPMath,
        ConstantFunctionMarketMakers.StableSwapMath,
        ConstantFunctionMarketMakers.Variables,
    ]
)
