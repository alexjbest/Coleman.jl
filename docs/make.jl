using Documenter, Coleman

makedocs(;
    modules=[Coleman],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/alexjbest/Coleman.jl/blob/{commit}{path}#L{line}",
    sitename="Coleman.jl",
    authors="Alex J. Best",
    assets=[],
)

deploydocs(;
    repo="github.com/alexjbest/Coleman.jl",
)
