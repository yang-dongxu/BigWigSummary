using BigWigSummary
using Documenter

DocMeta.setdocmeta!(BigWigSummary, :DocTestSetup, :(using BigWigSummary); recursive=true)

makedocs(;
    modules=[BigWigSummary],
    authors="Dongxu Yang",
    repo="https://github.com/yang-dongxu/BigWigSummary.jl/blob/{commit}{path}#{line}",
    sitename="BigWigSummary.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://yang-dongxu.github.io/BigWigSummary.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/yang-dongxu/BigWigSummary.jl",
    devbranch="main",
)
