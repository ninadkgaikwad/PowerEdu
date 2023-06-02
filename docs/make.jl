using PowerSystemsAnalysis
using Documenter

DocMeta.setdocmeta!(PowerSystemsAnalysis, :DocTestSetup, :(using PowerSystemsAnalysis); recursive=true)

makedocs(;
    modules=[PowerSystemsAnalysis],
    authors="Ninad Kiran Gaikwad",
    repo="https://github.com/ninadkgaikwad/PowerSystemsAnalysis.jl/blob/{commit}{path}#{line}",
    sitename="PowerSystemsAnalysis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ninadkgaikwad.github.io/PowerSystemsAnalysis.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ninadkgaikwad/PowerSystemsAnalysis.jl",
    devbranch="main",
)
