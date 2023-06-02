using PowerEdu
using Documenter

DocMeta.setdocmeta!(PowerEdu, :DocTestSetup, :(using PowerEdu); recursive=true)

makedocs(;
    modules=[PowerEdu],
    authors="Ninad Kiran Gaikwad",
    repo="https://github.com/ninadkgaikwad/PowerEdu.jl/blob/{commit}{path}#{line}",
    sitename="PowerEdu.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ninadkgaikwad.github.io/PowerEdu.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ninadkgaikwad/PowerEdu.jl",
    devbranch="main",
)
