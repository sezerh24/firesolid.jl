using firesolid
using Documenter

DocMeta.setdocmeta!(firesolid, :DocTestSetup, :(using firesolid); recursive=true)

makedocs(;
    modules=[firesolid],
    authors="Chris de Graaf, Invenia Technical Computing Corporation",
    sitename="firesolid.jl",
    format=Documenter.HTML(;
        canonical="https://sezerh24.github.io/firesolid.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/sezerh24/firesolid.jl",
    devbranch="main",
)
