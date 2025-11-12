using Documenter, CyclotomicNumbers

DocMeta.setdocmeta!(CyclotomicNumbers, :DocTestSetup, :(using CyclotomicNumbers); recursive=true)

makedocs(;
    modules=[CyclotomicNumbers],
    authors="Jean Michel <jean.michel@imj-prg.fr>",
    sitename="CyclotomicNumbers.jl",
    format=Documenter.HTML(;
        canonical="https://jmichel7.github.io/CyclotomicNumbers.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jmichel7/CyclotomicNumbers.jl",
    devbranch="main",
)
