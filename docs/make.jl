using Pkg
using TOML
using Documenter
using DocumenterCitations
using DocumenterInterLinks
using QuantumControlExamples: QuantumControlExamples, deploy_examples, collect_examples

PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
VERSION = PROJECT_TOML["version"]
NAME = PROJECT_TOML["name"]
AUTHORS = join(PROJECT_TOML["authors"], ", ") * " and contributors"
GITHUB = "https://github.com/JuliaQuantumControl/QuantumControlExamples.jl"

println("Starting makedocs")

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:numeric)


mkpath(joinpath(@__DIR__, "src", "examples"))
EXAMPLES = [
    example["title"] => joinpath("examples", splitext(example["file"])[1], "index.md")
    for (_, example) in collect_examples(joinpath(@__DIR__, "src", "examples"))
]
@show EXAMPLES

mkpath(joinpath(@__DIR__, "src", "tutorials"))
TUTORIALS = [
    tutorial["title"] =>
        joinpath("tutorials", splitext(tutorial["file"])[1], "index.md") for
    (_, tutorial) in collect_examples(joinpath(@__DIR__, "src", "tutorials"))
]
@show TUTORIALS

PAGES = [
    "Home" => "index.md",
    "Tutorials" => [
        "tutorials.md",
        [hide(title => index_md) for (title, index_md) in TUTORIALS]...
    ],
    "Examples" =>
        ["examples.md", [hide(title => index_md) for (title, index_md) in EXAMPLES]...],
    "References" => "references.md",
]


makedocs(;
    plugins=[bib],
    authors=AUTHORS,
    sitename="Examples",
    doctest=false,  # we have no doctests (but trying to run them is slow)
    format=Documenter.HTML(;
        prettyurls=true,
        canonical="https://juliaquantumcontrol.github.io/QuantumControlExamples.jl",
        assets=[
            "assets/citations.css",
            "assets/custom.css",
            "assets/custom.js",
            asset(
                "https://juliaquantumcontrol.github.io/QuantumControl.jl/dev/assets/topbar/topbar.css"
            ),
            asset(
                "https://juliaquantumcontrol.github.io/QuantumControl.jl/dev/assets/topbar/topbar.js"
            ),
        ],
        mathengine=KaTeX(),
        footer="[JuliaQuantumControl Tutorials and Examples]($GITHUB) powered by [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) and [Literate.jl](https://github.com/fredrikekre/Literate.jl)."
    ),
    pages=PAGES,
    warnonly=true,
)

println("Finished makedocs")

include("deploy.jl")
