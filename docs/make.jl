using Documenter
push!(LOAD_PATH, "../src/")
using CGcoefficient

makedocs(
    modules = [CGcoefficient],
    sitename = "CGcoefficient.jl",
    clean = false,
    pages = [
        "Home" => "index.md",
        "Formula" => "formula.md",
        "API" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/0382/CGcoefficient.jl.git",
    target = "build/"
)