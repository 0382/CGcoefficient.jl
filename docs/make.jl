using Documenter
push!(LOAD_PATH, "..")
using CGcoefficient

makedocs(
    modules = [CGcoefficient],
    sitename = "CGcoefficient.jl",
    clean = false,
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "Formula" => [
            "Wigner Symbols" => "wigner.md",
            "LS-jj recoupling" => "lsjj.md",
            "Moshinsky" => "moshinsky.md",
        ],
    ]
)

deploydocs(
    repo = "github.com/0382/CGcoefficient.jl.git",
    target = "build/"
)