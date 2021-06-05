using Documenter

push!(LOAD_PATH, "../src/")

using CGcoefficient

makedocs(
    sitename="CGcoefficient",
    pages = [
        "Home" => "index.md",
        "Formula" => "formula.md",
        "API" => "api.md"
    ]
)