using Documenter
using DyadicKDE

cp("../README.md", "src/index.md", force=true)
makedocs(sitename="DyadicKDE.jl")
