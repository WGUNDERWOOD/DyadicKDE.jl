#push!(LOAD_PATH,"../src/")

using Documenter
using DyadicKDE

makedocs(sitename="DyadicKDE.jl",
         modules=DyadicKDE,
         pages=["Home" => "index.md",
                "Documentation" => "documentation.md"])

deploydocs(repo="github.com/WGUNDERWOOD/DyadicKDE.jl.git")
