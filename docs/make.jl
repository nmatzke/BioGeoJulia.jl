using Documenter, Example, BioGeoJulia

makedocs(modules = [Example, BioGeoJulia], sitename = "BioGeoJulia.jl")

deploydocs(repo = "github.com/nmatzke/BioGeoJulia.jl.git")
