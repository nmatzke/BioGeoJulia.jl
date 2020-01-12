using Documenter, BGJ_Example, BioGeoJulia

makedocs(modules = [BGJ_Example, BioGeoJulia], sitename = "BioGeoJulia.jl")

deploydocs(repo = "github.com/nmatzke/BioGeoJulia.jl.git")
