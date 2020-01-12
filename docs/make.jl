using Documenter, BioGeoJulia

makedocs(modules = [BioGeoJulia], sitename = "BioGeoJulia.jl")

deploydocs(repo = "github.com/nmatzke/BioGeoJulia.jl.git")
