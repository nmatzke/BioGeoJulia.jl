
#######################################################
# Load the BioGeoJulia package
#######################################################
import Pkg
using Pkg
Pkg.rm("BioGeoJulia")
Pkg.add(PackageSpec(path="/GitHub/BioGeoJulia.jl"))
using BioGeoJulia

# Run the tests directory
Pkg.test("BioGeoJulia")

# Re-run the tests directory
# First, commit to Master (quick), then:
Pkg.add(PackageSpec(path="/GitHub/BioGeoJulia.jl"))
using BioGeoJulia
Pkg.test("BioGeoJulia")



include("/drives/Dropbox/_njm/__julia/julia4Rppl_v3.jl")

area_nums = collect(1:3)
states_list = areas_list_to_states_list(area_nums, 1, false)
states_list = areas_list_to_states_list(area_nums, 1, true)
states_list = areas_list_to_states_list(area_nums, 3, false)
states_list = areas_list_to_states_list(area_nums, 3, true)
areas_list_to_states_list()