
#######################################################
# Load dependencies of BioGeoJulia
#######################################################
import Pkg
using Pkg

# Add these packages
Pkg.add("Combinatorics")


# Then import or using as needed
import Combinatorics.combinations


#######################################################
# Load the BioGeoJulia package
#######################################################
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

# Try some functions
numstates_from_numareas(3,3,false)
numstates_from_numareas(3,3,true)
numstates_from_numareas(10,1,false)
numstates_from_numareas(10,2,false)
numstates_from_numareas(10,3,false)
numstates_from_numareas(10,10,false)
numstates_from_numareas(10,10,true)

area_nums = collect(1:3)
states_list = areas_list_to_states_list(area_nums, 1, false)
states_list = areas_list_to_states_list(area_nums, 1, true)
states_list = areas_list_to_states_list(area_nums, 3, false)
states_list = areas_list_to_states_list(area_nums, 3, true)
areas_list_to_states_list()


