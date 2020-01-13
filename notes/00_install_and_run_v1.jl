
#######################################################
# Add dependencies to BioGeoJulia package development
# https://discourse.julialang.org/t/how-to-manage-dependencies-of-developed-packages/25481/2
#######################################################
# To get your package to have correct-ish [deps] fields
# etc. in Project.toml and Manifest.toml:

cd("/GitHub/BioGeoJulia.jl")
Pkg.activate(".")
Pkg.add("Combinatorics")
Pkg.add("DataFrames")		# for DataFrame
Pkg.add("Dates"	)				# for e.g. Dates.now(), DateTime
Pkg.add("Distributed")  # for e.g. @spawn
Pkg.add("Random")  # for MersenneTwister()

# BioSequences before PhyloNetworks
# https://github.com/BioJulia/BioSequences.jl
]
registry add https://github.com/BioJulia/BioJuliaRegistry.git
add BioSequences
add PhyloNetworks

# de-activate with blank "activate"
Pkg.activate()
# ^C

#######################################################
# NOTE: You may also have to then manually go into
# Project.toml to lower the minimum version numbers of
# some packages.
#
# Commit to GitHub, then continue with below
#######################################################



#######################################################
# Load dependencies of BioGeoJulia (if needed)
#######################################################
import Pkg
using Pkg

# Add these packages
#Pkg.add("Combinatorics")
Pkg.add(Pkg.PackageSpec(;name="Combinatorics", version="1.0.0"))

# Then import or using as needed
import Combinatorics.combinations

Pkg.resolve()






#######################################################
# FROM FRESH JULIA: Load the BioGeoJulia package
#######################################################
import Pkg
using Pkg
Pkg.rm("BioGeoJulia")
Pkg.add(PackageSpec(path="/GitHub/BioGeoJulia.jl"))
using BioGeoJulia

# Run the tests directory
Pkg.test("BioGeoJulia")


#######################################################
# Re-run the tests directory
#######################################################
# First, commit to Master (quick), then:
Pkg.add(PackageSpec(path="/GitHub/BioGeoJulia.jl"))
using BioGeoJulia
Pkg.test("BioGeoJulia")


#######################################################
# Try some functions!
#######################################################
include("/drives/Dropbox/_njm/__julia/julia4Rppl_v3.jl")

# Try some functions
using BioGeoJulia.TrUtils
using BioGeoJulia.StateSpace
using BioGeoJulia.TreePass


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


using PhyloNetworks
using Random					# for MersenneTwister()
#using PhyloPlots


great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr

rootnodenum = tr.root
trdf = prt(tr, rootnodenum)
trdf


# Downpass
indexNum_table = get_nodeIndex_PNnumber(tr)

global res = construct_Res(tr);
current_nodeIndex = res.root_nodeIndex

res.likes_at_each_nodeIndex_branchTop
res.likes_at_each_nodeIndex_branchBot
res.thread_for_each_nodeOp
res.thread_for_each_branchOp
res.node_state

countloop_num_iterations = 10000000
y = countloop(countloop_num_iterations, 1)
@time y = countloop(countloop_num_iterations, 1)


# First time: compilation
# iterative_downpass! -- "!" means the function modifies its arguments

start_compilation = Dates.now()
calctime_in_sec1 = iterative_downpass_nonparallel!(res, max_iterations=Inf, num_iterations=countloop_num_iterations)
end_compilation = Dates.now()
compilation_time = (end_compilation-start_compilation).value / 1000








