#######################################################
# Likelihood calculations: input tree, tipdata, model,
#                          get likelihood back
#######################################################

module ModelLikes

using BenchmarkTools # for @time
using InvertedIndices # for Not
using Distributed
using Random					# for MersenneTwister()
using Dates						# for e.g. DateTime, Dates.now()
using PhyloNetworks

using DataFrames          # for DataFrame()
using BioGeoJulia.TrUtils # for flat2() (similar to unlist)
using BioGeoJulia.StateSpace
using BioGeoJulia.TreePass
using BioGeoJulia.SSEs

# (1) List all function names here:
export say_hello2

#######################################################
# Temporary file to store functions under development
#
# Start with:
# 
# Setup:

"""
cd("/GitHub/BioGeoJulia.jl/notes/")
include("tst2.jl")
"""
#######################################################


#######################################################
# (2) write the functions here
#######################################################

say_hello2() = println("Hello dude2!")


function setup_DEC_SSE(tr=readTopology("((chimp:1,human:1):1,gorilla:2);"), numareas=2)
	areas_list = [1:numareas]
	states_list = areas_list_to_states_list(areas_list, numareas, false)
	Cparams = default_Cparams()
	max_numareas = length(areas_list)
	maxent_constraint_01 = 0.0
	maxent01symp = relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
	maxent01sub = relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
	maxent01jump = relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
	maxent_constraint_01 = 0.0
	maxent01vic = relative_probabilities_of_vicariants(max_numareas, maxent_constraint_01)
	maxent01 = (maxent01symp=maxent01symp, maxent01sub=maxent01sub, maxent01vic=maxent01vic, maxent01jump=maxent01jump)

	Carray = setup_DEC_Cmat(areas_list, states_list, maxent01, Cparams)

	
	
	res = construct_Res(tr, n)
	rootnodenum = tr.root
	trdf = prt(tr, rootnodenum)
	
	birthRate = 0.222222
	deathRate = 0.1
	
	d_val = 0.01
	e_val = 0.001
	j_val = 0.0
	
	

end



function calclike_DEC_SSE(tr, tipdata)


end


end # end module ModelLikes

