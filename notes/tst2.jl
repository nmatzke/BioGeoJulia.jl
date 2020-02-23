#######################################################
# A test module to run Tmp.jl, not in scope Main
# 
# Development workflow recommended by:
#
# https://docs.julialang.org/en/v1/manual/workflow-tips/
#
# Setup:

"""
cd("/GitHub/BioGeoJulia.jl/notes/")
include("tst2.jl")

"""

# 
#######################################################


#######################################################
# Run all of the functions here
#######################################################

module Tst2
	include("ModelLikes.jl")
	import .ModelLikes
	#using .Tmp
	
	using DataFrames  # for DataFrame
	using PhyloNetworks
	using BioGeoJulia.TrUtils # for flat2() (similar to unlist)
	using BioGeoJulia.StateSpace
	using BioGeoJulia.TreePass
	using BioGeoJulia.SSEs

	ModelLikes.say_hello2()
	
	#######################################################
	# NEXT:
	# 1. multiply conditional probabilities by lambda
	# 2. load some geographic states, and a phylogeny
	# 3. input into likelihood
	# 4. maximum likelihood (or speed tests)
	#######################################################
	
	ModelLikes.setup_DEC_SSE(2)
#	ModelLikes.setup_DEC_SSE(4)
# 	ModelLikes.setup_DEC_SSE(10)
# 	ModelLikes.setup_DEC_SSE(20)
# 	ModelLikes.setup_DEC_SSE(100)
	
	# Update Qij_vals
	numareas = 3
	areas_list = collect(1:numareas)
	states_list = areas_list_to_states_list(areas_list, 3, true)
	numstates = length(states_list)
	amat = reshape(collect(1:(numareas^2)), (numareas,numareas))
	dmat = reshape(collect(1:(numareas^2)), (numareas,numareas)) ./ 100
	elist = repeat([0.123], numstates)
	allowed_event_types=["d","e"]

	Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["d","e"])
	Qarray_ivals = Qmat.Qarray_ivals
	Qarray_jvals = Qmat.Qarray_jvals
	Qij_vals = Qmat.Qij_vals
	event_type_vals = Qmat.event_type_vals
	Qmat1_df = hcat(Qarray_ivals, Qarray_jvals, Qij_vals, event_type_vals)

	# Update!
	dmat = reshape(repeat([0.5], numareas^2), (numareas,numareas))
	Qmat2 = ModelLikes.update_Qij_vals(Qmat, areas_list, states_list, dmat, elist, amat )
	Qmat2
	
	Qarray_ivals = Qmat2.Qarray_ivals
	Qarray_jvals = Qmat2.Qarray_jvals
	Qij_vals = Qmat2.Qij_vals
	event_type_vals = Qmat2.event_type_vals
	Qmat2_df = hcat(Qarray_ivals, Qarray_jvals, Qij_vals, event_type_vals)
	
	Qmat1_df
	Qmat2_df
end # End of module Tst2









