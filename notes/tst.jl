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
include("tst.jl")

"""

# 
#######################################################


#######################################################
# Run all of the functions here
#######################################################

module Tst
	include("Tmp.jl")
	import .Tmp
	#using .Tmp

	using BioGeoJulia.StateSpace 


	Tmp.say_hello()
	# say_hello()



	# Distribution of smaller daughters:
	numareas=6
	max_numareas=6
	maxent_constraint_01=0.0
	maxent_result = Tmp.discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas, maxent_constraint_01)
	
	max_numareas=6
	maxent_constraint_01=0.5
	maxent_result = Tmp.discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas, maxent_constraint_01)
	
	max_numareas=6
	maxent_constraint_01=1.0
	maxent_result = Tmp.discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas, maxent_constraint_01)

	max_numareas=6
	maxent_constraint_01 = 0.0001    # ranges from 0.0001 (all weight on ranges of size 1)
																	 # to 0.9999 (all weight on ranges of max size)
	NA_val=NaN
	maxent_constraint_01=0.0
	relprob_subsets_matrix = Tmp.relative_probabilities_of_subsets(max_numareas, maxent_constraint_01, NA_val)

	maxent_constraint_01=0.5
	relprob_subsets_matrix = Tmp.relative_probabilities_of_subsets(max_numareas, maxent_constraint_01, NA_val)

	maxent_constraint_01=1.0
	relprob_subsets_matrix = Tmp.relative_probabilities_of_subsets(max_numareas, maxent_constraint_01, NA_val)


	maxent_constraint_01=0.5
	relprob_subsets_matrix = Tmp.relative_probabilities_of_vicariants(max_numareas, maxent_constraint_01, NA_val)

	areas_list = [1,2,3]
	states_list = areas_list_to_states_list(areas_list, 3, true)
	params=(y=1.0,s=1.0,v=1.0,j=0.0)
	max_numareas = length(areas_list)
	maxent_constraint_01 = 0.0
	maxent01symp = Tmp.relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
	maxent01sub = Tmp.relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
	maxent01jump = Tmp.relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
	maxent_constraint_01 = 0.5
	maxent01vic = Tmp.relative_probabilities_of_vicariants(max_numareas, maxent_constraint_01)
	maxent01 = (maxent01symp=maxent01symp, maxent01sub=maxent01sub, maxent01vic=maxent01vic, maxent01jump=maxent01jump)
	predeclare_array_length=10000000
	Carray = Tmp.setup_DEC_Cmat(areas_list, states_list, maxent01, params)

	# your other test code here


end