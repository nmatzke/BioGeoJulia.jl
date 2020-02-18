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

	# your other test code here


end