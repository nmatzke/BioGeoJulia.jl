module Tmp
using BioGeoJulia.TrUtils # for flat2() (similar to unlist)


# (1) List all function names here:
export say_hello, setup_DEC_Cmat, relative_probabilities_of_subsets, relative_probabilities_of_vicariants, discrete_maxent_distrib_of_smaller_daughter_ranges


#######################################################
# Temporary file to store functions under development
#
# Start with:
# 
# Setup:

"""
cd("/GitHub/BioGeoJulia.jl/notes/")
include("tst.jl")
"""
#######################################################


#######################################################
# (2) write the functions here
#######################################################

say_hello() = println("Hello dude!")



#######################################################
# Example use of maximum entropy on discrete case
# https://github.com/JuliaOpt/Convex.jl/issues/64
#######################################################
using Distributions  # for quantile
using Convex				 # for Convex.entropy(), maximize()
using SCS						 # for SCSSolve, solve (maximize(entropy()))



"""
######################################################
Weight the possible range-sizes of the smaller daughter range

E.g., if you have subset sympatry from ancestor ABC,
the smaller daugher could be size 1, 2, or 3.

DEC says the weights of these would be 1 0 0
BayArea says the weights would be 0 0 1

With this function, input 
maxent_constraint_01=0.0001 = 1 0 0 
maxent_constraint_01=0.5    = 1 1 1 
maxent_constraint_01=0.9999 = 0 0 1 
######################################################

max_numareas=6
maxent_constraint_01 = 0.0001    # ranges from 0.0001 (all weight on ranges of size 1)
						                     # to 0.9999 (all weight on ranges of max size)
NA_val=NaN
relative_probabilities_of_subsets(max_numareas, maxent_constraint_01, NA_val)
"""

function relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01=0.5, NA_val=NaN)
	# Set up a matrix to hold the maxent distributions of relative prob of 
	# smaller daughter range sizes
	relprob_subsets_matrix = reshape(repeat([NA_val], max_numareas^2), (max_numareas,max_numareas))
	for i in 1:max_numareas
		ancestor_range_size = i
		maxent_result = discrete_maxent_distrib_of_smaller_daughter_ranges(ancestor_range_size, maxent_constraint_01)

		if (i == 1)
			relprob_subsets_matrix[i,1:i] .= maxent_result
		else
			relprob_subsets_matrix[i,1:i] = maxent_result
		end
	end
	
	return relprob_subsets_matrix
end



function relative_probabilities_of_vicariants(max_numareas=6, maxent_constraint_01=0.5, NA_val=NaN)
	# Set up a matrix to hold the maxent distributions of relative prob of 
	# smaller daughter range sizes
	relprob_subsets_matrix = reshape(repeat([NA_val], max_numareas^2), (max_numareas,max_numareas))
	relprob_subsets_matrix[1,:] .= NA_val
	
	for i in 2:max_numareas
		ancestor_range_size = i
		tmpstates = collect(1:i)# .+ 0.0
		tmpstates_Floats = collect(1:i) .+ 0.0
		max_smaller_rangesize = median(tmpstates_Floats)
		possible_vicariance_smaller_rangesizes = tmpstates[tmpstates_Floats .< max_smaller_rangesize]
		vic_max_numareas = maximum(possible_vicariance_smaller_rangesizes)

		maxent_result = flat2(discrete_maxent_distrib_of_smaller_daughter_ranges(vic_max_numareas, maxent_constraint_01))

		if (i <= 3)
# 			print("\n")
# 			print(i)
# 
# 			print("\n")
# 			print(vic_max_numareas)
# 			
# 			print("\n")
# 			print(relprob_subsets_matrix[i,1:vic_max_numareas])
# 
# 			print("\n")
# 			print(maxent_result)
# 
# 			print("\n")
			relprob_subsets_matrix[i,1:vic_max_numareas] .= maxent_result
		else
			relprob_subsets_matrix[i,1:vic_max_numareas] = maxent_result
		end
	end
	
	return relprob_subsets_matrix
end




"""
numareas = 6
max_numareas=6
maxent_constraint_01=0.0
discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas, maxent_constraint_01)

max_numareas=6
maxent_constraint_01=0.5
discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas, maxent_constraint_01)

max_numareas=6
maxent_constraint_01=1.0
discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas, maxent_constraint_01)

numareas = 2
max_numareas=6
maxent_constraint_01=0.0
discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas, maxent_constraint_01)

max_numareas=6
maxent_constraint_01=0.5
discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas, maxent_constraint_01)

max_numareas=6
maxent_constraint_01=1.0
discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas, maxent_constraint_01)

"""
function discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas=6, maxent_constraint_01=0.5)
	n = max_numareas
	x = 1:n

	#discrete_values_padded = cat(0, collect(x)[:], n+1; dims=1)
	discrete_values_padded = collect(x)
	maxent_constraint = quantile(discrete_values_padded, maxent_constraint_01)
	#maxent_constraint = 1

	# Vector of probabilities that must sum to 1
	p = Variable(length(x))
	probability_contraints = [0 <= p, p <= 1, sum(p) == 1];
	feature_constraints = sum(p'*x) == maxent_constraint

	# base or prior (e.g. Uniform)
	#h = pdf.(Uniform(1, n), x)
	
	# This solution updates p (look in p.values)
	problem = maximize(Convex.entropy(p), 0 <= p, p <= 1, sum(p) == 1, feature_constraints)
	sol = Convex.solve!(problem, SCS.SCSSolver(verbose=0))
	maxent_result = abs.(round.(p.value; digits=3))
	return maxent_result
end


#######################################################
# Set up a Cmat (cladogenesis event-weights matrix)
#######################################################
function setup_DEC_Cmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["d","e"])
	x=1


end




end # end module