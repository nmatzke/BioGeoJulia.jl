module Tmp
using BioGeoJulia.TrUtils # for flat2() (similar to unlist)


# (1) List all function names here:
export say_hello, setup_DEC_Cmat, relative_probabilities_of_subsets, relative_probabilities_of_vicariants, discrete_maxent_distrib_of_smaller_daughter_ranges, array_in_array, setup_DEC_Cmat


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

function array_in_array(items, collection)
	TF = collect(repeat([false], length(items)))
	for i in 1:length(items)
		TF[i] = in(items[i], collection)
	end
	if (sum(TF) == length(TF))
		return true
	else
		return false
	end
end


"""
areas_list = [1,2,3]
states_list = areas_list_to_states_list(areas_list, 3, true)
predeclare_array_length=10000000
Carray = setup_DEC_Cmat(areas_list, states_list)

i = 8
j = i
k = i
ancstate = states_list[i]
ancsize = length(ancstate)
lstate = states_list[j]
lsize = length(lstate)
rstate = states_list[k]
rsize = length(rstate)	

ancstate == lstate == rstate

areas_list = [1,2,3]
states_list = areas_list_to_states_list(areas_list, 3, true)
params=(y=1.0,s=1.0,v=1.0,j=0.0)
max_numareas = length(areas_list)
maxent_constraint_01 = 0.0
maxent01symp = relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
maxent01sub = relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
maxent01jump = relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
maxent_constraint_01 = 0.5
maxent01vic = relative_probabilities_of_vicariants(max_numareas, maxent_constraint_01)
maxent01 = (maxent01symp=maxent01symp, maxent01sub=maxent01sub, maxent01vic=maxent01vic, maxent01jump=maxent01jump)
predeclare_array_length=10000000
Carray = setup_DEC_Cmat(areas_list, states_list, params)

"""

function setup_DEC_Cmat(areas_list, states_list, maxent01, params=(y=1.0,s=1.0,v=1.0,j=0.0); predeclare_array_length=10000000)
	numstates = length(states_list)
	
	maxent01symp = maxent01.maxent01symp
	maxent01sub = maxent01.maxent01sub
	maxent01vic = maxent01.maxent01vic
	maxent01jump = maxent01.maxent01jump
	
	y = params.y
	s = params.s
	v = params.v
	j = params.j
	
	# MAKE SURE the states_list is sorted by size
	range_sizes = length.(states_list)
	rangesize_diffs = collect(repeat([0.0], numstates-1))
	for i in 2:length(range_sizes)
	 rangesize_diffs[i-1] = length(states_list[i]) - length(states_list[i-1]) + 0.0
	end
	TF = rangesize_diffs .< 0
	if sum(TF) > 0
		txt = "\nSTOP ERROR in setup_DEC_Cmat(): Your states_list is not ordered in ascending rangesize. Printing states_list:\n"
		print(txt)
		print("\nstates_list")
		print(states_list)
		error(txt)
	end
	
	
	# Pre-declare arrays
	# This contains the sum of the weights for each ancestor row
	row_weightvals = collect(repeat([0.0], numstates))
	
	# Int32 ranges from -32768 to 32767, this suggests a maximum number of states
	Carray_ivals = collect(repeat([0], predeclare_array_length))
	Carray_jvals = collect(repeat([0], predeclare_array_length))
	Carray_kvals = collect(repeat([0], predeclare_array_length))
	Cijk_vals = collect(repeat([0.0], predeclare_array_length))
	Carray_event_types = collect(repeat([""], predeclare_array_length))
	numC = 0 # counter of the number of allow cladogenesis events
	
	
	# Go through:
	# i = ancestor state index
	# j = left state index
	# k = right state index
	for i in 1:numstates
		ancstate = states_list[i]
		ancsize = length(ancstate)
		# Exclude the null range as an ancestor for cladogenetic range-changing events
		# (ancestor of range NULL give 0 likelihood to the observed data)
		if (ancsize == 0)
			continue
		end
		
		# Loop through left and right descendant states
		# Left states
		for j in 1:numstates
			lstate = states_list[j]
			lsize = length(lstate)
			if (lsize == 0)
				continue
			end

			# Right states
			for k in j:numstates # We only have to do half of the possible right events;
													 # reverse each to make a pair
				rstate = states_list[k]
				rsize = length(rstate)
				if (rsize == 0)
					continue
				end
				
				# Sympatry (range-copying)
				if (ancstate == lstate == rstate)
					# Check if the weight > 0.0
					# ancsize, lsize, rsize are the same so we don't have to 
					# choose the smaller daugher
					tmp_weightval = y * maxent01symp[ancsize, lsize] * 1.0 * 1.0
					if tmp_weightval > 0.0
						# Record the range-copying event
						numC += 1
						Carray_event_types[numC] = "y"
						Carray_ivals[numC] = i
						Carray_jvals[numC] = j
						Carray_kvals[numC] = k
						row_weightvals[i] += tmp_weightval
					end # end if tmp_weightval > 0.0
				end # end if (ancstate == lstate == rstate)
				
				# If one of the descendants is identical to the ancestor, 
				# (after we've excluded sympatry)
				# we can have jump dispersal or subset speciation
				if ( (ancstate == rstate) )
				
					# Check for subset sympatry: lstate smaller than rstate, lstate inside rstate
					if ((array_in_array(lstate, rstate) == true) && (lsize < rsize))
						# Check if the weight > 0.0
						# lsize is smaller by definition
						# choose the smaller daughter
						tmp_weightval = s * maxent01sub[ancsize, lsize] * 1.0 * 1.0
						if tmp_weightval > 0.0
							# Record the range-copying event
							numC += 1
							Carray_event_types[numC] = "s"
							Carray_ivals[numC] = i
							Carray_jvals[numC] = j
							Carray_kvals[numC] = k
							row_weightvals[i] += tmp_weightval

							# Same event, flip left/right descendant states
							numC += 1
							Carray_event_types[numC] = "s"
							Carray_ivals[numC] = i
							Carray_jvals[numC] = k
							Carray_kvals[numC] = j
							row_weightvals[i] += tmp_weightval
						end # end if tmp_weightval > 0.0
					end # if ((array_in_array(lstate, rstate) == true) && (lsize < rsize))
					
					# If the left descendant is of size 1, & different from right, then
					# jump dispersal
					if ( (lsize == 1) && (array_in_array(lstate, rstate) == false) )
						tmp_weightval = j * maxent01jump[ancsize, lsize] * 1.0 * 1.0
					end
					
				end # if ( (ancstate == rstate) )
			end # end i (right state indices)
		end # end j (left state indices)
	end # end i (ancestor state indices)
	
	Carray = (Carray_ivals=Carray_ivals, Carray_jvals=Carray_jvals, Carray_kvals=Carray_kvals, Carray_event_types=Carray_event_types)
	
	"""
	# Extract the values
	Carray_ivals = Carray.Carray_ivals;
	Carray_jvals = Carray.Carray_jvals;
	Carray_kvals = Carray.Carray_kvals;
	Carray_event_types = Carray.Carray_event_types;
	hcat(Carray_ivals, Carray_jvals, Carray_kvals, Carray_event_types)
	"""
	
	return Carray
end




end # end module