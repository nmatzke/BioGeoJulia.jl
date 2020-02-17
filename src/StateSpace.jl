module StateSpace
using Combinatorics  # for e.g. combinations()
using DataFrames     # for e.g. DataFrame()
using PhyloNetworks
export numstates_from_numareas, areas_list_to_states_list, get_default_inputs, run_model, setup_MuSSE, setup_DEC_DEmat


"""
numstates_from_numareas(3,3,false)
numstates_from_numareas(3,3,true)
numstates_from_numareas(10,1,false)
numstates_from_numareas(10,2,false)
numstates_from_numareas(10,3,false)
numstates_from_numareas(10,10,false)
numstates_from_numareas(10,10,true)
"""


"""
	numstates_from_numareas(numareas[, maxareas, include_null_range])

Calculates the number of possible states in the state space (e.g., the number of possible geographic 
ranges). The inputs are the number of areas (`numareas`), and the maximum range size (`maxareas`).

If `numareas` = `maxareas`, then the size of the state space is 2^`numareas`. This is
because, for each area, a species could be absent (0) or present (1). Therefore, if
there is one area, there are 2x1=2 possible ranges (0 and 1). If there are two areas, there 
are 2x2=4 possible ranges (00, 01, 10, 11). If there are three areas, there are 2x2x2=8 
possible ranges (000, 001, 010, 011, 100, 101, 110, 111), etc.

The `include_null_range` input, if set to `true` allows the "all absent" range (e.g. 000). 
If set to `false`, this range is disallowed, decreasing the size of the state space by 1.

*NOTE:* The size of the state space is a fundamental constraint on the computational speed
of likelihood calculations for biogeographical models using matrix exponentiation, ODEs 
(ordinary differential equations), etc. Researchers must think carefully about how large 
of a state space they require to test their hypotheses, and whether their analysis will 
run fast enough (or at all!).  If `numareas`=20 and `maxareas`=20, the size of the state
space is 1,048,576. Trying to calculate likelihoods for such a huge state space will likely
just fill up computer memory and then crash the program.

Researchers (and reviewers and editors!) should clearly recognize that any computational
inference is inherently a compromise between a very complex reality, and the memory and 
speed limitations of computers. Different people might reach different conclusions about
where exactly this compromise should end up. I tend to think that running somewhat simpler and 
faster models, thus allowing more time to run variant models, model checks, etc., is more 
valuable than setting up the most complex model you can think of, devoting months of 
computing time to it, and then either (a) losing all of that work when it crashes, or (b)
treating the output as gospel because you don't have the time or money available to do
anything else.


# Examples
```julia-repl
julia> numstates_from_numareas(3,3,false)
7

julia> numstates_from_numareas(3,3,true)
8

julia> numstates_from_numareas(10,1,false)
10

julia> numstates_from_numareas(10,2,false)
55

julia> numstates_from_numareas(10,3,false)
175

julia> numstates_from_numareas(10,10,false)
1023

julia> numstates_from_numareas(10,10,true)
1024

julia> numstates_from_numareas(20,20,true)
1048576
```
"""
function numstates_from_numareas(numareas=3, maxareas=1, include_null_range=false)
	# The formula for the number of geographic states, based on the number of areas,
	# is:
	#
	# sum_from_k=1_to_m (N choose k)
	# 
	numstates = 0
	
	# to avoid "does not accept keyword arguments"
	n = numareas + 0
	
	for k in 1:maxareas
		tmp_numstates = binomial(n, k)   # see also Rchoose
		numstates = numstates + tmp_numstates
	end
	
	if include_null_range == true
		numstates = numstates + 1
	end
	return numstates
end



# areas_list (1-based) to states_list (1-based)
"""
areas_list = collect(1:3)
states_list = areas_list_to_states_list(areas_list, 1, false)
states_list = areas_list_to_states_list(areas_list, 1, true)
states_list = areas_list_to_states_list(areas_list, 3, false)
states_list = areas_list_to_states_list(areas_list, 3, true)
areas_list_to_states_list()
"""


"""
	areas_list_to_states_list(areas_list[, maxareas, include_null_range])

Provides the list of possible states (e.g., geographic ranges) in the state space. The 
inputs are:

* `areas_list` - The list of areas. Each area is described with a number. This is done
with the `collect` function, e.g. collect(1:3).

* `maxareas` - The maximum number of areas occupied per geographic range. See
`numstates_from_numareas` for a discussion of how the state space grows (rapidly!) 
with `numareas` and `maxareas`.

* `include_null_range` - if set to `true`, allows the "all absent" range (e.g. 000). 
If set to `false`, this range is disallowed, decreasing the size of the state space by 1.

NOTE: The size of the state space is a fundamental constraint on the computational speed
of likelihood calculations for biogeographical models using matrix exponentiation, ODEs 
(ordinary differential equations), etc. Researchers must think carefully about how large 
of a state space they require to test their hypotheses, and whether their analysis will 
run fast enough (or at all!).  If `numareas`=20 and `maxareas`=20, the size of the state
space is 1,048,576. Trying to calculate likelihoods for such a huge state space will likely
just fill up computer memory and then crash the program.

Researchers (and reviewers and editors!) should clearly recognize that any computational
inference is inherently a compromise between a very complex reality, and the memory and 
speed limitations of computers. Different people might reach different conclusions about
where exactly this compromise should end up. I tend to think that running somewhat simpler and 
faster models, thus allowing more time to run variant models, model checks, etc., is more 
valuable than setting up the most complex model you can think of, devoting months of 
computing time to it, and then either (a) losing all of that work when it crashes, or (b)
treating the output as gospel because you don't have the time or money available to do
anything else.


# Examples
```julia-repl
julia> areas_list = collect(1:3)
3-element Array{Int64,1}:
 1
 2
 3

julia> states_list = areas_list_to_states_list(areas_list, 1, false)
3-element Array{Array{Any,1},1}:
 [1]
 [2]
 [3]

julia> states_list = areas_list_to_states_list(areas_list, 1, true)
4-element Array{Array{Any,1},1}:
 [] 
 [1]
 [2]
 [3]

julia> states_list = areas_list_to_states_list(areas_list, 3, false)
7-element Array{Array{Any,1},1}:
 [1]      
 [2]      
 [3]      
 [1, 2]   
 [1, 3]   
 [2, 3]   
 [1, 2, 3]

julia> states_list = areas_list_to_states_list(areas_list, 3, true)
8-element Array{Array{Any,1},1}:
 []       
 [1]      
 [2]      
 [3]      
 [1, 2]   
 [1, 3]   
 [2, 3]   
 [1, 2, 3]
```
"""
function areas_list_to_states_list(areas_list=collect(1:3), maxareas=3, include_null_range=false)
	
	# Initialize the states_list to the correct size
	numareas = length(areas_list)
	if maxareas > numareas
		maxareas = numareas
	end
	numstates = numstates_from_numareas(numareas, maxareas, include_null_range)
	# Empty list with numstates states
	states_list = repeat([[]], numstates)
	
	# Populate the states_list
	# If include_null_range=true, the first state is NULL, i.e. [] (empty list)
	if include_null_range == true
		state_num = 1
	else
		state_num = 0
	end
	
	# Fill in the states_list
	for k in 1:maxareas
		tmp_states_list = collect(Combinatorics.combinations(areas_list, k))
		for j in 1:length(tmp_states_list)
			state_num = state_num + 1
			states_list[state_num] = tmp_states_list[j]
		end
	end
	
	return states_list
end






#######################################################
# Run models
#######################################################
function get_default_inputs(n=2)
	# Load tree
	great_ape_newick_string = "((chimp:1,human:1):1,gorilla:2);"
	tr = readTopology(great_ape_newick_string)
	rootnodenum = tr.root
	trdf = prt(tr, rootnodenum)
	#trdf
	
	# Set up a simple MuSSE model
	p_Es_v5 = setup_MuSSE(n; birthRate=0.222222, deathRate=0.1, q01=0.01, q10=0.001)
	
	# Solutions to the E vector
	u0_Es = repeat([0.0], 1*n)
	uE = repeat([0.0], n)
	tspan = (0.0, 1.2*trdf[tr.root,:node_age]) # 110% of tree root age

	prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, u0_Es, tspan, p_Es_v5)
	sol_Es_v5 = solve(prob_Es_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9);

	#######################################################
	# Downpass with ClaSSE
	#######################################################
	# Solve for the Ds
	du = repeat([0.0], n)
	u0 = repeat([0.0], n)
	u0[2] = 1.0  # Starting likelihood
	
	
	#tspan = (0.0, 2.0*trdf[tr.root,:node_age]) # 110% of tree root age
	current_nodeIndex = 1
	res = construct_Res(tr, n)

	# Populate the tip likelihoods -- MODIFY THIS
	res.likes_at_each_nodeIndex_branchTop[1] = u0;
	res.likes_at_each_nodeIndex_branchTop[2] = u0;
	res.likes_at_each_nodeIndex_branchTop[4] = u0;

	solver_options = construct_SolverOpt()
	solver_options.solver=Tsit5()
	solver_options.abstol = 1.0e-6
	solver_options.reltol = 1.0e-6

	inputs = (res=res, trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options)
	return inputs
end


function run_model(inputs=default_inputs())
	res = inputs.res
	trdf = inputs.trdf
	p_Ds_v5 = inputs.p_Ds_v5
	solver_options = inputs.solver_options
	(total_calctime_in_sec, iteration_number) = iterative_downpass_nonparallel_ClaSSE_v5!(res, trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10)
	return (total_calctime_in_sec, iteration_number)
	# return res?
end








#######################################################
# Set up models
#######################################################

# Set up a MuSSE (or BiSSE) model with n states
# This default model assumes:
#   - same birthRate & deathRate across all states
#   - transitions only possible to adjacent states
# 
"""
p_Es_v5 = setup_MuSSE(2; birthRate=0.222222, deathRate=0.1, q01=0.01, q10=0.001)
p_Es_v5 = setup_MuSSE(3; birthRate=0.222222, deathRate=0.1, q01=0.01, q10=0.001)
p_Es_v5 = setup_MuSSE(4, birthRate=0.222222, deathRate=0.1, q01=0.01, q10=0.001)

# Anagenetic transition matrix
hcat(p_Es_v5.p_indices.Qarray_ivals, p_Es_v5.p_indices.Qarray_jvals, p_Es_v5.params.Qij_vals)
DataFrame(p_Es_v5.p_TFs.Qi_eq_i)
p_Es_v5.p_TFs.Qi_sub_i
p_Es_v5.p_TFs.Qj_sub_i

# Cladogenetic transition matrix
hcat(p_Es_v5.p_indices.Carray_ivals, p_Es_v5.p_indices.Carray_jvals, p_Es_v5.p_indices.Carray_kvals, p_Es_v5.params.Cijk_vals)
DataFrame(p_Es_v5.p_TFs.Ci_eq_i)
p_Es_v5.p_TFs.Ci_sub_i
p_Es_v5.p_TFs.Cj_sub_i
p_Es_v5.p_TFs.Ck_sub_i
"""

function setup_MuSSE(n=2; birthRate=0.222222, deathRate=0.1, q01=0.01, q10=0.001)
	# Define Qarray - zeros
	Qarray_ivals = collect(1:(n-1))
	Qarray_jvals = collect(2:n)
	Qarray_ivals = vcat(Qarray_ivals, collect(2:n))
	Qarray_jvals = vcat(Qarray_jvals, collect(1:(n-1)))
	Qij_vals = vcat(repeat([q01],(n-1)), repeat([q10],(n-1)))
	
	# Carray: Cladogenetic parameters
	# column 1: state i
	# column 2: state j
	# column 3: state k
	# column 4: lambda_ijk (a parameter that is at least possibly nonzero, under the given model)
	Carray_ivals = collect(1:n)
	Carray_jvals = collect(1:n)
	Carray_kvals = collect(1:n)
	Cijk_vals = repeat([birthRate], n)

# 	Qij_vals[((Qarray_ivals .== 1) .+ (Qarray_jvals .!= 1)) .== 2]
# 	hcat(Qarray_ivals, Qarray_jvals, Qij_vals)
# 	hcat(Carray_ivals, Carray_jvals, Carray_kvals, Cijk_vals)
	
	# Extinction rates
	mu_vals = repeat([deathRate], n)
	
	# Assemble a "params" tuple
	
	# Possibly varying parameters
	params = (mu_vals=mu_vals, Qij_vals=Qij_vals, Cijk_vals=Cijk_vals)

	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	p_indices = (Qarray_ivals=Qarray_ivals, Qarray_jvals=Qarray_jvals, Carray_ivals=Carray_ivals, Carray_jvals=Carray_jvals, Carray_kvals=Carray_kvals)

	# True/False statements by index
	# The calculation of dEi and dDi for state i involves many
	# ==i and !=i operations across Q and C. These only need to be 
	# done once per problem (may or may not save time to 
	# pre-calculate).
	# 
	# Pre-allocating the Carray_ivals .== i, Qarray_jvals[Qarray_ivals .== i
	# Reduces GC (Garbage Collection) from 40% to ~5%
	# 10+ times speed improvement (!)
	Qi_eq_i = Any[]
	Ci_eq_i = Any[]

	Qi_sub_i = Any[]
	Qj_sub_i = Any[]
	
	# These are the (e.g.) j state-indices (left descendant) when the ancestor==state i
	Ci_sub_i = Any[]
	Cj_sub_i = Any[]
	Ck_sub_i = Any[]
	
	# The push! operation may get slow at huge n
	# This will have to change for non-Mk models
	for i in 1:n
		push!(Qi_eq_i, Qarray_ivals .== i)
		push!(Qi_sub_i, Qarray_ivals[Qarray_ivals .== i])
		push!(Qj_sub_i, Qarray_jvals[Qarray_ivals .== i])

		push!(Ci_eq_i, Carray_ivals .== i)
		push!(Ci_sub_i, Carray_ivals[Carray_ivals .== i])
		push!(Cj_sub_i, Carray_jvals[Carray_ivals .== i])
		push!(Ck_sub_i, Carray_kvals[Carray_ivals .== i])
	end
	
	# Inputs to the Es calculation
	p_TFs = (Qi_eq_i=Qi_eq_i, Ci_eq_i=Ci_eq_i, Qi_sub_i=Qi_sub_i, Qj_sub_i=Qj_sub_i, Ci_sub_i=Ci_sub_i, Cj_sub_i=Cj_sub_i, Ck_sub_i=Ck_sub_i)
	p_Es_v5 = (n=n, params=params, p_indices=p_indices, p_TFs=p_TFs)
	
	tmptxt="""
	n = p_Es_v5.n
	params = p_Es_v5.params
	p_indices = p_Es_v5.p_indices
	p_TFs = p_Es_v5.p_TFs
	"""
	
	return(p_Es_v5)
end

# Set up a sparse Qmat for the DEC model
# It will contain references to the parameters
# dmat: a numareas x numareas matrix of "d" values (range-expansion dispersal)
# elist: a numareas array of "e" values (range-contraction, extirpation, "local extinction")
# amat: a numareas x numareas matrix of "a" values
# mult_mat: a numareas x numareas matrix of dispersal multipliers
# e_mult: a numareas list of extirpation multipliers
# exclude_zeros if true, exclude from the matrix 

"""
numareas = 3
areas_list = collect(1:numareas)
states_list = areas_list_to_states_list(areas_list, 1, false)
amat = reshape(collect(1:(numareas^2)), (numareas,numareas))
dmat = reshape(collect(1:(numareas^2)), (numareas,numareas)) ./ 100
elist = repeat([0.123], numstates)
Qmat = setup_DEC_DEmat(areas_list, states_list, )


states_list = areas_list_to_states_list(areas_list, 1, true)

"""

function setup_DEC_DEmat(areas_list, states_list, dmat=reshape(repeat([0.1], numstates),(length(areas_list),length(areas_list))), elist=repeat([0.01],length(areas_list)), amat=reshape(repeat([0.0], numstates),(length(areas_list),length(areas_list))); allowed_event_types=["d","e"])
	# Set up items to iterate over
	numstates = length(states_list)
	statenums = collect(1:numstates)
	range_sizes = length.(states_list)
	#areas_list = sort(unique(flat2(states_list)))
	numareas = length(areas_list)
	
	# Get the approx size of the nonzero rates (for pre-allocation)
	num_e_rates = numstates
	
	# Count the approximate number of d events by
	# iterating counts upwards (assumes states_list is size-ordered)
	#num_d_rates = ceil((numstates^2)/2)
	num_d_rates = 0
	lengths_states_list = length.(states_list)
	for (i in 1:(length(states_list)-1))
		for (j in (i+1):length(states_list))
			if ((lengths_states_list[i]+1) == lengths_states_list[j])
				num_d_rates = num_d_rates + 1
			end
		end
	end
	
	num_a_rates = 2*numstates
	num_nonzero_rates = num_e_rates + num_d_rates + num_a_rates
	
	# Initialize empty arrays
	Qarray_ivals = Array{Int64, num_nonzero_rates}
	Qarray_jvals = Array{Int64, num_nonzero_rates}
	Qij_vals =  Array{Float64, num_nonzero_rates}  # 0-element Array{Any,1}  
	                     # This is populated by calculating through the others
	#base_vals = Array{Float64, num_nonzero_rates}  # base rates -- d, e, a
	#mod_vals = Array{Float64, num_nonzero_rates}  # base modifiers -- e.g., "2" for AB -> ABC
	#mult_vals = Array{Float64, num_nonzero_rates} # multipliers from mult_mat, e_mult 
	# (ie modification by distance, area, etc.)
	event_type_vals = Array{String, num_nonzero_rates}
	index = 0
	
	
	# Events of "a" type: anagenetic range-switching
	# restricted to single-area ranges (!)
	if (setdiff(["a"], allowed_event_types) == ["a"])
		rangesize_eq1_TF = range_sizes .== 1
		statenums_of_size1 = statenums[rangesize_eq1_TF]
		
		for i in 1:(length(statenums_of_size1)-1)			# starting states
			for j in (i+1):length(statenums_of_size1)		# ending states
				statenum_ival = statenums_of_size1[i]
				statenum_jval = statenums_of_size1[j]
				starting_state = states_list[statenum_ival]
				ending_state = states_list[statenum_jval]
				size_i = length(starting_state)
				size_j = length(ending_state)
			
				# "a" events -- anagenetic transitions between single-area states
				# Only use a's, if a>0, or 0 a's are desired
				# Do "amat[i,j][]" instead of "amat[i,j]" in case it's a Ref()
				if (starting_state != ending_state) && (size_i == 1) && (size_j == 1) # single-areas, not null
					starting_areanum = starting_state[] # only works with length 1
					ending_areanum = ending_state[]     # only works with length 1
					if (amat[starting_areanum, ending_areanum][] > 0) || (exclude_zeros == false)
						# Add to arrays
						# Forward event
						index = index + 1
						event_type_vals[index] = "a"
						Qarray_ivals[index] = statenum_ival
						Qarray_jvals[index] = statenum_jval
						Qij_vals[index] = amat[starting_areanum,ending_areanum]
						
						# Reverse event
						index = index + 1
						event_type_vals[index] = "a"
						Qarray_ivals[index] = statenum_jval
						Qarray_jvals[index] = statenum_ival
						Qij_vals[index] = amat[ending_areanum, starting_areanum]

	# 					push!(event_type_vals, "a")
	# 					push!(Qarray_ivals, i)
	# 					push!(Qarray_jvals, j)
	# 					push!(base_vals, amat[starting_areanum,ending_areanum])
	# 					push!(mod_vals, 1)
	# 					push!(mult_vals, 1)
					end # end if (amat
				end # ending if a[] > 0
			end # ending for j in (i+1):length(statenums_of_size1
		end # ending for i in 1:(length(statenums_of_size1)-1)
	end # ending if (setdiff(["a"]
	
	
	# Events of "d" type: anagenetic range-expansion dispersal
	if (setdiff(["d"], allowed_event_types) == ["d"])
		for i in 1:(numstates-1)			# starting states
			for j in (i+1):numstates		# ending states
				starting_state = states_list[i]
				ending_state = states_list[j]
		
				# "d" events -- anagenetic range-expansion events
				# Is the ending range 1 area more than the starting range?
				if (starting_state != ending_state) && ((size_i+1) == size_j) && (size_i != 0) # state i is 1 smaller; not null
					starting_areanums = starting_state
					ending_areanums = ending_state
					end_areanums_not_found_in_start_areas = setdiff(ending_areanums, starting_areanums)
					if length(end_areanums_not_found_in_start_areas) == 1
						# Add to arrays
						index = index + 1
						event_type_vals[index] = "d"
						Qarray_ivals[index] = i
						Qarray_jvals[index] = j
				
						# Add up the d events
						tmp_d_sum = 0.0
						for k in 1:size_i
							tmp_d_sum = tmp_d_sum + dmat[starting_areanums[k], ending_areanums[end_areanums_not_found_in_start_areas]][]
						end

						Qij_vals[index] = tmp_d_sum
	# 					push!(event_type_vals, "d")
	# 					push!(Qarray_ivals, i)
	# 					push!(Qarray_jvals, j)
	# 					push!(base_vals, tmp_d_sum)
	# 					push!(mod_vals, size_i)
	# 					push!(mult_vals, 1)
						
					end # ending if length(end_areanums_not_found_in_start_areas) == 1
				end # ending if (starting_state != ending_state)...
			end # ending j loop
		end # ending i loop
	end # ending if (setdiff(["d"]
		
	# Events of "e" type: anagenetic range-loss/extirpation
	if (setdiff(["e"], allowed_event_types) == ["e"])
		for i in 2:numstates			# starting states
			for j in 1:(i-1)		# ending states
				if (starting_state != ending_state) && ((size_i-1) == size_j) && (size_i != 0) # state i is 1 bigger; not null
					starting_areanums = starting_state
					ending_areanums = ending_state
					start_areanums_not_found_in_end_areas = setdiff(starting_areanums, ending_areanums)
					if length(start_areanums_not_found_in_end_areas) == 1
						# Add to arrays
						index = index + 1
						event_type_vals[index] = "e"
						Qarray_ivals[index] = i
						Qarray_jvals[index] = j
						Qij_vals[index] = elist[start_areanums_not_found_in_end_areas]
					end # ending if length(start_areanums_not...
				end # ending if (starting_state != ending_state)...
			end # ending j loop
		end # ending i loop
	end # ending if (setdiff(["e"]
	
	# Return results
	Qmat = (Qarray_ivals, Qarray_jvals, Qij_vals, event_type_vals)
	
	"""
	Qarray_ivals = Qmat.Qarray_ivals
	Qarray_jvals = Qmat.Qarray_jvals
	Qij_vals = Qmat.Qij_vals
	event_type_vals = Qmat.event_type_vals
	
	hcat(Qarray_ivals, Qarray_jvals, Qij_vals, event_type_vals)
	"""
	
	return Qmat
end







end # Ends the module command
