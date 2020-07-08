module StateSpace
using BioGeoJulia.TrUtils # for e.g. flat2
using Combinatorics  # for e.g. combinations()
using DataFrames     # for e.g. DataFrame()
using PhyloNetworks

#######################################################
# Example use of maximum entropy on discrete case
# https://github.com/JuliaOpt/Convex.jl/issues/64
#######################################################
using Distributions  # for quantile
using Convex				 # for Convex.entropy(), maximize()
using SCS						 # for SCSSolve, solve (maximize(entropy()))


export CparamsStructure, default_Cparams, sumy, sums, sumv, sumj, numstates_from_numareas, areas_list_to_states_list, get_default_inputs, run_model, setup_MuSSE, setup_DEC_DEmat, update_Qij_vals, setup_DEC_Cmat, relative_probabilities_of_subsets, relative_probabilities_of_vicariants, discrete_maxent_distrib_of_smaller_daughter_ranges, array_in_array, is_event_vicariance, setup_DEC_Cmat, update_Cijk_vals





# Cladogenetic parameter weights structure
# "mutable" means you can change the values referred to by the keys
mutable struct CparamsStructure
	y::Float64
	s::Float64
	v::Float64
	j::Float64
end

"""
Cparams = default_Cparams()
"""
function default_Cparams()
	y = 1.0
	s = 1.0
	v = 1.0
	j = 0.0
	Cparams = CparamsStructure(y, s, v, j)
end



function sumy(x)
	sum(x .== "y")
end

function sums(x)
	sum(x .== "s")
end

function sumv(x)
	sum(x .== "v")
end

function sumj(x)
	sum(x .== "j")
end




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
		numstates += 1
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
			state_num += 1
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








#######################################################
# Set up a sparse Qmat for the DEC model
#######################################################
# It will contain references to the parameters
# dmat: a numareas x numareas matrix of "d" values (range-expansion dispersal)
# elist: a numareas array of "e" values (range-contraction, extirpation, "local extinction")
# amat: a numareas x numareas matrix of "a" values
### mult_mat: a numareas x numareas matrix of dispersal multipliers
### e_mult: a numareas list of extirpation multipliers
### exclude_zeros if true, exclude from the matrix 

"""
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
hcat(Qarray_ivals, Qarray_jvals, Qij_vals, event_type_vals)


states_list = areas_list_to_states_list(areas_list, 3, false)
Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["d","e"])
Qarray_ivals = Qmat.Qarray_ivals
Qarray_jvals = Qmat.Qarray_jvals
Qij_vals = Qmat.Qij_vals
event_type_vals = Qmat.event_type_vals
hcat(Qarray_ivals, Qarray_jvals, Qij_vals, event_type_vals)

states_list = areas_list_to_states_list(areas_list, 3, false)
Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["a"])
Qarray_ivals = Qmat.Qarray_ivals
Qarray_jvals = Qmat.Qarray_jvals
Qij_vals = Qmat.Qij_vals
event_type_vals = Qmat.event_type_vals
hcat(Qarray_ivals, Qarray_jvals, Qij_vals, event_type_vals)

"""

# states_list=areas_list_to_states_list(areas_list, length(areas_list), true), dmat=reshape(repeat([0.1], numstates),(length(areas_list),length(areas_list))), elist=repeat([0.01],length(areas_list)), amat=reshape(repeat([0.0], numstates),(length(areas_list),length(areas_list))); allowed_event_types=["d","e"]

function setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["d","e"])
	# Set up items to iterate over
	numstates = length(states_list)
	statenums = collect(1:numstates)
	range_sizes = length.(states_list)
	#areas_list = sort(unique(flat2(states_list)))
	numareas = length(areas_list)
	
	# Get the approx size of the nonzero rates (for pre-allocation)
	if (in("e", allowed_event_types))
		num_e_rates = numstates
	else
		num_e_rates = 0
	end

	
	# Count the approximate number of d events by
	# iterating counts upwards (assumes states_list is size-ordered)
	#num_d_rates = ceil((numstates^2)/2)
	if (in("d", allowed_event_types))
		num_d_rates = 0
		lengths_states_list = length.(states_list)
		for i in 1:(length(states_list)-1)
			for j in (i+1):length(states_list)
				if ((lengths_states_list[i]+1) == lengths_states_list[j])
					num_d_rates += 1
				end
			end
		end
	else
		num_d_rates = 0
	end

	if (in("a", allowed_event_types))
		num_a_rates = 2*numstates
	else
		num_a_rates = 0
	end
	
	num_nonzero_rates = num_e_rates + num_d_rates + num_a_rates
	
	# Initialize empty arrays
	Qarray_ivals = repeat([0], num_nonzero_rates)
	Qarray_jvals = repeat([0], num_nonzero_rates)
	Qij_vals =  repeat([0.0], num_nonzero_rates)  # 0-element Array{Any,1}  
	                     # This is populated by calculating through the others
	#base_vals = Array{Float64, num_nonzero_rates}  # base rates -- d, e, a
	#mod_vals = Array{Float64, num_nonzero_rates}  # base modifiers -- e.g., "2" for AB -> ABC
	#mult_vals = Array{Float64, num_nonzero_rates} # multipliers from mult_mat, e_mult 
	# (ie modification by distance, area, etc.)
	event_type_vals = repeat([""], num_nonzero_rates)
	index = 0
	
	
	# Events of "a" type: anagenetic range-switching
	# restricted to single-area ranges (!)
	if (in("a", allowed_event_types))
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
						index += 1
						event_type_vals[index] = "a"
						Qarray_ivals[index] = statenum_ival
						Qarray_jvals[index] = statenum_jval
						Qij_vals[index] = amat[starting_areanum,ending_areanum]
						
						# Reverse event
						index += 1
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
	end # ending if (in("a", allowed_event_types))
	
	
	# Events of "d" type: anagenetic range-expansion dispersal
	if (in("d", allowed_event_types))
		for i in 1:(numstates-1)			# starting states
			for j in (i+1):numstates		# ending states
				starting_state = states_list[i]
				ending_state = states_list[j]
				size_i = length(starting_state)
				size_j = length(ending_state)
						
				# "d" events -- anagenetic range-expansion events
				# Is the ending range 1 area more than the starting range?
				if (starting_state != ending_state) && ((size_i+1) == size_j) && (size_i != 0) # state i is 1 smaller; not null
					starting_areanums = starting_state
					ending_areanums = ending_state
					end_areanums_not_found_in_start_areas = setdiff(ending_areanums, starting_areanums)
					if length(end_areanums_not_found_in_start_areas) == 1
						# Add to arrays
						index += 1
						event_type_vals[index] = "d"
						Qarray_ivals[index] = i
						Qarray_jvals[index] = j
				
						# Add up the d events
						tmp_d_sum = 0.0
						for k in 1:size_i
							# Because there is only 1 end_areanums_not_found_in_start_areas
							tmp_d_sum += dmat[starting_areanums[k], end_areanums_not_found_in_start_areas[1]][]
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
	end # ending if (in("d", allowed_event_types)
		
	# Events of "e" type: anagenetic range-loss/extirpation
	# NOTE: we could combine with "d" with some effort
	if (in("e", allowed_event_types))
		for i in 2:numstates			# starting states
			for j in 1:(i-1)		# ending states
				starting_state = states_list[i]
				ending_state = states_list[j]
				size_i = length(starting_state)
				size_j = length(ending_state)

				if (starting_state != ending_state) && ((size_i-1) == size_j) && (size_i != 0) # state i is 1 bigger; not null
					starting_areanums = starting_state
					ending_areanums = ending_state
					start_areanums_not_found_in_end_areas = setdiff(starting_areanums, ending_areanums)
					if length(start_areanums_not_found_in_end_areas) == 1
						# Add to arrays
						index += 1
						event_type_vals[index] = "e"
						Qarray_ivals[index] = i
						Qarray_jvals[index] = j
						# Because there is only 1 area in start_areanums_not_found_in_end_areas
						Qij_vals[index] = elist[start_areanums_not_found_in_end_areas[1]]
					end # ending if length(start_areanums_not...
				end # ending if (starting_state != ending_state)...
			end # ending j loop
		end # ending i loop
	end # ending if (in("e", allowed_event_types)
	
	keepTF = event_type_vals .!= ""
	Qarray_ivals = Qarray_ivals[keepTF]
	Qarray_jvals = Qarray_jvals[keepTF]
	Qij_vals = Qij_vals[keepTF]
	event_type_vals = event_type_vals[keepTF]
	
	# Return results
	Qmat = (Qarray_ivals=Qarray_ivals, Qarray_jvals=Qarray_jvals, Qij_vals=Qij_vals, event_type_vals=event_type_vals)
	
	"""
	Qarray_ivals = Qmat.Qarray_ivals
	Qarray_jvals = Qmat.Qarray_jvals
	Qij_vals = Qmat.Qij_vals
	event_type_vals = Qmat.event_type_vals
	
	hcat(Qarray_ivals, Qarray_jvals, Qij_vals, event_type_vals)
	"""
	
	return Qmat
end # end setup_DEC_DEmat()






"""
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
Qmat2 = update_Qij_vals(Qmat, areas_list, states_list, dmat, elist, amat )
Qmat2

Qarray_ivals = Qmat2.Qarray_ivals
Qarray_jvals = Qmat2.Qarray_jvals
Qij_vals = Qmat2.Qij_vals
event_type_vals = Qmat2.event_type_vals
Qmat2_df = hcat(Qarray_ivals, Qarray_jvals, Qij_vals, event_type_vals)

Qmat1_df
Qmat2_df

"""
function update_Qij_vals(Qmat, areas_list, states_list, dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))), elist=repeat([1.0], length(areas_list)), amat=dmat )

	numstates = length(states_list)
	statenums = collect(1:numstates)
	range_sizes = length.(states_list)
	#areas_list = sort(unique(flat2(states_list)))
	numareas = length(areas_list)

	Qarray_ivals = Qmat.Qarray_ivals
	Qarray_jvals = Qmat.Qarray_jvals
	Qij_vals = Qmat.Qij_vals
	event_type_vals = Qmat.event_type_vals

	# Update the "d" events (anagenetic range expansion)
	TF = event_type_vals .== "d"
	if (sum(TF) > 0)
		ivals = Qarray_ivals[TF]
		jvals = Qarray_jvals[TF]
		rates = Qij_vals[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			decstate = states_list[j]
			decsize = length(decstate)
			
			starting_areanums = ancstate
			ending_areanums = decstate
			end_areanums_not_found_in_start_areas = setdiff(ending_areanums, starting_areanums)
			
			# Add up the d events
			tmp_d_sum = 0.0
			for k in 1:ancsize
				# Because there is only 1 end_areanums_not_found_in_start_areas
				tmp_d_sum += dmat[starting_areanums[k], end_areanums_not_found_in_start_areas[1]][]
			end
			# Store the result
			rates[z] = tmp_d_sum
		end
		Qij_vals[TF] = rates
	end # End update of d event weights


	# Update the "a" events (anagenetic range expansion)
	TF = event_type_vals .== "a"
	if (sum(TF) > 0)
		ivals = Qarray_ivals[TF]
		jvals = Qarray_jvals[TF]
		rates = Qij_vals[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			decstate = states_list[j]
			decsize = length(decstate)
			
			starting_areanum = ancstate # because it has size=1 by def.
			ending_areanum = decstate   # because it has size=1 by def.

			# Store the result
			rates[z] = amat[starting_areanum,ending_areanum]
		end
		Qij_vals[TF] = rates
	end # End update of a event weights

	# Update the "e" events (anagenetic range extinction/
	# local extirpation/range loss)
	TF = event_type_vals .== "e"
	if (sum(TF) > 0)
		ivals = Qarray_ivals[TF]
		jvals = Qarray_jvals[TF]
		rates = Qij_vals[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			decstate = states_list[j]
			decsize = length(decstate)
			
			starting_areanums = ancstate
			ending_areanums = decstate
			start_areanums_not_found_in_end_areas = setdiff(starting_areanums, ending_areanums)

			# Store the result
			rates[z] = elist[start_areanums_not_found_in_end_areas[1]]
		end
		Qij_vals[TF] = rates
	end # End update of a event weights

	# Return results
	# Return results
	Qmat2 = (Qarray_ivals=Qarray_ivals, Qarray_jvals=Qarray_jvals, Qij_vals=Qij_vals, event_type_vals=event_type_vals)
	return Qmat2
end # end function update_Qij_vals











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

total_numareas=6
maxent_constraint_01 = 0.0001    # ranges from 0.0001 (all weight on ranges of size 1)
						                     # to 0.9999 (all weight on ranges of max size)
NA_val=NaN
relative_probabilities_of_subsets(total_numareas, maxent_constraint_01, NA_val)
"""

function relative_probabilities_of_subsets(total_numareas=6, maxent_constraint_01=0.5, NA_val=NaN)
	# Set up a matrix to hold the maxent distributions of relative prob of 
	# smaller daughter range sizes
	relprob_subsets_matrix = reshape(repeat([NA_val], total_numareas^2), (total_numareas,total_numareas))
	for i in 1:total_numareas
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



function relative_probabilities_of_vicariants(total_numareas=6, maxent_constraint_01=0.5, NA_val=NaN)
	# Set up a matrix to hold the maxent distributions of relative prob of 
	# smaller daughter range sizes
	relprob_subsets_matrix = reshape(repeat([NA_val], total_numareas^2), (total_numareas,total_numareas))
	relprob_subsets_matrix[1,:] .= NA_val
	
	for i in 2:total_numareas
		ancestor_range_size = i
		tmpstates = collect(1:i)# .+ 0.0
		tmpstates_Floats = collect(1:i) .+ 0.0
		max_smaller_rangesize = median(tmpstates_Floats)
		possible_vicariance_smaller_rangesizes = tmpstates[tmpstates_Floats .< max_smaller_rangesize]
		vic_total_numareas = maximum(possible_vicariance_smaller_rangesizes)

		maxent_result = flat2(discrete_maxent_distrib_of_smaller_daughter_ranges(vic_total_numareas, maxent_constraint_01))

		if (i <= 3)
# 			print("\n")
# 			print(i)
# 
# 			print("\n")
# 			print(vic_total_numareas)
# 			
# 			print("\n")
# 			print(relprob_subsets_matrix[i,1:vic_total_numareas])
# 
# 			print("\n")
# 			print(maxent_result)
# 
# 			print("\n")
			relprob_subsets_matrix[i,1:vic_total_numareas] .= maxent_result
		else
			relprob_subsets_matrix[i,1:vic_total_numareas] = maxent_result
		end
	end
	
	return relprob_subsets_matrix
end




"""
numareas = 6
total_numareas=6
maxent_constraint_01=0.0
discrete_maxent_distrib_of_smaller_daughter_ranges(total_numareas, maxent_constraint_01)

total_numareas=6
maxent_constraint_01=0.5
discrete_maxent_distrib_of_smaller_daughter_ranges(total_numareas, maxent_constraint_01)

total_numareas=6
maxent_constraint_01=1.0
discrete_maxent_distrib_of_smaller_daughter_ranges(total_numareas, maxent_constraint_01)

numareas = 2
total_numareas=6
maxent_constraint_01=0.0
discrete_maxent_distrib_of_smaller_daughter_ranges(total_numareas, maxent_constraint_01)

total_numareas=6
maxent_constraint_01=0.5
discrete_maxent_distrib_of_smaller_daughter_ranges(total_numareas, maxent_constraint_01)

total_numareas=6
maxent_constraint_01=1.0
discrete_maxent_distrib_of_smaller_daughter_ranges(total_numareas, maxent_constraint_01)

"""
function discrete_maxent_distrib_of_smaller_daughter_ranges(total_numareas=6, maxent_constraint_01=0.5)
	n = total_numareas
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
	#sol = Convex.solve!(problem, SCS.Optimizer())
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
ancstate = [1, 2, 3,4];
lstate = [1, 2];
rstate = [4];
is_event_vicariance(ancstate, lstate, rstate)


ancstate = [1, 2, 3,4];
lstate = [1, 2];
rstate = [2, 4];
is_event_vicariance(ancstate, lstate, rstate)

ancstate = [1, 2, 3,4];
lstate = [1, 2];
rstate = [3, 4];
is_event_vicariance(ancstate, lstate, rstate)

"""


function is_event_vicariance(ancstate, lstate, rstate)
	ancsize = length(ancstate)
	lsize = length(lstate)
	rsize = length(rstate)
	
	if (ancsize != lsize+rsize)
		return false
	end
	
	# If any lstate areas are in rstate, then not vicariance
	for i in 1:lsize
		if in(lstate[i], rstate)
			return false
		end
	end
	
	# If any lstate areas not in ancstate, then not vicariance
	for i in 1:lsize
		if in(lstate[i], ancstate) == false
			return false
		end
	end

	# If any rstate areas not in ancstate, then not vicariance
	for i in 1:rsize
		if in(rstate[i], ancstate) == false
			return false
		end
	end
	
	# Otherwise, vicariance
	return true
	
end


"""
# Default maxent multipliers for different kinds of events
include("/GitHub/BioGeoJulia.jl/src/StateSpace.jl")
import .StateSpace

total_numareas = 3
maxent01 = StateSpace.maxent01_defaults(total_numareas);
maxent01symp = maxent01.maxent01symp
maxent01sub = maxent01.maxent01sub
maxent01vic = maxent01.maxent01vic
maxent01jump = maxent01.maxent01jump


total_numareas = 4
maxent01 = StateSpace.maxent01_defaults(total_numareas);
maxent01symp = maxent01.maxent01symp
maxent01sub = maxent01.maxent01sub
maxent01vic = maxent01.maxent01vic
maxent01jump = maxent01.maxent01jump

total_numareas = 4
maxent01 = StateSpace.maxent01_defaults(total_numareas; maxent_constraint_01v=0.5);
maxent01symp = maxent01.maxent01symp
maxent01sub = maxent01.maxent01sub
maxent01vic = maxent01.maxent01vic
maxent01jump = maxent01.maxent01jump

"""
function maxent01_defaults(total_numareas=4; maxent_constraint_01=0.0, maxent_constraint_01v=0.0)
	maxent_constraint_01 = 0.0
	maxent01symp = relative_probabilities_of_subsets(total_numareas, maxent_constraint_01)
	maxent01sub = relative_probabilities_of_subsets(total_numareas, maxent_constraint_01)
	maxent01jump = relative_probabilities_of_subsets(total_numareas, maxent_constraint_01)
	maxent_constraint_01v = 0.0
	maxent01vic = relative_probabilities_of_vicariants(total_numareas, maxent_constraint_01v)
	maxent01 = (maxent01symp=maxent01symp, maxent01sub=maxent01sub, maxent01vic=maxent01vic, maxent01jump=maxent01jump)
	return maxent01
end


"""
# Print a Carray to a data.frame
"""
function prtC(Carray)
	Cdf = DataFrame(event=Carray.Carray_event_types, i=Carray.Carray_ivals, j=Carray.Carray_jvals, k=Carray.Carray_kvals, wt=Carray.Cijk_weights, val=Carray.Cijk_vals)
	return Cdf
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
Cparams=(y=1.0,s=1.0,v=1.0,j=0.0)
total_numareas = length(areas_list)
maxent_constraint_01 = 0.0
maxent01symp = relative_probabilities_of_subsets(total_numareas, maxent_constraint_01)
maxent01sub = relative_probabilities_of_subsets(total_numareas, maxent_constraint_01)
maxent01jump = relative_probabilities_of_subsets(total_numareas, maxent_constraint_01)
maxent_constraint_01 = 0.5
maxent01vic = relative_probabilities_of_vicariants(total_numareas, maxent_constraint_01)
maxent01 = (maxent01symp=maxent01symp, maxent01sub=maxent01sub, maxent01vic=maxent01vic, maxent01jump=maxent01jump)
predeclare_array_length=10000000
Carray = setup_DEC_Cmat(areas_list, states_list, maxent01, Cparams)

"""

function setup_DEC_Cmat(areas_list, states_list, maxent01=NaN, Cparams=default_Cparams(), dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))); predeclare_array_length=Integer(min(length(states_list)*length(states_list)*round((length(states_list)/2)), 10000000)), min_precision=1e-9)
	numareas = length(areas_list)
	total_numareas = numareas
	numstates = length(states_list)
	
	if (isnan(maxent01) == true)
		maxent01 = maxent01_defaults(total_numareas)
	end
	
	maxent01symp = maxent01.maxent01symp
	maxent01sub = maxent01.maxent01sub
	maxent01vic = maxent01.maxent01vic
	maxent01jump = maxent01.maxent01jump
	
	# Weights
	y_wt = Cparams.y
	s_wt = Cparams.s
	v_wt = Cparams.v
	j_wt = Cparams.j
	
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
	Cijk_weights = collect(repeat([0.0], predeclare_array_length))
	Carray_event_types = collect(repeat([""], predeclare_array_length))
	numC = 0 # counter of the number of allow cladogenesis events
	
	
	# Go through:
	# i = ancestor state index
	# j = left state index
	# k = right state index
	
	# Preallocate this vector ONCE, size = numareas * 2
	#tmp_merged_vec = repeat([0], 2*numareas)
	
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
				continue # go to the next j
			end
		
			# Right states
			for k in 1:numstates # We only have to do half of the possible right events;
													 # reverse each to make a pair
				rstate = states_list[k]
				rsize = length(rstate)
				if (rsize == 0)
					continue # go to the next k
				end
				
				if (y_wt > min_precision)
					# Sympatry (range-copying)
# 204:
					if (all([lsize == ancsize, rsize==ancsize, lstate==ancstate, rstate==ancstate]))

#2244					if (all([lstate==ancstate, rstate==ancstate]))
						# Check if the weight > 0.0
						# ancsize, lsize, rsize are the same so we don't have to 
						# choose the smaller daugher
						tmp_weightval = y_wt * maxent01symp[ancsize, lsize] * 1.0 * 1.0
						if tmp_weightval > min_precision
							# Record the range-copying event
							numC += 1
							Carray_event_types[numC] = "y"
							Carray_ivals[numC] = i
							Carray_jvals[numC] = j
							Carray_kvals[numC] = k
							Cijk_weights[numC] = tmp_weightval
							row_weightvals[i] += tmp_weightval
							continue
						end # end if tmp_weightval > 0.0
					end # end if (all([lsize == ancsize, rsize==ancsize, lstate==ancstate, rstate==ancstate])
				end # end if (y_wt > min_precision)
				
				# If one of the descendants is identical to the ancestor, 
				# (after we've excluded sympatry)
				# we can have jump dispersal or subset speciation
				if ( all([ancsize==rsize, ancstate==rstate]) )
					# Subset sympatry
					if (s_wt > min_precision)
						# Check for subset sympatry: lstate smaller than rstate, lstate inside rstate
						if ((lsize < rsize) && (array_in_array(lstate, rstate) == true) )
							# Check if the weight > 0.0
							# lsize is smaller by definition
							# choose the smaller daughter
							tmp_weightval = s_wt * maxent01sub[ancsize, lsize] * 1.0 * 1.0
							if tmp_weightval > min_precision
								# Record the range-copying event
								numC += 1
								Carray_event_types[numC] = "s"
								Carray_ivals[numC] = i
								Carray_jvals[numC] = j
								Carray_kvals[numC] = k
								Cijk_weights[numC] = tmp_weightval
								row_weightvals[i] += tmp_weightval

								# Same event, flip left/right descendant states
								numC += 1
								Carray_event_types[numC] = "s"
								Carray_ivals[numC] = i
								Carray_jvals[numC] = k
								Carray_kvals[numC] = j
								Cijk_weights[numC] = tmp_weightval
								row_weightvals[i] += tmp_weightval
								continue
							end # end if tmp_weightval > 0.0
						end # end if ((array_in_array(lstate, rstate) == true) && (lsize < rsize))
					end # end if (s_wt > min_precision)
					
					# Jump dispersal
					if (j_wt > min_precision)
						# If the left descendant is of size 1, & different from right, then
						# jump dispersal
						if ( (lsize == 1) && (array_in_array(lstate, rstate) == false) )
							# Historically, on analogy with other DEC cladogenesis events
							# the weight of each j event was not influenced by the number of
							# ranges of the ancestor. That is followed here. Obviously, 
							# other models could be imagined! (e.g., d events *are* influenced
							# by ancestor rangesize).
							try_jump_dispersal_based_on_dist = true
							normalize_by_number_of_dispersal_events = true
							jweight_for_cell_based_on_distances = 0.0
							if (try_jump_dispersal_based_on_dist == true)
								for anc_area in ancstate
									for left_area in lstate
										 jweight_for_cell_based_on_distances += dmat[anc_area,left_area]
									end
								end
								# Normalize by number of possible jump dispersals
								if (normalize_by_number_of_dispersal_events == true)
									jweight_for_cell_based_on_distances = jweight_for_cell_based_on_distances / (ancsize * lsize)
								end
							else
								# 
								jweight_for_cell_based_on_distances = 1.0
							end # end if (try_jump_dispersal_based_on_dist == true)
					
							# Calculate the final weight of this jump dispersal
							tmp_weightval = j_wt * maxent01jump[ancsize, lsize] * 1.0 * 1.0 * jweight_for_cell_based_on_distances
					
							print("\n")
							print([i, j, k, tmp_weightval])
							if (tmp_weightval > min_precision)
								print("yes")
								# Record the jump-dispersal event
								numC += 1
								Carray_event_types[numC] = "j"
								Carray_ivals[numC] = i
								Carray_jvals[numC] = j
								Carray_kvals[numC] = k
								Cijk_weights[numC] = tmp_weightval
								row_weightvals[i] += tmp_weightval

								# Same event, flip left/right descendant states
								numC += 1
								Carray_event_types[numC] = "j"
								Carray_ivals[numC] = i
								Carray_jvals[numC] = k
								Carray_kvals[numC] = j
								Cijk_weights[numC] = tmp_weightval
								row_weightvals[i] += tmp_weightval
								continue
							end # if (tmp_weightval > 0.0)
							# end of jump dispersal
						end # end if ( (lsize == 1) && (array_in_array(lstate, rstate) == false) )
					end # end if (j_wt > min_precision)
				end # end if ( (ancstate == rstate) )
			
				# Vicariance
				if (v_wt > min_precision)
					# Check if the combined vector equals the ancestor vector					
					#tmp_merged_vec .= repeat([0], 2*numareas)
					#combined_size = length(lstate)+length(rstate)
					
					#tmp_merged_vec[1:length(lstate)] = lstate
					#tmp_merged_vec[(length(lstate)+1):(length(lstate)+length(rstate))] = rstate
					#combined_vector = sort(tmp_merged_vec)
					if ( is_event_vicariance(ancstate, lstate, rstate) )
						smaller_range_size = min(lsize, rsize)
						tmp_weightval = v_wt * maxent01vic[ancsize,smaller_range_size] * 1.0 * 1.0
						if (tmp_weightval > min_precision)
							# Record the jump-dispersal event
							numC += 1
							Carray_event_types[numC] = "v"
							Carray_ivals[numC] = i
							Carray_jvals[numC] = j
							Carray_kvals[numC] = k
							Cijk_weights[numC] = tmp_weightval
							row_weightvals[i] += tmp_weightval
							continue
							# Same event, flip left/right descendant states
							# You won't hit it again, as k >= i
# 							numC += 1
# 							Carray_event_types[numC] = "v"
# 							Carray_ivals[numC] = i
# 							Carray_jvals[numC] = k
# 							Carray_kvals[numC] = j
# 							Cijk_weights[numC] = tmp_weightval
# 							row_weightvals[i] += tmp_weightval
						end # end if (tmp_weightval > 0.0)
					end # end if ( combined_vector == sort(ancstate) )
				end # end if ( (v>0.0) && ..
				# End of Vicariance
			end # end i (right state indices)
		end # end j (left state indices)
	end # end i (ancestor state indices)
	
	TF = Carray_event_types .!= ""
	Carray_ivals = Carray_ivals[TF]
	Carray_jvals = Carray_jvals[TF]
	Carray_kvals = Carray_kvals[TF]
	Cijk_weights = Cijk_weights[TF]
	Carray_event_types = Carray_event_types[TF]
	
	# Convert the weights to conditional event probabilities
	num_clado_events = length(Cijk_weights)
	Cijk_vals = collect(repeat([0.0], num_clado_events))
	for i in 1:length(states_list)
		TF = Carray_ivals .== i
		Cijk_vals[TF] = Cijk_weights[TF] ./ row_weightvals[i]
	end

	
	Carray = (Carray_event_types=Carray_event_types, Carray_ivals=Carray_ivals, Carray_jvals=Carray_jvals, Carray_kvals=Carray_kvals, Cijk_weights=Cijk_weights, Cijk_vals=Cijk_vals, row_weightvals=row_weightvals)
	
	"""
	# Extract the values
	Carray_event_types = Carray.Carray_event_types;
	Carray_ivals = Carray.Carray_ivals;
	Carray_jvals = Carray.Carray_jvals;
	Carray_kvals = Carray.Carray_kvals;
	Cijk_weights = Carray.Cijk_weights;
	Cijk_vals = Carray.Cijk_vals;
	row_weightvals = Carray.row_weightvals
	DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals)
	row_weightvals
	"""
	
	return Carray
end # end setup_DEC_Cmat()



#######################################################
# Update the Cijk_vals
#######################################################
function update_Cijk_vals(Carray, areas_list, states_list, maxent01, Cparams=default_Cparams(), dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))) )

	Carray_ivals = Carray.Carray_ivals;
	Carray_jvals = Carray.Carray_jvals;
	Carray_kvals = Carray.Carray_kvals;
	Carray_event_types = Carray.Carray_event_types;
	Cijk_weights = Carray.Cijk_weights;
	Cijk_vals = Carray.Cijk_vals;
	row_weightvals = Carray.row_weightvals

	numstates = length(states_list)
	
	maxent01symp = maxent01.maxent01symp
	maxent01sub = maxent01.maxent01sub
	maxent01vic = maxent01.maxent01vic
	maxent01jump = maxent01.maxent01jump
	
	# Weights
	y_wt = Cparams.y
	s_wt = Cparams.s
	v_wt = Cparams.v
	j_wt = Cparams.j
	
	# Update the "y" events (narrow sympatry / range-copying)
	TF = Carray_event_types .== "y"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			weights[z] = y_wt * maxent01symp[ancsize, lsize] * 1.0 * 1.0
		end
		Cijk_weights[TF] = weights
	end # End update of y event weights


	# Update the "s" events (subset sympatry)
	TF = Carray_event_types .== "s"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			weights[z] = s_wt * maxent01sub[ancsize, lsize] * 1.0 * 1.0
		end
		Cijk_weights[TF] = weights
	end # End update of s event weights

	# Update the "v" events (vicariance)
	TF = Carray_event_types .== "v"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			smaller_range_size = min(lsize, rsize)
			weights[z] = v_wt * maxent01vic[ancsize,smaller_range_size] * 1.0 * 1.0
		end
		Cijk_weights[TF] = weights
	end # End update of s event weights
	
	# Update the "j" events (jump dispersal)
	TF = Carray_event_types .== "v"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			
			# j events, modified by distance / multipliers (via input "dmat") if needed
			try_jump_dispersal_based_on_dist = true
			normalize_by_number_of_dispersal_events = true
			jweight_for_cell_based_on_distances = 0.0
			if (try_jump_dispersal_based_on_dist == true)
				for anc_area in ancstate
					for left_area in lstate
						 jweight_for_cell_based_on_distances += dmat[anc_area,left_area]
					end
				end
				# Normalize by number of possible jump dispersals
				if (normalize_by_number_of_dispersal_events == true)
					jweight_for_cell_based_on_distances = jweight_for_cell_based_on_distances / (ancsize * lsize)
				end
			else
				# 
				jweight_for_cell_based_on_distances = 1.0
			end # end if (try_jump_dispersal_based_on_dist == true)
	
			# Calculate the final weight of this jump dispersal
			tmp_weightval = j_wt * maxent01jump[ancsize, lsize] * 1.0 * 1.0 * jweight_for_cell_based_on_distances
			weights[z] = tmp_weightval
		end
		Cijk_weights[TF] = weights
	end # End update of s event weights

	df1 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
	
	row_weightvals_df = by(df1, :i, :weight => sum)
	row_weightvals = row_weightvals_df[!,2]
	
	# If i=1 is missing from row_weightvals_df, add it to row_weightvals
	if (in(1, row_weightvals_df[!,1]) == false)
		row_weightvals = repeat([0], 1+length(row_weightvals_df[!,2]))
		row_weightvals[1] = 1
		row_weightvals[2:length(row_weightvals)] = row_weightvals_df[!,2]
	end
	row_weightvals
	
	# Convert the weights to conditional event probabilities
	num_clado_events = length(Cijk_weights)
	Cijk_vals = collect(repeat([0.0], num_clado_events))
	for i in 1:length(states_list)
		TF = Carray_ivals .== i
		Cijk_vals[TF] = Cijk_weights[TF] ./ row_weightvals[i]
	end

	df2 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
	row_weightvals_df = by(df2, :i, :weight => sum)
	
	# Finally, return updated Carray:
	Carray = (Carray_event_types=Carray_event_types, Carray_ivals=Carray_ivals, Carray_jvals=Carray_jvals, Carray_kvals=Carray_kvals, Cijk_weights=Cijk_weights, Cijk_vals=Cijk_vals, row_weightvals=row_weightvals)

	"""
	# Extract the values
	Carray_event_types = Carray.Carray_event_types;
	Carray_ivals = Carray.Carray_ivals;
	Carray_jvals = Carray.Carray_jvals;
	Carray_kvals = Carray.Carray_kvals;
	Cijk_weights = Carray.Cijk_weights;
	Cijk_vals = Carray.Cijk_vals;
	row_weightvals = Carray.row_weightvals;
	DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals)
	row_weightvals
	"""

	return Carray
end # end update_Cijk_vals()












end # Ends the module command
