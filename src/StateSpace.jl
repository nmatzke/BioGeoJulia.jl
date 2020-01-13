module StateSpace
using Combinatorics  # for e.g. combinations()
export numstates_from_numareas, areas_list_to_states_list


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
there is one area, there are 2*1=2 possible ranges (0 and 1). If there are two areas, there 
are 2*2=4 possible ranges (00, 01, 10, 11). If there are three areas, there are 2*2*2=8 
possible ranges (000, 001, 010, 011, 100, 101, 110, 111), etc.

The `include_null_range` input, if set to `true` allows the "all absent" range (e.g. 000). 
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
area_nums = collect(1:3)
states_list = areas_list_to_states_list(area_nums, 1, false)
states_list = areas_list_to_states_list(area_nums, 1, true)
states_list = areas_list_to_states_list(area_nums, 3, false)
states_list = areas_list_to_states_list(area_nums, 3, true)
areas_list_to_states_list()
"""


"""
	areas_list_to_states_list(area_nums[, maxareas, include_null_range])

Provides the list of possible states (e.g., geographic ranges) in the state space. The 
inputs are:

* `area_nums` - The list of areas. Each area is described with a number. This is done
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
julia> area_nums = collect(1:3)
3-element Array{Int64,1}:
 1
 2
 3

julia> states_list = areas_list_to_states_list(area_nums, 1, false)
3-element Array{Array{Any,1},1}:
 [1]
 [2]
 [3]

julia> states_list = areas_list_to_states_list(area_nums, 1, true)
4-element Array{Array{Any,1},1}:
 [] 
 [1]
 [2]
 [3]

julia> states_list = areas_list_to_states_list(area_nums, 3, false)
7-element Array{Array{Any,1},1}:
 [1]      
 [2]      
 [3]      
 [1, 2]   
 [1, 3]   
 [2, 3]   
 [1, 2, 3]

julia> states_list = areas_list_to_states_list(area_nums, 3, true)
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
function areas_list_to_states_list(area_nums=collect(1:3), maxareas=3, include_null_range=false)
	
	# Initialize the states_list to the correct size
	numareas = length(area_nums)
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
		tmp_states_list = collect(Combinatorics.combinations(area_nums, k))
		for j in 1:length(tmp_states_list)
			state_num = state_num + 1
			states_list[state_num] = tmp_states_list[j]
		end
	end
	
	return states_list
end



end # Ends the module command
