using Test, BioGeoJulia, DataFrames

using Dates									# for e.g. Dates.now(), DateTime
using PhyloNetworks					# most maintained, emphasize; for HybridNetwork
using Distributed						# for e.g. @spawn
using Combinatorics					# for e.g. combinations()
using DataFrames

# List each BioGeoJulia code file prefix here
using BioGeoJulia.Example
using BioGeoJulia.StateSpace
using BioGeoJulia.TreePass
using BioGeoJulia.TrUtils
using BioGeoJulia.SSEs

@testset "Example" begin
	@test hello("Julia") == "Hello, Julia"
	@test domath(2.0) ≈ 7.0
end


@testset "BioGeoJulia" begin
	@test hello_BioGeoJulia("Julia") == "BioGeoJulia says, hi Julia"
	@test add_one_BioGeoJulia(2.0) ≈ 3.0   # note the approximate equals (compare ≈=)
end

@testset "StateSpace" begin
	@test numstates_from_numareas(3,3,false) == 7
	@test numstates_from_numareas(3,3,true) == 8
	@test numstates_from_numareas(10,1,false) == 10
	@test numstates_from_numareas(10,2,false) == 55
	@test numstates_from_numareas(10,3,false) == 175
	@test numstates_from_numareas(10,10,false) == 1023
	@test numstates_from_numareas(10,10,true) == 1024
	@test numstates_from_numareas(20,20,true) == 1048576
	
	# Set up list of areas
	area_nums = collect(1:3)
	
	tmpstr = "Array{Any,1}[[1], [2], [3]]"
	states_list_answer = eval(Meta.parse(tmpstr))
	states_list = areas_list_to_states_list(area_nums, 1, false)
	@test states_list == states_list_answer

	tmpstr = "Array{Any,1}[[], [1], [2], [3]]"
	states_list_answer = eval(Meta.parse(tmpstr))
	states_list = areas_list_to_states_list(area_nums, 1, true)
	@test states_list == states_list_answer

	tmpstr = "Array{Any,1}[[1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]"
	states_list_answer = eval(Meta.parse(tmpstr))
	states_list = areas_list_to_states_list(area_nums, 3, false)
	@test states_list == states_list_answer
	
	tmpstr = "Array{Any,1}[[], [1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]"
	states_list_answer = eval(Meta.parse(tmpstr))
	states_list = areas_list_to_states_list(area_nums, 3, true)
	@test states_list == states_list_answer

	# still need to be done:
	"""
	@test default_Cparams() 
	Cparams == CparamsStructure(1.0, 1.0, 1.0, 0.0)
	# default_Cparams

	note I'm getting 'Cparams not defined' which makes no sense at that
	literally in the definition of the function? the output of default_cparams() is
	still CparamsStructure(1.0, 1.0, 1.0, 0.0) (however when this is tested, still comes back false
	when compared to CparamsStructure(1.0, 1.0, 1.0, 0.0), so not sure what's going on there)
	
	# sumy
	# sums
	# sumv
	# sumj
	unsure what these do
	"""
	# get_default_inputs prt not defined?
	# run_model
	
	tmpstr = "(mu_vals = [0.1, 0.1], Qij_vals = [0.01, 0.001], Cijk_vals = [0.222222, 0.222222])"
	p_Es_v5 = setup_MuSSE(2; birthRate=0.222222, deathRate=0.1, q01=0.01, q10=0.001)
	p_Es_v5_tmp =  eval(Meta.parse(tmpstr))
	@test p_Es_v5.params == p_Es_v5_tmp
	
	tmpstr = "(mu_vals = [0.1, 0.1], Qij_vals = [0.01, 0.001], Cijk_vals = [0.333333, 0.333333])"
	p_Es_v5 = setup_MuSSE(2; birthRate=0.333333, deathRate=0.1, q01=0.01, q10=0.001)
	p_Es_v5_tmp =  eval(Meta.parse(tmpstr))
	@test p_Es_v5.params == p_Es_v5_tmp

	tmpstr = "(mu_vals = [0.2, 0.2], Qij_vals = [0.01, 0.001], Cijk_vals = [0.222222, 0.222222])"
	p_Es_v5 = setup_MuSSE(2; birthRate=0.222222, deathRate=0.2, q01=0.01, q10=0.001)
	p_Es_v5_tmp =  eval(Meta.parse(tmpstr))
	@test p_Es_v5.params == p_Es_v5_tmp

	tmpstr = "(mu_vals = [0.1, 0.1], Qij_vals = [0.2, 0.001], Cijk_vals = [0.222222, 0.222222])"
	p_Es_v5 = setup_MuSSE(2; birthRate=0.222222, deathRate=0.1, q01=0.2, q10=0.001)
	p_Es_v5_tmp =  eval(Meta.parse(tmpstr))
	@test p_Es_v5.params == p_Es_v5_tmp

	tmpstr = "(mu_vals = [0.1, 0.1], Qij_vals = [0.01, 0.2], Cijk_vals = [0.222222, 0.222222])"
	p_Es_v5 = setup_MuSSE(2; birthRate=0.222222, deathRate=0.1, q01=0.01, q10=0.2)
	p_Es_v5_tmp =  eval(Meta.parse(tmpstr))
	@test p_Es_v5.params == p_Es_v5_tmp

	# setup_DEC_DEmat
	numareas = 3
	areas_list = collect(1:numareas)
	states_list = areas_list_to_states_list(areas_list, 3, true)
	numstates = length(states_list)
	amat = reshape(collect(1:(numareas^2)), (numareas,numareas))
	dmat = reshape(collect(1:(numareas^2)), (numareas,numareas)) ./ 100
	elist = repeat([0.123], numstates)
	allowed_event_types=["d","e"]

	Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat)
	@test Qmat.Qarray_jvals[4] == 7

	states_list = areas_list_to_states_list(areas_list, 3, false)
	Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["d","e"])
	@test Qmat.Qarray_jvals[4] == 6

	states_list = areas_list_to_states_list(areas_list, 3, false)
	Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["a"])
	@test Qmat.Qarray_jvals[4] == 1

	# update_Qij_vals
	Qmat1_df = hcat(Qarray_ivals, Qarray_jvals, Qij_vals, event_type_vals)
	dmat = reshape(repeat([0.5], numareas^2), (numareas,numareas))
	"""
	Qmat2 = update_Qij_vals(Qmat, areas_list, states_list, dmat, elist, amat )
	ERROR: MethodError: Cannot `convert` an object of type Array{Int64,2} to an object of type Float64
	I SEE NO CONVERT-LIKE MOVES WITHIN THIS FUNCTION, I AM UNSURE WHERE IT IS TRYING TO CONVERT ANYTHING
	"""
	# relative_probabilities_of_subsets
	# relative_probabilities_of_vicariants()
	# discrete_maxent_distrib_of_smaller_daughter_ranges()
	"""
	relative_probabilities_of_subsets()
	ERROR: ArgumentError: MathProgBase solvers like `solve!(problem, SCSSolver())` are no longer supported. Use instead e.g. `solve!(problem, SCS.Optimizer)`.
	
	Same issue for 
	relative_probabilities_of_vicariants()
	discrete_maxent_distrib_of_smaller_daughter_ranges()
	
	"""

	# array_in_array
	

	# is_event_vicariance
	ancstate = [1, 2, 3,4];
	lstate = [1, 2];
	rstate = [4];
	@test is_event_vicariance(ancstate, lstate, rstate) == false

	ancstate = [1, 2, 3,4];
	lstate = [1, 2];
	rstate = [2, 4];
	@test is_event_vicariance(ancstate, lstate, rstate) == false

	ancstate = [1, 2, 3,4];
	lstate = [1, 2];
	rstate = [3, 4];
	@test is_event_vicariance(ancstate, lstate, rstate) == true

	# setup_DEC_Cmat
	"""
	areas_list = [1,2,3]
	states_list = areas_list_to_states_list(areas_list, 3, true)
	Cparams=(y=1.0,s=1.0,v=1.0,j=0.0)
	max_numareas = length(areas_list)
	maxent_constraint_01 = 0.0
	maxent01symp = relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
	maxent01sub = relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
	maxent01jump = relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
	maxent_constraint_01 = 0.5
	maxent01vic = relative_probabilities_of_vicariants(max_numareas, maxent_constraint_01)
	maxent01 = (maxent01symp=maxent01symp, maxent01sub=maxent01sub, maxent01vic=maxent01vic, maxent01jump=maxent01jump)
	predeclare_array_length=10000000
	Carray = setup_DEC_Cmat(areas_list, states_list, Cparams)

	Cannot be run because relative_probabilities_of_subsets cannot be run!
	"""

	# update_Cijk_vals
end

@testset "TrUtils" begin
	tmpstr = "HybridNetwork"
	answer = eval(Meta.parse(tmpstr))
	great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
	tr = readTopology(great_ape_newick_string)
	@test type(tr) == answer
	
	setwd("/Users/")
	@test setwd("/Users/") == cd("/Users/")
	@test getwd() == pwd()
	@test getwd() == "/Users"
	@test Rgetwd() == pwd()
	@test Rgetwd() == "/Users"

	tmparray = recursive_find("/GitHub/BioGeoJulia.jl")
	@test tmparray[1] == "/GitHub/BioGeoJulia.jl/src/BioGeoJulia.jl"

	tmparray = include_jls("/GitHub/BioGeoJulia.jl")
	@test tmparray[1] == "/GitHub/BioGeoJulia.jl/src/BioGeoJulia.jl"

	A = ones(3,3)
	B = ones(3,3)
	@test dim(A) == size(A)
	@test Rdim(A) == size(A)

	C = Int64[1,2,3,4,5,6,7,8,9,10]
	@test seq(1, 10, 1) == C
	@test Rchoose(10,5) == 252

	D = ones(3,6)
	@test Rcbind(A, B) == hcat(A,B)
	@test Rcbind(A, B) == D
	E = ones(6,3)
	@test Rrbind(A, B) == vcat(A,B)
	@test Rrbind(A, B) == E

	@test paste(1:12, delim="") == "123456789101112"
	@test paste0(1:12) == "123456789101112"

	F = "tester"
	@test type(F) == String
	@test class(F) == "String"
	@test Rclass(F) == "String"
	
	#is there a reason slashslash() has the internal code repeated several times?
	@test slashslash("//GitHub/BioGeoJulia.jl//src//BioGeoJulia.jl") == "/GitHub/BioGeoJulia.jl/src/BioGeoJulia.jl"
	@test addslash("/GitHub/BioGeoJulia.jl/src/BioGeoJulia.jl") == "/GitHub/BioGeoJulia.jl/src/BioGeoJulia.jl/"

	tmpmatrix = [3 1; 3 2; 5 3; 5 4; 7 5; 7 6]
	tmpstr = repr(tmpmatrix)
	tmpstr2 = eval(Meta.parse(tmpstr))
	@test Reval(tmpstr) == tmpstr2

	tmpstr = "Array{Any,1}[[1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]"
	tmpstr2 = tmpstr
	states_list = Reval(tmpstr)
	@test Rdput(states_list) == tmpstr2

	great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
	tr = readTopology(great_ape_newick_string)
	@test Rnames(tr)[1] == :numTaxa
	@test Rtypes(tr)[1] == Int64
	@test ont(tr) == Any[:numTaxa Int64; :numNodes Int64; :numEdges Int64; :node Array{PhyloNetworks.Node,1}; :edge Array{PhyloNetworks.Edge,1}; :root Int64; :names Array{String,1}; :hybrid Array{PhyloNetworks.Node,1}; :numHybrids Int64; :cladewiseorder_nodeIndex Array{Int64,1}; :visited Array{Bool,1}; :edges_changed Array{PhyloNetworks.Edge,1}; :nodes_changed Array{PhyloNetworks.Node,1}; :leaf Array{PhyloNetworks.Node,1}; :ht Array{Float64,1}; :numht Array{Int64,1}; :numBad Int64; :hasVeryBadTriangle Bool; :index Array{Int64,1}; :loglik Float64; :blacklist Array{Int64,1}; :partition Array{PhyloNetworks.Partition,1}; :cleaned Bool; :isRooted Bool]


	@test Rnrow(A) == 3
	@test Rncol(A) == 3
	@test Rsize(A) == (3,3)

	tmpDF = DataFrame(A = 1:4, B = ["M", "F", "F", "M"], C = 5:8)
	tmparray = [1,2,3,4]
	@test Rorder(tmpDF) == tmparray
	
	# tmpDF2 = DataFrame(A = 1:4, C = 5:8)
	# headLR(tmpDF, 1, 1) == tmpDF2
	# false?
	
"""
	A note! Even when the dataframes are identical, 
	pulling the left and right collumns seems to twist them?
	In this case there are ONLY 2 collumns? so the outcome should be identical to itself?

	Something to do with it being flattened? 
	Unsure on how to test atm

	tmptmp = DataFrame(A = 1:4, B = ["M", "F", "F", "M"])
	tmptmp2 = DataFrame(A = 1:4, B = ["M", "F", "F", "M"])
	headLR(tmptmp, 1, 1) == tmptmp2

	FALSE
"""
	
	tmparray = (1)
	@test single_element_array_to_scalar(tmparray) == 1

	tmpline = print("1 module BioGeoJulia")
	@test headf("/GitHub/BioGeoJulia.jl/src/BioGeoJulia.jl"; numlines=1) == tmpline

	# How to test these?
	# @test df_to_Rdata
	# @test source("/GitHub/BioGeoJulia.jl/src/BioGeoJulia.jl") == include("/GitHub/BioGeoJulia.jl/src/BioGeoJulia.jl")
	# @test headLR(df, num_startcols=4, num_endcols=4) ==
	# @test flat2(arr) ==
	# @test moref(fn) ==
	# @test scr2str ==

end

@testset "TreePass" begin
	# Test if the printed tree table from prt()
	# gets the node ages correct
	tmpstr = "[0.0, 0.0, 6.0, 0.0, 7.0, 0.0, 12.0]"
	node_ages_in_prt = eval(Meta.parse(tmpstr))
	
	great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
	tr = readTopology(great_ape_newick_string)
	rootnodenum = tr.root
	trdf = prt(tr, rootnodenum)
	@test trdf[!, :node_age] == node_ages_in_prt
end


@testset "SSEs" begin
	# parameterized_ClaSSE
	# parameterized_ClaSSE_Es
	# parameterized_ClaSSE_Ds
	# parameterized_ClaSSE_v5
	# parameterized_ClaSSE_Es_v5
	# parameterized_ClaSSE_Ds_v5

end
