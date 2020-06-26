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
	
	tmparray = (1)
	@test single_element_array_to_scalar(tmparray) == 1
	
	# How to test these?
	# @test df_to_Rdata
	# @test source("/GitHub/BioGeoJulia.jl/src/BioGeoJulia.jl") == include("/GitHub/BioGeoJulia.jl/src/BioGeoJulia.jl")
	# @test Rorder(A) ==
	# @test headLR(df, num_startcols=4, num_endcols=4) ==
	# @test flat2(arr) ==
	# @test single_element_array_to_scalar(tmparray) ==  

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

end
