using Test, BioGeoJulia, Combinatorics, DataFrames

# List each BioGeoJulia code file prefix here
using BioGeoJulia.Example
using BioGeoJulia.StateSpace
using BioGeoJulia.TreePass
using BioGeoJulia.PhyDFutils

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
#	@test prt
end

@testset "TreePass" begin
#	@test prt
end