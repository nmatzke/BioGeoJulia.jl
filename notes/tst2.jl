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
include("tst2.jl")

"""

# 
#######################################################


#######################################################
# Run all of the functions here
#######################################################

module Tst2
	include("ModelLikes.jl")
	import .ModelLikes
	#using .Tmp
	
	using DataFrames  # for DataFrame
	using PhyloNetworks
	using BioGeoJulia.TrUtils # for flat2() (similar to unlist)
	using BioGeoJulia.StateSpace
	using BioGeoJulia.TreePass
	using BioGeoJulia.SSEs

	ModelLikes.say_hello2()
	
	#######################################################
	# NEXT:
	# 1. multiply conditional probabilities by lambda
	# 2. load some geographic states, and a phylogeny
	# 3. input into likelihood
	# 4. maximum likelihood (or speed tests)
	#######################################################
	
	ModelLikes.setup_DEC_SSE(2)
	ModelLikes.setup_DEC_SSE(4)
# 	ModelLikes.setup_DEC_SSE(10)
# 	ModelLikes.setup_DEC_SSE(20)
# 	ModelLikes.setup_DEC_SSE(100)
	
	
	
	
	
end # End of module Tst2









