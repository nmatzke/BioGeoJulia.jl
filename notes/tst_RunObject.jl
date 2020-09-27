#######################################################
# A test module to run RunObject.jl, not in scope Main
# 
# Development workRunObject recommended by:
#
# https://docs.julialang.org/en/v1/manual/workRunObject-tips/
#
# Setup:

"""
cd("/GitHub/BioGeoJulia.jl/notes/")
include("tst_RunObject.jl")

"""

module Tst_RunObject
	include("RunObject.jl")
	import .RunObject

	using BioGeoJulia.StateSpace 
	using DataFrames  # for DataFrame

	RunObject.say_hello4()
	
end # End of module