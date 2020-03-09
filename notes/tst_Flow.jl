#######################################################
# A test module to run Flow.jl, not in scope Main
# 
# Development workflow recommended by:
#
# https://docs.julialang.org/en/v1/manual/workflow-tips/
#
# Setup:

"""
cd("/GitHub/BioGeoJulia.jl/notes/")
include("tst_Flow.jl")

"""

module Tst_Flow
	include("ModelLikes.jl")
	import .ModelLikes
	#using .Tmp

	include("Flow.jl")
	import .Flow


	inputs = ModelLikes.setup_DEC_SSE(numareas=2, tr=readTopology("((chimp:1,human:1):1,gorilla:2);"))
	
	
end

