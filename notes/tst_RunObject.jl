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
	
	
	# construct the sub-objects, then the full one
	lgdata_fn = "/GitHub/BioGeoJulia.jl/ex/homs1/homs_geog.txt"
	geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

	tr_bkup = readTopology("(((human:0.5):0.5,(chimp:0.5):0.5):1.0,(gorilla:1.5):0.5);")
	trfn = "/GitHub/BioGeoJulia.jl/ex/homs1/homs_w_breakpoint.newick"
	tr = readTopology(trfn)
	trdf = prt(tr)

	ro = RunObject.construct_RunObject(trfn, lgdata_fn)

	
end # End of module