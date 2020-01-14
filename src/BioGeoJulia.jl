module BioGeoJulia
export hello_BioGeoJulia, add_one_BioGeoJulia

# List each BioGeoJulia code file here
include("Example.jl")			# default examples
include("TrUtils.jl")			# basic utility functions for trees, prt() tree tables (DFs), etc.
include("StateSpace.jl")	# set up lists of areas and states (geographic ranges)
include("SSEs.jl")				# SSE calculations with various amounts of speed optimization
include("TreePass.jl")		# downpass and uppass through the phylogeny; prt() etc.


#######################################################
# Put functions here
#######################################################
"""
    hello(who::String)

Return "BioGeoJulia says, hi `who`".
"""
hello_BioGeoJulia(who::String) = "BioGeoJulia says, hi $who"

"""
    add_one_BioGeoJulia(x::Number)

Return `x + 1`.
"""
add_one_BioGeoJulia(x::Number) = x + 1

end # Ends the module command
