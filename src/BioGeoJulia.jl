module BioGeoJulia
export hello_BioGeoJulia, add_one_BioGeoJulia

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

end # Ends the module or export command?
