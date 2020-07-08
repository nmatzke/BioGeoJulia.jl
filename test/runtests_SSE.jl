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


#######################################################
# Do a bunch of tests of the SSE calculation of 
# Ds, Es, and likelihoods, on
# branches, nodes, and trees,
# under a variety of simple and more complex models
#######################################################

@testset "SSEs" begin

#######################################################
# Calculation of Es and Ds on a single branch
# Example BiSSE calculation
# result_EsDs_biSSE_1branch_pureBirth_bl1
# (1 branch, pure birth, no Q transitions, branchlength=1)
#
# Run with:
# source("/GitHub/BioGeoJulia.jl/test/biSSE_1branch_v1.R")
#######################################################
result_EsDs = [1 0.0860322055692215 0.0860322055692215 0 0.739232931655601]











end
