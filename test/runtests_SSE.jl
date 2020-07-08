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
include("/GitHub/BioGeoJulia.jl/notes/ModelLikes.jl")
import .ModelLikes
inputs = ModelLikes.setup_MuSSE_biogeo()
Rnames(inputs)
Rnames(inputs.p_Ds_v5.p_indices)
Rcbind(inputs.p_Ds_v5.p_indices.Qarray_ivals, inputs.p_Ds_v5.p_indices.Qarray_jvals, inputs.p_Ds_v5.params.Qij_vals, inputs.p_Ds_v5.params.Qarray_event_types)
Rcbind(inputs.p_Ds_v5.p_indices.Carray_ivals, inputs.p_Ds_v5.p_indices.Carray_jvals, inputs.p_Ds_v5.p_indices.Carray_kvals, inputs.p_Ds_v5.params.Cijk_vals)


result_EsDs = [1 0.0860322055692215 0.0860322055692215 0 0.739232931655601]

# Repeat calculation in Julia
in_params = (birthRate=0.2, deathRate=0.0, d_val=0.0, e_val=0.0, a_val=0.0, j_val=0.0)
inputs = ModelLikes.setup_MuSSE_biogeo(numareas=2, tr=readTopology("((chimp:1,human:1):1,gorilla:2);"); root_age_mult=1.5, max_range_size=1, include_null_range=true, in_params=in_params)
Es_interpolator = inputs.p_Ds_v5.sol_Es_v5

# Initialize evolving state vector
u = repeat([0.0], 2*n)

# Starting values
du = repeat([0.0], 2*n)
u0 = repeat([0.0], 2*n)
u0[n+1] = 1.0
tspan = (0.0, 10.0)

prob = ODEProblem(parameterized_ClaSSE, u0, tspan, p)
sol = solve(prob, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)










end
