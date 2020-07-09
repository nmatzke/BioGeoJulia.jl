using Test, BioGeoJulia, DataFrames

using Dates									# for e.g. Dates.now(), DateTime
using PhyloNetworks					# most maintained, emphasize; for HybridNetwork
using Distributed						# for e.g. @spawn
using Combinatorics					# for e.g. combinations()
using DataFrames

using LinearAlgebra  # for "I" in: Matrix{Float64}(I, 2, 2)
										 # https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using DataFrames  # for DataFrame
using DifferentialEquations
using OrdinaryDiffEq, Sundials, DiffEqDevTools, Plots, ODEInterfaceDiffEq, ODE, LSODA


# List each BioGeoJulia code file prefix here
using BioGeoJulia.Example
using BioGeoJulia.StateSpace
using BioGeoJulia.TreePass
using BioGeoJulia.TrUtils
using BioGeoJulia.SSEs


"""
# Run with:
include("/GitHub/BioGeoJulia.jl/test/runtests_BiSSE_tree_n1.jl")
"""

@testset "Example" begin
	@test hello("runtests_BiSSE_tree_n1.jl") == "Hello, runtests_BiSSE_tree_n1.jl"
#	@test domath(2.0) ≈ 7.0
end


#######################################################
# Do a bunch of tests of the SSE calculation of 
# Ds, Es, and likelihoods, on
# branches, nodes, and trees,
# under a variety of simple and more complex models
#######################################################

@testset "runtests_BiSSE_tree_n1.jl" begin

#######################################################
# Calculation of Es and Ds on a single branch
# Example BiSSE calculation
# result_EsDs_biSSE_1branch_pureBirth_bl1
# (1 branch, pure birth, no Q transitions, branchlength=1)
#
# Run with:
# source("/GitHub/BioGeoJulia.jl/test/BiSSE_branchlikes_w_BD_v4_WORKING_n1.R")
# Truth:
R_result_branch_lnL = -3.128581
R_result_total_lnL = -4.937608
#######################################################


include("/GitHub/BioGeoJulia.jl/src/TreePass.jl")
import .TreePass

# Repeat calculation in Julia
include("/GitHub/BioGeoJulia.jl/notes/ModelLikes.jl")
import .ModelLikes
tr = readTopology("((chimp:1,human:1):1,gorilla:2);")
in_params = (birthRate=0.222222222, deathRate=0.1, d_val=0.0, e_val=0.0, a_val=0.1, j_val=0.0)
numstates = 2
n = 2
inputs = ModelLikes.setup_MuSSE_biogeo(numstates, tr; root_age_mult=1.5, in_params=in_params)
inputs.res.likes_at_each_nodeIndex_branchTop
inputs.res.normlikes_at_each_nodeIndex_branchTop

res = inputs.res
trdf = inputs.trdf
p_Ds_v5 = inputs.p_Ds_v5
root_age = maximum(trdf[!, :node_age])
Es_interpolator = inputs.p_Ds_v5.sol_Es_v5;

# Parameters
prtQi(inputs)
prtCi(inputs)
inputs.p_Ds_v5.params.mu_vals
p_Ds_v5.sol_Es_v5(1.0)

# Do downpass
(total_calctime_in_sec, iteration_number) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=construct_SolverOpt(), max_iterations=10^10)

Rnames(res)
res.likes_at_each_nodeIndex_branchTop
res.normlikes_at_each_nodeIndex_branchTop
res.likes_at_each_nodeIndex_branchBot
res.normlikes_at_each_nodeIndex_branchBot

sum.(res.likes_at_each_nodeIndex_branchTop)
log.(sum.(res.likes_at_each_nodeIndex_branchTop))
sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop)))
Julia_sum_lq = sum(res.lq_at_branchBot[1:(length(res.lq_at_branchBot)-1)])

# Does the total of the branch log-likelihoods match?
@test round(R_result_branch_lnL; digits=5) == round(Julia_sum_lq; digits=5)

# Add the root probabilities
# Assuming diversitree options:
# root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
# i.e., the root state probs are just the root_Ds/sum(root_Ds)
d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]
root_stateprobs = d_root_orig/sum(d_root_orig)
rootstates_lnL = log(sum(root_stateprobs .* d_root_orig))
Julia_total_lnL = lq + rootstates_lnL

# Does the total lnL match R?
@test round(R_result_total_lnL; digits=5) == round(Julia_total_lnL; digits=5)


print("\nDifferences between Julia and R lnLs for _compare_ClaSSE_calcs_v3_compare2julia.R calculation:\n")
print("R_result_branch_lnL (lq) - Julia_sum_lq: ")
print(R_result_branch_lnL - Julia_sum_lq)
print("\n")
print("R_result_total_lnL (lq) - Julia_total_lnL: ")
print(R_result_total_lnL - Julia_total_lnL)
print("\n")

end # END @testset "runtests_BiSSE_tree_n1" begin




