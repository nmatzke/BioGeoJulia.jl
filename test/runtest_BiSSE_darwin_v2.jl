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

# 
# """
# # Run with:
# include("/GitHub/BioGeoJulia.jl/test/runtests_BiSSE_tree_n1.jl")
# """
# 
# @testset "Example" begin
# 	@test hello("runtests_BiSSE_tree_n1.jl") == "Hello, runtests_BiSSE_tree_n1.jl"
# #	@test domath(2.0) ≈ 7.0
# end
# 
# 
# #######################################################
# # Do a bunch of tests of the SSE calculation of 
# # Ds, Es, and likelihoods, on
# # branches, nodes, and trees,
# # under a variety of simple and more complex models
# #######################################################
# 
# @testset "runtests_BiSSE_tree_n1.jl" begin
# 
#######################################################
# Calculation of Es and Ds on a single branch
# Example BiSSE calculation
# result_EsDs_biSSE_1branch_pureBirth_bl1
# (1 branch, pure birth, no Q transitions, branchlength=1)
#
# Run with:
# source("/GitHub/BioGeoJulia.jl/test/BiSSE_branchlikes_w_BD_v4_WORKING_n1.R")
# Truth:
R_result_branch_lnL = -0.4583174
R_result_total_lnL = -0.4583174
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = 12.04537
#######################################################


include("/GitHub/BioGeoJulia.jl/src/TreePass.jl")
import .TreePass

# Repeat calculation in Julia
include("/GitHub/BioGeoJulia.jl/notes/ModelLikes.jl")
import .ModelLikes
nexustr = readNexusTrees("/GitHub/BioGeoJulia.jl/wallis/Geospiza.nex")
tr = single_element_array_to_scalar(nexustr)

"""
Original tr style is: 
HybridNetwork, Rooted Network
4 edges
5 nodes: 3 tips, 0 hybrid nodes, 2 internal tree nodes.
tip labels: chimp, human, gorilla
((chimp:1.0,human:1.0):1.0,gorilla:2.0);

the nexus tree turns it into:
1-element Array{HybridNetwork,1}:
 HybridNetwork, Rooted Network
26 edges
27 nodes: 14 tips, 0 hybrid nodes, 13 internal tree nodes.
tip labels: 1, 2, 3, 4, ...
((((((((((1:0.055,2:0.055):0.055,3:0.11):0.073,4:0.183):0.009,5:0.192):0.036,6:0.228):0.103,(7:0.087,((8:0.02,9:0.02):0.015,10:0.035):0.052):0.245):0.134,11:0.466):0.069,12:0.534):0.049,13:0.583):0.297,14:0.881);

We need this to not be an array...
"""

in_params = (birthRate=3.682184, deathRate=2.263549, d_val=0.0, e_val=0.0, a_val=0.1, j_val=0.0)
numstates = 2
n = 2
inputs = ModelLikes.setup_MuSSE_biogeo(numstates, tr; root_age_mult=1.5, in_params=in_params)
(res, trdf, solver_options, p_Ds_v5, Es_tspan) = inputs
inputs.res.likes_at_each_nodeIndex_branchTop
inputs.res.normlikes_at_each_nodeIndex_branchTop

res = inputs.res
trdf = inputs.trdf
p_Ds_v5 = inputs.p_Ds_v5
root_age = maximum(trdf[!, :node_age])

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5)
# This solution is a linear interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5)

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

# Does the total of the branch log-likelihoods (lq) match?
@test round(R_result_branch_lnL; digits=5) == round(Julia_sum_lq; digits=5)

# Add the root probabilities
# Assuming diversitree options:
# root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
# i.e., the root state probs are just the root_Ds/sum(root_Ds)
d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]
root_stateprobs = d_root_orig/sum(d_root_orig)
rootstates_lnL = log(sum(root_stateprobs .* d_root_orig))
Julia_total_lnL = Julia_sum_lq + rootstates_lnL

# Does the total lnL match R?
@test round(R_result_total_lnL; digits=5) == round(Julia_total_lnL; digits=5)

# Does the total of branch likelihoods (lq) + node likelihoods match R?
# 
# R: R_result_sum_log_computed_likelihoods_at_each_node_x_lambda 
#    = sum(log(computed_likelihoods_at_each_node_x_lambda))
Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
R_sum_lq_nodes = R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
@test round(Julia_sum_lq_nodes; digits=5) == round(R_sum_lq_nodes; digits=5)


# The standard diversitree lnL calculation sums:
# 1. log-likelihoods at branch-bottoms (lq)
# 2. root-state likelihoods (e.g. rootstates_lnL)

# We can get this in R also:
# 1. sum(lq)
# 2. res1t = bisse_2areas(pars=bisse_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
#    (is sum(lq) + sum(rootstates_lnL)

# Weirdly, it seems like this R calculation extracts the total likelihood at the 
# branch bottoms (lq), but not the total likelihood at the nodes, after 
# the normlikes from above branches are combined at the node, and multiplied by
# the birthRate.
#
# This suggests another possible likelihood:
#
# 1. Julia: sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
#
# ...which matches...
#
# 2. sum(log(computed_likelihoods_at_each_node_x_lambda))
#
# (And which might be correct!  BioGeoJulia will produce both.)
#




print("\nDifferences between Julia and R lnLs for\n/GitHub/BioGeoJulia.jl/test/BiSSE_branchlikes_w_BD_v4_WORKING_n1.R\ncalculation:\n")
print("R_result_branch_lnL (lq) - Julia_sum_lq: ")
print(R_result_branch_lnL - Julia_sum_lq)
print("\n")
print("R_result_total_lnL (lq) - Julia_total_lnL: ")
print(R_result_total_lnL - Julia_total_lnL)
print("\n")
print("R_result_total_lnL (lq) - Julia_sum_lq_nodes: ")
print(R_sum_lq_nodes - Julia_sum_lq_nodes)
print("\n")

# end # END @testset "runtests_BiSSE_tree_n1" begin




