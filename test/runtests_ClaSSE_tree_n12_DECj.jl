using Test, BioGeoJulia, DataFrames

using Dates									# for e.g. Dates.now(), DateTime
using PhyloNetworks					# most maintained, emphasize; for HybridNetwork
using Distributed						# for e.g. @spawn
using Combinatorics					# for e.g. combinations()
using DataFrames
using Optim                 # for e.g. L-BFGS-B Maximum Likelihood optimization

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
# include("/GitHub/BioGeoJulia.jl/test/runtests_ClaSSE_tree_n9_DECj.jl")
# """
# 
# @testset "Example" begin
# 	@test hello("runtests_ClaSSE_tree_n9_DECj.jl") == "runtests_ClaSSE_tree_n9_DECj.jl"
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
# @testset "runtests_ClaSSE_tree_n9_DECj.jl" begin
# 
#######################################################
# DEMONSTRATES MATCHING BETWEEN DIVERSITREE, BIOGEOBEARS, AND JULIA
# ON HAWAIIAN PSYCHOTRIA, 16-STATE DEC MODEL
#
# Run with:
# source("/GitHub/BioGeoJulia.jl/Rsrc/compare_BGB_diversitree_DEC+J_v1.R")
# Truth:
DEC_lnL = -34.54313
DEC_R_result_branch_lnL = -67.6295
DEC_R_result_total_LnLs1 = -72.60212
DEC_R_result_total_LnLs1t = -71.48986
DEC_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -120.1545

# DEC+J
DECj_lnL = -20.94759
R_result_branch_lnL = -55.37332
R_result_total_LnLs1 = -58.83758
R_result_total_LnLs1t = -57.72533
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -96.34151

#######################################################


include("/GitHub/BioGeoJulia.jl/src/TreePass.jl")
import .TreePass

# Repeat calculation in Julia
include("/GitHub/BioGeoJulia.jl/notes/ModelLikes.jl")
import .ModelLikes

# Read in the island areas using Parsers.jl (v10) rather than just manual input (v9)
include("/GitHub/BioGeoJulia.jl/notes/Parsers.jl")
import .Parsers

# Island numbers (KOMH = 1234) in Rnodenums order:
#island_nums = [3, 3, 2, 2, 3, 3, 2, 1, 1, 3, 4, 2, 1, 1, 1, 1, 1, 1, 2]
lgdata_fn = "/GitHub/BioGeoJulia.jl/Rsrc/Psychotria_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

# Psychotria tree
tr = readTopology("((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);")
in_params = (birthRate=0.3288164, deathRate=0.0, d_val=1e-12, e_val=1e-12, a_val=0.0, j_val=0.1142057)
numareas = 4
n = 16            # 4 areas, 16 states

# CHANGE PARAMETERS BEFORE E INTERPOLATOR
inputs = ModelLikes.setup_DEC_SSE(numareas, tr; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, in_params=in_params)
(setup, res, trdf, solver_options, p_Ds_v5, Es_tspan) = inputs

# Update the tip likelihoods, with the geography data
inputs.res.likes_at_each_nodeIndex_branchTop
size(inputs.res.likes_at_each_nodeIndex_branchTop)
numstates = length(inputs.res.likes_at_each_nodeIndex_branchTop[1])

inputs.setup.observed_statenums
inputs.res.likes_at_each_nodeIndex_branchTop

inputs = Parsers.tipranges_to_tiplikes(inputs, geog_df);
inputs.setup.observed_statenums
inputs.res.likes_at_each_nodeIndex_branchTop


inputs.res.likes_at_each_nodeIndex_branchTop
res = inputs.res

trdf = inputs.trdf
p_Ds_v5 = inputs.p_Ds_v5
root_age = maximum(trdf[!, :node_age])

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5)
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5)

prtQi(inputs)
prtCi(inputs)
inputs.p_Ds_v5.params.mu_vals
p_Ds_v5.sol_Es_v5(1.0)

Es_interpolator(1.0)


# Parameters

# Do downpass
res.likes_at_each_nodeIndex_branchTop
(total_calctime_in_sec, iteration_number) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=construct_SolverOpt(), max_iterations=10^10)
res.likes_at_each_nodeIndex_branchTop

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
@test round(R_result_branch_lnL; digits=4) == round(Julia_sum_lq; digits=4)

# Add the root probabilities

# Assuming diversitree options:
# root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
# i.e., the root state probs are just the root_Ds/sum(root_Ds)
d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]
root_stateprobs = d_root_orig/sum(d_root_orig)
rootstates_lnL = log(sum(root_stateprobs .* d_root_orig))
Julia_total_lnLs1 = Julia_sum_lq + rootstates_lnL


# Assuming diversitree options:
# root=ROOT.OBS, root.p=NULL, condition.surv=TRUE (these seem to be the defaults)
# i.e., the root state probs are just the root_Ds/sum(root_Ds)
d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]
root_stateprobs = d_root_orig/sum(d_root_orig)
lambda = in_params.birthRate
e_root = Es_interpolator(root_age)


#d_root = d_root_orig ./ sum(root_stateprobs .* inputs.p_Ds_v5.params.Cijk_vals .* (1 .- e_root).^2)
# diversitree::rootfunc.classe
# lambda <- colSums(matrix(pars[1:(nsum * k)], nrow = nsum))
sum_of_lambdas = collect(repeat([0.0], n))
for i in 1:n
	sum_of_lambdas[i] = sum(inputs.p_Ds_v5.params.Cijk_vals[inputs.p_Ds_v5.p_TFs.Ci_eq_i[i]])
end
sum_of_lambdas
d_root = d_root_orig ./ sum(root_stateprobs .* sum_of_lambdas .* (1 .- e_root).^2)
rootstates_lnL = log(sum(root_stateprobs .* d_root))
# The above all works out to [0,1] for the Yule model with q01=q02=0.0

Julia_total_lnLs1t = Julia_sum_lq + rootstates_lnL

# Does the total lnL match R?
# root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
@test round(R_result_total_LnLs1; digits=4) == round(Julia_total_lnLs1; digits=4)

# root=ROOT.OBS, root.p=NULL, condition.surv=TRUE
@test round(R_result_total_LnLs1t; digits=4) == round(Julia_total_lnLs1t; digits=4)

# Does the total of branch likelihoods (lq) + node likelihoods match R?
# 
# R: R_result_sum_log_computed_likelihoods_at_each_node_x_lambda 
#    = sum(log(computed_likelihoods_at_each_node_x_lambda))
res.likes_at_each_nodeIndex_branchTop
log.(sum.(res.likes_at_each_nodeIndex_branchTop))

sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop)))

Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
R_sum_lq_nodes = R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
@test round(Julia_sum_lq_nodes; digits=4) == round(R_sum_lq_nodes; digits=4)








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




print("\nDifferences between Julia and R lnLs for\n/GitHub/BioGeoJulia.jl/Rsrc/_compare_ClaSSE_calcs_v3_compare2julia.R\n calculation:\n")
print("R_result_branch_lnL (lq) - Julia_sum_lq: ")
print(R_result_branch_lnL - Julia_sum_lq)
print("\n")
print("R_result_total_LnLs1 (lq) - Julia_total_lnLs1: ")
print(R_result_total_LnLs1 - Julia_total_lnLs1)
print("\n")
print("R_result_total_LnLs1t (lq) - Julia_total_lnLs1t: ")
print(R_result_total_LnLs1t - Julia_total_lnLs1t)
print("\n")
print("R_result_total_lnL (lq) - Julia_sum_lq_nodes: ")
print(R_sum_lq_nodes - Julia_sum_lq_nodes)
print("\n")


#######################################################
# OK, let's do an ML inference
#######################################################
# 1. Figure out how ML inference works
# 
# https://julianlsolvers.github.io/Optim.jl/stable/#examples/generated/maxlikenlm/
#
# 
# 2. Write function to update parameters

pars = [0.3, 0.2, 1.0, 1.0, 1.0, 0.0]
parnames = ["d", "e", "y", "s", "v", "j"]

# Set inputs to starting values
Rnames(p_Ds_v5)

# Yule birthRate
birthRate = 0.3288164

prtQp(p_Ds_v5)
prtCp(p_Ds_v5)
Qdf = prtQi(inputs)
Cdf = prtCi(inputs)

for i in 1:length(parnames)
	# Update the Q values
	TF1 = parnames[i] .== Qdf[:,:event]
	Qdf[:,:val][TF1] .= pars[i]
	p_Ds_v5.params.Qij_vals[TF1] .= pars[i]
	#print(sum(TF1))
	
	# Update the C values
	TF2 = parnames[i] .== Cdf[:,:event]
	Cdf[:,:val][TF2] .= pars[i]
	p_Ds_v5.params.Cijk_weights[TF2] .= pars[i]
	#print(sum(TF2))
end

# Convert the Cijk_weights to birthRates

# for each possible ancestral state
num_possible_ancestral_states = length(p_Ds_v5.p_TFs.Ci_sub_i)
sum_of_weights_per_ancstate = collect(repeat([0.0], num_possible_ancestral_states))
for i in 1:num_possible_ancestral_states
	sum_of_weights_per_ancstate[i] = sum(p_Ds_v5.params.Cijk_weights[p_Ds_v5.p_TFs.Ci_eq_i[i]])
	# Convert to rates
	p_Ds_v5.params.Cijk_vals[p_Ds_v5.p_TFs.Ci_eq_i[i]] = birthRate .* p_Ds_v5.params.Cijk_weights[p_Ds_v5.p_TFs.Ci_eq_i[i]] ./ sum_of_weights_per_ancstate[i]
end # end for-loop
sum_of_weights_per_ancstate

prtQp(p_Ds_v5)
prtCp(p_Ds_v5)



function func_to_optimize(pars, parnames, inputs, p_Ds_v5; returnval="lnL")
	# Get the Q, C
	res = inputs.res
	trdf = inputs.trdf
	
	# DON'T take p_Ds_v5 from inputs, because that nukes the sol_Es that you modified!
	#p_Ds_v5=inputs.p_Ds_v5
	
	Qdf = prtQi(inputs)
	Cdf = prtCi(inputs)
	# Rnames(inputs.p_Ds_v5.params)
	
	# Update the parameters
	# (this will only work for simple cases where no formula is needed)
	for i in 1:length(parnames)
		# Update the Q values
# 		TF1 = parnames[i] .== Qdf[:,:event]
# 		Qdf[:,:val][TF1] .= pars[i]
# 		p_Ds_v5.params.Qij_vals[TF1] .= pars[i]
		#print(sum(TF1))
		
		# Update the C values
		TF2 = parnames[i] .== Cdf[:,:event]
		Cdf[:,:val][TF2] .= pars[i]
		p_Ds_v5.params.Cijk_vals[TF2] .= pars[i]
		#print(sum(TF2))
	end
	
	# Update the Qmat
	Qmat = (Qarray_ivals=p_Ds_v5.p_indices.Qarray_ivals, Qarray_jvals=p_Ds_v5.p_indices.Qarray_jvals, Qij_vals=p_Ds_v5.params.Qij_vals, Qarray_event_types=p_Ds_v5.p_indices.Qarray_event_types)
	
	areas_list = inputs.setup.areas_list[:]
	states_list = inputs.setup.states_list[:]
	
	dmat = reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list)))
	d_val = 0.0
	TF = parnames .== "d"
	d_val = pars[TF][1]
	dmat[:] = dmat .* d_val
	amat = deepcopy(dmat)

	elist=repeat([1.0], length(areas_list))
	e_val = 0.0
	TF = parnames .== "e"
	e_val = pars[TF][1]
	elist[:] = elist .* e_val
	
	# Update
	Qmat2 = update_Qij_vals(Qmat, areas_list, states_list, dmat, elist, amat)
	p_Ds_v5.params.Qij_vals[:] = Qmat2.Qij_vals
	
	#res2 = deepcopy(res)
	(total_calctime_in_sec, iteration_number, Julia_sum_lqA, rootstates_lnLA, Julia_total_lnLs1A) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

	
	txt = paste0(["pars[1]=", pars[1], ", pars[2]=", pars[2], ",	Julia_sum_lqA=", round(Julia_sum_lqA; digits=3), ", rootstates_lnLA=", round(rootstates_lnLA; digits=3), ",	Julia_total_lnLs1A=", Julia_total_lnLs1A])
	print(txt) 
	print("\n")
	
	if returnval == "lnL"
		return(-Julia_total_lnLs1A)
	end
	if returnval == "inputs"
		return(inputs)
	end
	# Shouldn't get here
	return(NaN)
end # END function func_to_optimize(pars, parnames)



func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="lnL")

pars = [0.3, 0.2]
parnames = ["d", "e"]
lnL = func_to_optimize(pars, parnames, inputs, p_Ds_v5; returnval="lnL")
func(pars)

pars = [0.03, 0.02]
parnames = ["d", "e"]
lnL = func_to_optimize(pars, parnames, inputs, p_Ds_v5; returnval="lnL")
func(pars)


pars = [0.3, 0.2]
parnames = ["d", "e"]
lower = [0.0, 0.0]
upper = [5.0, 5.0]
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="lnL")
MLres = optimize(func, lower, upper, pars)#, Fminbox(LBFGS()))

p_Ds_v5.params.Qij_vals .= 0.1
iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

end # END @testset "runtests_BiSSE_tree_n3" begin

# Update Qmat
Qmat = (Qarray_ivals=p_Ds_v5.p_indices.Qarray_ivals, Qarray_jvals=p_Ds_v5.p_indices.Qarray_jvals, Qij_vals=p_Ds_v5.params.Qij_vals, Qarray_event_types=p_Ds_v5.p_indices.Qarray_event_types)
Qmat2 = update_Qij_vals(Qmat, areas_list, states_list, dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))), elist=repeat([1.0], length(areas_list)), amat=dmat )

