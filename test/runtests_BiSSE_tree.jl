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
include("/GitHub/BioGeoJulia.jl/test/runtests_SSE.jl")
"""

@testset "runtests_SSE.jl" begin
	@test hello("runtests_SSE.jl") == "Hello, runtests_SSE.jl"
#	@test domath(2.0) ≈ 7.0
end


#######################################################
# Do a bunch of tests of the SSE calculation of 
# Ds, Es, and likelihoods, on
# branches, nodes, and trees,
# under a variety of simple and more complex models
#######################################################

@testset "biSSE_tree_n1" begin

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
lq = sum(res.lq_at_branchBot[1:(length(res.lq_at_branchBot)-1)])

# Does the total of the branch log-likelihoods match?
@test round(R_result_branch_lnL; digits=5) == round(lq; digits=5)

# Add the root probabilities
# Assuming diversitree options:
# root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
# i.e., the root state probs are just the root_Ds/sum(root_Ds)
d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]
root_stateprobs = d_root_orig/sum(d_root_orig)
rootstates_lnL = log(sum(root_stateprobs .* d_root_orig))
total_lnL = lq + rootstates_lnL

# Does the total lnL match R?
@test round(R_result_total_lnL; digits=5) == round(total_lnL; digits=5)


#######################################################
# Check Es at t=1.0
#######################################################
Julia_result_Es = Es_interpolator(1.0)
R_Es = R_result_EsDs[(2:(n+1))]
@test all(round.(Julia_result_Es; digits=4) .== round.(R_Es; digits=4))


#######################################################
# Check Ds at t=1.0
#######################################################
n = inputs.p_Ds_v5.n
p_Ds_v5 = inputs.p_Ds_v5
u0 = collect(repeat([0.0], n))
u0[2] = 1.0
tspan = (0.0, 1.2*root_age)
prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p_Ds_v5)

ground_truth_Ds_interpolatorT = solve(prob_Ds_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_Ds_interpolatorG = solve(prob_Ds_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_Ds_interpolatorL = solve(prob_Ds_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)

ground_truth_Ds_interpolatorT.u[length(ground_truth_Ds_interpolatorT.u)]
ground_truth_Ds_interpolatorG.u[length(ground_truth_Ds_interpolatorG.u)]
ground_truth_Ds_interpolatorL.u[length(ground_truth_Ds_interpolatorL.u)]

# Check if all of the results are equal, using standard Ds SSE calcs
Tsit5_results = round.(ground_truth_Ds_interpolatorT.u[length(ground_truth_Ds_interpolatorT.u)]; digits=4)
GMRES_results = round.(ground_truth_Ds_interpolatorG.u[length(ground_truth_Ds_interpolatorG.u)]; digits=4)
LSODA_results = round.(ground_truth_Ds_interpolatorL.u[length(ground_truth_Ds_interpolatorL.u)]; digits=4)

@test all(Tsit5_results .== GMRES_results)
@test all(Tsit5_results .== LSODA_results)


#######################################################
# Check if all of the results are equal, using "Flow" Ds SSE calcs
#######################################################
include("/GitHub/BioGeoJulia.jl/notes/Flow.jl")
import .Flow

# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2)
A = reshape(tmpzero, (n,n))

# Map the likelihood "flow" of Ds, G (or Gmap or Psi).
# Start with an identity matrix
# The "I" requires "include NumericAlgebra"
G0 = Matrix{Float64}(I, n, n) 

pG = (n=n, p_Ds_v5=p_Ds_v5, A=A)
tspan = (0.0, 1.1*root_age)
prob_Gs_v5 = DifferentialEquations.ODEProblem(Flow.calc_Gs_SSE!, G0, tspan, pG)

Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
Gflow_to_01_Tsit5  = solve(prob_Gs_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
Gflow_to_01_Lsoda  = solve(prob_Gs_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)

# Check that the different interpolators match
Gmat_GMRES = round.(Gflow_to_01_GMRES(1.0); digits=5)
Gmat_Tsit5 = round.(Gflow_to_01_Tsit5(1.0); digits=5)
@test all(Gmat_GMRES .== Gmat_Tsit5)

# Lsoda requires a less precise match (!)
Gmat_GMRES = round.(Gflow_to_01_GMRES(1.0); digits=3)
Gmat_Tsit5 = round.(Gflow_to_01_Tsit5(1.0); digits=3)
Gmat_Lsoda = round.(Gflow_to_01_Lsoda(1.0); digits=3)
@test all(Gmat_GMRES .== Gmat_Lsoda)
@test all(Gmat_Tsit5 .== Gmat_Lsoda)


# Calculate the flow, on a single branch, starting from u0
# (tip values, so no fakeX0 calculation needed)
X0 = u0
Xc_GMRES = Gflow_to_01_GMRES(1.0) * X0
Xc_Tsit5 = Gflow_to_01_Tsit5(1.0) * X0
Xc_Lsoda = Gflow_to_01_Lsoda(1.0) * X0

# Compare standard to Flow
@test all(round.(ground_truth_Ds_interpolatorG(1.0); digits=6) .== round.(Xc_GMRES; digits=6))
@test all(round.(ground_truth_Ds_interpolatorT(1.0); digits=6) .== round.(Xc_Tsit5; digits=6))
@test all(round.(ground_truth_Ds_interpolatorL(1.0); digits=3) .== round.(Xc_Lsoda; digits=3)) # Much worse match

ground_truth_Ds_interpolatorG(1.0)
ground_truth_Ds_interpolatorT(1.0)
ground_truth_Ds_interpolatorL(1.0)

Ds_indices = 1 .+ collect((n+1):(2*n));
R_Ds = R_result_EsDs[Ds_indices]

# Test the R Ds, against the Julia Ds
# Standard calc of Ds
@test all(round.(ground_truth_Ds_interpolatorG(1.0); digits=4) .== round.(R_Ds; digits=4))
@test all(round.(ground_truth_Ds_interpolatorT(1.0); digits=4) .== round.(R_Ds; digits=4))
@test all(round.(ground_truth_Ds_interpolatorL(1.0); digits=3) .== round.(R_Ds; digits=3)) # Much worse match
# Flow calc of Ds
@test all(round.(R_Ds; digits=4) .== round.(Xc_GMRES; digits=4))
@test all(round.(R_Ds; digits=4) .== round.(Xc_Tsit5; digits=4))
@test all(round.(R_Ds; digits=3) .== round.(Xc_Lsoda; digits=3)) # Much worse match

R_EDs = vcat(R_Es, R_Ds)
Julia_EDs_GMRES = vcat(Julia_result_Es, Xc_GMRES)
Julia_EDs_Tsit5 = vcat(Julia_result_Es, Xc_Tsit5)
Julia_EDs_Lsoda = vcat(Julia_result_Es, Xc_Lsoda)

GMRES_Rlsoda_diffs = R_EDs .- Julia_EDs_GMRES
Tsit5_Rlsoda_diffs = R_EDs .- Julia_EDs_Tsit5
Lsoda_Rlsoda_diffs = R_EDs .- Julia_EDs_Lsoda

print("\nDifferences between Julia and R biSSE_1branch_n1 calculation:\n")
print("GMRES: ")
print(GMRES_Rlsoda_diffs)
print("\n")
print("Tsit5: ")
print(Tsit5_Rlsoda_diffs)
print("\n")
print("LSODA: ")
print(Lsoda_Rlsoda_diffs)
print("\n")

end # END @testset "biSSE_1branch_n1" begin




