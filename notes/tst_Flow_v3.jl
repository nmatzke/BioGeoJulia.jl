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

	using LinearAlgebra  # for "I" in: Matrix{Float64}(I, 2, 2)
											 # https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
	using Profile     # for @profile
	using DataFrames  # for DataFrame
	using PhyloNetworks
	using BioGeoJulia.TrUtils # for flat2() (similar to unlist)
	using BioGeoJulia.StateSpace
	using BioGeoJulia.TreePass
	using BioGeoJulia.SSEs
	
	using DifferentialEquations
	using OrdinaryDiffEq, Sundials, DiffEqDevTools, Plots, ODEInterfaceDiffEq, ODE, LSODA
	#Pkg.add(PackageSpec(url="https://github.com/JuliaDiffEq/deSolveDiffEq.jl"))
	#using deSolveDiffEq 
	# https://docs.juliadiffeq.org/stable/solvers/ode_solve/index.html

	using Profile     # for @profile
	using DataFrames  # for DataFrame
	using PhyloNetworks
#	using RCall       # for df_to_Rdata, reval, g = globalEnv

	# Set up a DEC-like model; will calculate Es over 120% of root depth
	#inputs = ModelLikes.setup_DEC_SSE(2, readTopology("((chimp:10,human:10):10,gorilla:20);"))
	inputs = ModelLikes.setup_MuSSE(2, readTopology("((chimp:10,human:10):10,gorilla:20);"))
	res = inputs.res
	trdf = inputs.trdf
	root_age = maximum(trdf[!, :node_age])
	solver_options = inputs.solver_options
	solver_options.save_everystep
	p_Ds_v5 = inputs.p_Ds_v5
	Rnames(p_Ds_v5)
	
	# Check that the E interpolator flatted out at large times
	sol_Es = p_Ds_v5.sol_Es_v5
	sol_Es(seq(1.0, root_age, 1.0) )
	
	# Look at the model parameters (Q and C matrix)
	Rcbind(p_Ds_v5.p_indices.Qarray_ivals, p_Ds_v5.p_indices.Qarray_jvals, p_Ds_v5.params.Qij_vals)
	Rcbind(p_Ds_v5.p_indices.Carray_ivals, p_Ds_v5.p_indices.Carray_jvals, p_Ds_v5.p_indices.Carray_kvals, p_Ds_v5.params.Cijk_vals)
	

	# Run numerical integration of the Ds, for ground truth
	# Creates an interpolator, ground_truth_Ds_interpolator, to produce
	# Ds at any timepoint
	n = inputs.p_Ds_v5.n
	u0 = collect(repeat([0.0], n))
	u0[2] = 1.0
	tspan = (0.0, 1.2*root_age)
	prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p_Ds_v5)
	ground_truth_Ds_interpolator = solve(prob_Ds_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
	ground_truth_Ds_interpolator = solve(prob_Ds_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
	ground_truth_Ds_interpolator[length(ground_truth_Ds_interpolator)]
	
	# Now, do this calculation with the "Flow" calculation
	include("Flow.jl")
	import .Flow
	
	# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
	tmpzero = repeat([0.0], n^2)
	A = reshape(tmpzero, (n,n))
	
	
	tvals = [0.01, 0.1, 1.0, 10.0]
	sol_Es.(tvals)
	
	# A(t) is just a function of t, given the parameters and the Es
	# Using  A[:,:] is required, so that the starting "A" doesn't mutate
	for t in tvals
		A_at_t = Flow.parameterized_ClaSSE_As_v5(t, A[:,:], p_Ds_v5)
		display(A_at_t)
#		print("\n")
	end
	
	# Using  A[:,:] is required, so that the starting "A" doesn't mutate
	t = 0.01
	Flow.parameterized_ClaSSE_As_v5(t, A[:,:], p_Ds_v5)
	t = 0.1
	Flow.parameterized_ClaSSE_As_v5(t, A[:,:], p_Ds_v5)
	t = 1.0
	Flow.parameterized_ClaSSE_As_v5(t, A[:,:], p_Ds_v5)


	# Check the linear dynamics (matrix A) for timespans where the kappa rate
	# exceeds log(max_condition_number). max_condition_number is usually between
	# 1e4 (slower) and 1e8 (faster)
	# 
	# "exponential growth rate of the condition number ("kappa_rate") based on the 
	#  largest singular value of the dynamics A" 
	tvals = collect(0:0.1:root_age)
	kappa_Arates_df = Flow.check_linearDynamics_of_As(tvals, p_Ds_v5; max_condition_number=1e8)
	
	# Are there any condition numbers too big?
	seq(1,nrow(kappa_Arates_df),1)[kappa_Arates_df[!,:condbigTF]]

	# Map the likelihood "flow" of Ds, G (or Gmap or Psi).
	# Start with an identity matrix
	# The "I" requires "include NumericAlgebra"
	G0 = Matrix{Float64}(I, n, n) 

	pG = (n=n, p_Ds_v5=p_Ds_v5, A=A)
	tspan = (0.0, 20.0)
	prob_Gs_v5 = DifferentialEquations.ODEProblem(Flow.calc_Gs_SSE!, G0, tspan, pG)

	Gflow_to_01g  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
	Gflow_to_01t  = solve(prob_Gs_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
	Gflow_to_01l  = solve(prob_Gs_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
	
	# Check that the different interpolators match
	Gflow_to_01g(0.001)
	Gflow_to_01t(0.001)
	Gflow_to_01l(0.001)

	Gflow_to_01g(0.01)
	Gflow_to_01t(0.01)
	Gflow_to_01l(0.01)

	Gflow_to_01g(0.1)
	Gflow_to_01t(0.1)
	Gflow_to_01l(0.1)

	Gflow_to_01g(1.0)
	Gflow_to_01t(1.0)
	Gflow_to_01l(1.0)

	Gflow_to_01g(10.0)
	Gflow_to_01t(10.0)
	Gflow_to_01l(10.0)

	
	mean(Gflow_to_01.u[1])
	mean(Gflow_to_01.u[2])

	
	# OK, we now have an equation that calculates the flow, G, down any timespan
	
	#######################################################
	# TEST the flow calculation, on a single branch
	#######################################################
	
	# Calculate Ds at t=10
	factored_G = factorize(Gflow_to_01l(10.0))
	# Solve for imaginary X0 at t=0 that would produce
	Xc = ground_truth_Ds_interpolator(10.0)
	X0 = factored_G \ Xc
	
	X0_v2 = similar(Xc)
	ldiv!(X0_v2, factored_G, Xc)
	X0
	X0_v2
	
	
	# Matrix x vector product, stored in Xc_v2
	Xc_v2 = Gflow_to_01l(10.0) * X0
	Xc_v2
	Xc
	
	
	X0_v3 = similar(Xc)
	mul!(X0_v3, Gflow_to_01l(10.0), X0)
	
	
	
	
	#######################################################
	# KEY IDEA
	
	# To get from the likelihoods at:
	# Xchild  at time tc
	# ...to...
	# Xparent at time tp
	# ...normally you would integrate the SSE equation down 
	# the branch, with the starting value of Xchild.
	
	# HOWEVER, if you have the flow (G or Gmap or psi(tc, tp))
	# then you can just go 
	# psi(t0, tp) / psi(t0,tc) * Xc
	# which is equivalent to 
	# psi(t0, tp) * fakeX0
	# Gmap(tp) * fakeX0
	#
	# ...because Gmap(tc) * fakeX0 = Xc
	# Solve for fakeX0 with G(tc) \ Xc = fakeX0
	#######################################################
	
	# Let's calculate Xc, starting from a tip
	X0 = [0.0, 1.0, 0.0]
	
	tc = 1.0
	Xc_from_flow = similar(X0)
	mul!(Xc_from_flow, Gflow_to_01l(tc), X0)
	Xc_from_flow2 = Gflow_to_01g(tc) * X0
	ground_truth_Ds_interpolator(tc)

	tc = 2.0
	Xc_from_flow = similar(X0)
	mul!(Xc_from_flow, Gflow_to_01l(tc), X0)
	Xc_from_flow2 = Gflow_to_01g(tc) * X0
	ground_truth_Ds_interpolator(tc)

	tc = 3.0
	Xc_from_flow = similar(X0)
	mul!(Xc_from_flow, Gflow_to_01l(tc), X0)
	Xc_from_flow2 = Gflow_to_01g(tc) * X0
	ground_truth_Ds_interpolator(tc)
	Xc_from_flow2 ./ ground_truth_Ds_interpolator(tc) 


	tc = 10.0
	Xc_from_flow = similar(X0)
	mul!(Xc_from_flow, Gflow_to_01l(tc), X0)
	Xc_from_flow2 = Gflow_to_01g(tc) * X0
	ground_truth_Ds_interpolator(tc)
	Xc_from_flow2 ./ ground_truth_Ds_interpolator(tc) 
	log.(Xc_from_flow2 ./ ground_truth_Ds_interpolator(tc))



	factored_G = factorize(Gflow_to_01l(0.0))
	fakeX0 = factored_G \ X0
	Xc_from_flow3 = Gflow_to_011(10.0) * fakeX0

	
	factored_G = factorize(Gflow_to_01l(10.0))
	fakeX0 = factored_G \ X0
	Xc_from_flow3 = Gflow_to_01g(10.0) * fakeX0
	
	
	#@benchmark sol = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=false, abstol = 1e-9, reltol = 1e-9)
	u1 = Gflow.u
	normalized_u1 = u1 ./ (sum.(u1))
	display(normalized_u1[1])
	display(normalized_u1[2])
	
	for i in 1:length(u1)
		print(sum(all.(u1[i] .>= 0.0)))
		print("\n")
	end
	all.(u1[1] .>= 0.0)

	# Doesn't seem to work?
	prob_Gs_v5 = DifferentialEquations.ODEProblem(Flow.calc_Gs_SSE!, G0, (0.0, 0.1), pG)
	sol2 = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
	u2 = sol2.u
	#@benchmark sol = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=false, abstol = 1e-9, reltol = 1e-9)
	
	Gflow(0.01)
	sol2(0.01)
	
	
	u1[length(u1)]
	u2[length(u2)]
	

#	@benchmark sol = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
#	@benchmark sol = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=false, abstol = 1e-9, reltol = 1e-9)
	
# 	u1 ./ (sum.(u1))
# 	u2 ./ (sum.(u2))
	
	
# 	# SLOWWW
# 	sol = solve(prob_Gs_v5, BS5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
# 	u2 = sol.u
# 	
# 	@benchmark sol = solve(prob_Gs_v5, BS5(), save_everystep=false, abstol = 1e-9, reltol = 1e-9)
# 	
# 	sol = solve(prob_Gs_v5, Rosenbrock23(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
# 	sol = solve(prob_Gs_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
# 	u3 = sol.u

	hcat(sum.(u1), sum.(u2))
	
	tmpzero = repeat([0.0], n^2)
	A = reshape(tmpzero, (n,n))

	G0 = reshape(tmpzero, (n,n))
	for i in 1:n
		G0[i,i] = 1.0
	end
	G = G0
	
	for i in 1:100
		t = 0.01
		A = Flow.parameterized_ClaSSE_As_v5(t, A, p_Ds_v5)
		G = A * G
		print(G)
	end
	
	include("Flow.jl")
	import .Flow
	Flow.run_Gs(p_Ds_v5)
	
	
	Es = inputs.p_Ds_v5.sol_Es_v5
	
	Rnames(inputs.p_Ds_v5.p_indices)
	
	
  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  Qij_vals = p.params.Qij_vals
  Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Qarray_ivals = p.p_indices.Qarray_ivals
	Qarray_jvals = p.p_indices.Qarray_jvals
	Carray_ivals = p.p_indices.Carray_ivals
	Carray_jvals = p.p_indices.Carray_jvals
	Carray_kvals = p.p_indices.Carray_kvals
	
	# Pre-calculated solution of the Es
	sol_Es = p.sol_Es_v5
	uE = p.uE
	uE = sol_Es(t)
	
	two = 1.0
  @inbounds for i in 1:n
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]

		# Calculation of "D" (likelihood of tip data)
		du[i] = -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i])*u[i] +  # case 1: no event
			(sum(Qij_vals[Qi_sub_i] .* u[Qj_sub_i])) + 	# case 2	
			(sum(Cijk_vals[Ci_sub_i] .*                                               # case 34: change + eventual extinction
				 (u[Ck_sub_i].*uE[Cj_sub_i] 
			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))
  end
  	
	# Number of lambdas and Qs
	numLambdas = length(Cijk_vals)
	numQs = length(Qij_vals)
	
	# build a Qmat-sized A matrix to hold the linear dynamics at a particular timepoint
	numstates = inputs.p_Ds_v5.n
	tmpzero = repeat([0.0], numstates^2)
	A = reshape(tmpzero, (numstates,numstates))

	# For Es at a given t, the rate of change A[i,j] is:
	
	# The Q matrix
	@inbounds for m in 1:numQs
		A[Qarray_ivals[m], Qarray_jvals[m]] = A[Qarray_ivals[m], Qarray_jvals[m]] + Qij_vals[m]
	end

	# Do the diagonal of A
	# This is the diagonal of the Q transition matrix
	# This accounts for the -sum(Q {i != j}) in the "no events occured" category
	@inbounds for m in 1:numstates
		A[m,m] = 0.0
		A[m,m] = -sum(A[m,:]) # sum the rows, ie standard Q matrix
	end
	
	# Add the lamda no-change by ancestor i
  @inbounds for i in 1:n
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]

		A[i,i] = A[i,i] - (sum(Cijk_vals[Ci_sub_i]) + mu[i]) # case 1: no event occured (Q-diag is above)
	end
	
	# Instead of assuming all (speciation + extinction of 1) events are range-copying events,
	# we have to iterate through the ancestral states i and add to the Q (actually A) matrix
	
	# Add the lambda transition events (speciational)
	# (this is not found in castor's getLinearDynamics, as MuSSE doesn't have cladogenesis events)
	@inbounds for m in 1:numLambdas
		A[Carray_ivals[m], Carray_jvals[m]] = A[Carray_ivals[m], Carray_jvals[m]] + Cijk_vals[m]
		A[Carray_ivals[m], Carray_kvals[m]] = A[Carray_ivals[m], Carray_kvals[m]] + Cijk_vals[m]
	end
	
	# Do the diagonal of A
	@inbounds for m in 1:numstates
	A[m,m] = 0.0
	
	
	end	
	
	
end

