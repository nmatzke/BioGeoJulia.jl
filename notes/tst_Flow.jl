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


	using Profile     # for @profile
	using DataFrames  # for DataFrame
	using PhyloNetworks
	using BioGeoJulia.TrUtils # for flat2() (similar to unlist)
	using BioGeoJulia.StateSpace
	using BioGeoJulia.TreePass
	using BioGeoJulia.SSEs

	using OrdinaryDiffEq, Sundials, DiffEqDevTools, Plots, ODEInterfaceDiffEq, ODE, LSODA
	#Pkg.add(PackageSpec(url="https://github.com/JuliaDiffEq/deSolveDiffEq.jl"))
	#using deSolveDiffEq 
	# https://docs.juliadiffeq.org/stable/solvers/ode_solve/index.html

	using Profile     # for @profile
	using DataFrames  # for DataFrame
	using PhyloNetworks
	using RCall       # for df_to_Rdata, reval, g = globalEnv


	inputs = ModelLikes.setup_DEC_SSE(2, readTopology("((chimp:1,human:1):1,gorilla:2);"))
# 	for i in 1:length(inputs.p_Ds_v5.params.Qij_vals)
# 		inputs.p_Ds_v5.params.Qij_vals[i] = inputs.p_Ds_v5.params.Qij_vals[i] / 100
# 	end
	
	
	res = inputs.res
	trdf = inputs.trdf
	solver_options = inputs.solver_options
	p_Ds_v5 = inputs.p_Ds_v5
	
	
	p_Ds_v5.sol_Es_v5
	
	mu_vals = p_Ds_v5.params.mu_vals
	Qarray_ivals = p_Ds_v5.p_indices.Qarray_ivals
	Qarray_jvals = p_Ds_v5.p_indices.Qarray_jvals
	Carray_ivals = p_Ds_v5.p_indices.Carray_ivals
	Carray_jvals = p_Ds_v5.p_indices.Carray_jvals
	Carray_kvals = p_Ds_v5.p_indices.Carray_kvals
  Qij_vals = p_Ds_v5.params.Qij_vals
  Cijk_vals = p_Ds_v5.params.Cijk_vals
	
	Rcbind(Qarray_ivals, Qarray_jvals, Qij_vals)
	Rcbind(Carray_ivals, Carray_jvals, Carray_kvals, Cijk_vals)
	
	u0 = collect(repeat([0.0], n))
	u0[2] = 1.0
	tspan = (0, 2.5)
	prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p_Ds_v5)

	
	
	
	# build an A
	n = p_Ds_v5.n
	tmpzero = repeat([0.0], n^2)
	A = reshape(tmpzero, (n,n))
	#parameterized_ClaSSE_As_v5 = (dA,p,t)
	
	sol_Es = p_Ds_v5.sol_Es_v5
	
	# A(t) is just a function of t, given the parameters and the Es
	
	tvals = [0.01, 0.1, 1.0, 10.0, 100.0]
	sol_Es.(tvals)
	
	for t in tvals
		show(Flow.parameterized_ClaSSE_As_v5(t, A, p_Ds_v5))
	end
	
	
	include("Flow.jl")
	import .Flow
	t = 0.01
	Flow.parameterized_ClaSSE_As_v5(t, A, p_Ds_v5)
	t = 0.1
	Flow.parameterized_ClaSSE_As_v5(t, A, p_Ds_v5)
	t = 1.0
	Flow.parameterized_ClaSSE_As_v5(t, A, p_Ds_v5)
	
	pG = (p_Ds_v5=p_Ds_v5, A=A)

	tmpzero = repeat([0.0], n^2)
	G0 = reshape(tmpzero, (n,n))
	for i in 1:n
		G0[i,i] = 1.0
	end

	prob_Gs_v5 = DifferentialEquations.ODEProblem(Flow.calc_Gs_SSE, G0, tspan, pG)
	sol = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
	sol = solve(prob_Gs_v5, BS5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
	sol = solve(prob_Gs_v5, Rosenbrock23(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
	sol = solve(prob_Gs_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
	

	
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

