module Flow
print("\n\nStarting module 'Flow'...loading dependencies...\n")
using LinearAlgebra  # for mul! (matrix multiplication)
using BenchmarkTools # for @time
using InvertedIndices # for Not
using LSODA
using DifferentialEquations
using Distributed
using Random					# for MersenneTwister()
using Dates						# for e.g. DateTime, Dates.now()
using PhyloNetworks
#using Plots						# for plot
using DataFrames          # for DataFrame()
using BioGeoJulia.TrUtils # for flat2() (similar to unlist)
using BioGeoJulia.StateSpace
using BioGeoJulia.TreePass
using BioGeoJulia.SSEs

export parameterized_ClaSSE_As_v5, calc_Gs_SSE, calc_Gs_SSE!, run_Gs


# Construct interpolation function for calculating linear dynamics A, 
# at any timepoint t

#function get_LinearDynamics_A(inputs)

# This version excludes Xi (and Xj), just like castor's get_LinearDynamics_A
# calculation of A
parameterized_ClaSSE_As_v5 = (t, A, p) -> begin

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
	# Iterate through the ancestral states
  @inbounds for i in 1:n
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]

		# Calculation of "A" (the from-to probabilities between every pair of states)
		# Pull out the Q transitions - diagonal
		# case 1: no event
		A[i,i] = A[i,i]  + -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i]) # *u[i]  
		
		# case 2: anagenetic change
		@inbounds for m in 1:length(Qi_sub_i)
			A[Qi_sub_i[m],Qj_sub_i[m]] = A[Qi_sub_i[m],Qj_sub_i[m]] + Qij_vals[Qi_sub_i[m]] #* u[Qj_sub_i[m]])
		end
		
		# case 34: change + eventual extinction
		@inbounds for m in 1:length(Ci_sub_i)
			# each cladogenesis event puts probability in 2 places
			# excluding the u[], i.e. the Ds, i.e. the Xs, just as is done in 
			# 2*speciation_rates[r]*current_E[r]
			#rate_sp_then_ex = Cijk_vals[Ci_sub_i] * ((u[Ck_sub_i] * uE[Cj_sub_i]) + (u[Cj_sub_i] * uE[Ck_sub_i]))
			rate_sp_then_ex = Cijk_vals[Ci_sub_i[m]] * (uE[Cj_sub_i[m]] + uE[Ck_sub_i[m]])
			A[Ci_sub_i[m],Cj_sub_i[m]] = A[Ci_sub_i[m],Cj_sub_i[m]] + rate_sp_then_ex
			A[Ci_sub_i[m],Ck_sub_i[m]] = A[Ci_sub_i[m],Ck_sub_i[m]] + rate_sp_then_ex
		end
		
# 		du[i] = -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i])*u[i] +  # case 1: no event
# 			(sum(Qij_vals[Qi_sub_i] .* u[Qj_sub_i])) + 	# case 2	
# 			(sum(Cijk_vals[Ci_sub_i] .*                                               # case 34: change + eventual extinction
# 				 (u[Ck_sub_i].*uE[Cj_sub_i] 
# 			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))

  end # End @inbounds for i in 1:n
 	return(A)
end


calc_Gs_SSE = (dG, G, pG, t) -> begin
	A = pG.A
	p_Ds_v5 = pG.p_Ds_v5
	A = parameterized_ClaSSE_As_v5(t, A, p_Ds_v5)
	#display(A)
	#dG = A * G
	#display(G)
#	mul!(dG, A, G)
	mul!(dG, G, A)
	display(cond(G))
	print("\n")
	
	# https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/
	# When p=2, the operator norm is the spectral norm, equal to the largest singular value of A
	# When p=1, much faster, seems to always be bigger
	display(opnorm(A,1))
	display(opnorm(A,2))
	#display(opnorm(A,Inf))
	
	#display(dG)
	return(dG)
end # End calc_Gs_SSE

# Doesn't match, risky
calc_Gs_SSE! = (dG, G, pG, t) -> begin
	A = pG.A
	p_Ds_v5 = pG.p_Ds_v5
	A = parameterized_ClaSSE_As_v5(t, A, p_Ds_v5)
	#display(A)
	#dG = A * G
	mul!(dG, A, G)
	#display(dG)
	return(dG)
end # End calc_Gs_SSE



run_Gs = (p_Ds_v5) -> begin
	n = p_Ds_v5.n
	tmpzero = repeat([0.0], n^2)
	A = reshape(tmpzero, (n,n))

	G0 = reshape(tmpzero, (n,n))
	for i in 1:n
		G0[i,i] = 1.0
	end
	G = G0
	
	for i in 1:100
		t = 0.01
		A = parameterized_ClaSSE_As_v5(t, A, p_Ds_v5)
		G = A * G
		#Base.showarray(G)
		display(G)
	end

end


# This version includes Xi, Xj in the A equation (but this can't be more efficient I don't think,
# since these will be different on every branch)
parameterized_ClaSSE_A_v5xx = (du,u,A,p,t) -> begin

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
	# Iterate through the ancestral states
  @inbounds for i in 1:n
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]

		# Calculation of "A" (the from-to probabilities between every pair of states)
		# Pull out the Q transitions - diagonal
		# case 1: no event
		A[i,i] = A[i,i] + -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i])#*u[i]  
		
		# case 2: anagenetic change
		@inbounds for m in 1:length(Qi_sub_i)
			A[Qi_sub_i[m],Qj_sub_i[m]] = A[Qi_sub_i[m],Qj_sub_i[m]] + Qij_vals[Qi_sub_i[m]]# * u[Qj_sub_i[m]]
		end
		
		# case 34: change + eventual extinction
		@inbounds for m in 1:length(Ci_sub_i)
			# each cladogenesis event puts probability in 2 places
			rate_sp_then_ex = (u[Ck_sub_i] * uE[Cj_sub_i]) + ( uE[Ck_sub_i])
			A[Ci_sub_i[m],Cj_sub_i[m]] = A[Ci_sub_i[m],Cj_sub_i[m]] + rate_sp_then_ex
			A[Ci_sub_i[m],Ck_sub_i[m]] = A[Ci_sub_i[m],Ck_sub_i[m]] + rate_sp_then_ex
		end
		
# 		du[i] = -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i])*u[i] +  # case 1: no event
# 			(sum(Qij_vals[Qi_sub_i] .* u[Qj_sub_i])) + 	# case 2	
# 			(sum(Cijk_vals[Ci_sub_i] .*                                               # case 34: change + eventual extinction
# 				 (u[Ck_sub_i].*uE[Cj_sub_i] 
# 			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))

  end # End @inbounds for i in 1:n
end


end # end module Flow

