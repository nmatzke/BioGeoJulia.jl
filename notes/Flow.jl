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

# 	modelD.getLinearDynamics(ages[a], dynamics); // get linear dynamics at this age
# 		// calculate largest singular value (sigma1) of the dynamics at this age
# 		// then kappa_rate <= 2*sigma1, since sigma1(exp(t*A))/sigma2(exp(t*A)) <= exp(t*2*sigma1(A)) [So and Thompson (2000) Singular values of Matrix Exponentials. Theorem 2.1]
#
# NOTE: In Julia, 
# When p=2, the operator norm is the spectral norm, equal to the largest singular value of A.
# 
# 
#
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

#
calc_Gs_SSE = (dG, G, pG, t) -> begin
	A = pG.A
	p_Ds_v5 = pG.p_Ds_v5
	A = parameterized_ClaSSE_As_v5(t, A, p_Ds_v5)
	#display(A)
	#dG = A * G
	#display(G)
#	mul!(dG, A, G)
	mul!(dG, G, A)
	condG1 = cond(G,1)
	condG2 = cond(G,2)
	condGInf = cond(G,Inf)
	tmpstr = paste0(["\nAt time t=", string(round(t, digits=6)), ", Condition number=", string(round(condG1, digits=6)), ", ", string(round(condG2, digits=6)), ", ", string(round(condGInf, digits=6))])
	print(tmpstr)
	#display(cond(G))
	
	# From Louca & Pennell code:
	# phylogenetics_cpp_routines.cpp
	# max_condition_number,			// (INPUT) unitless number, 
	# the maximum acceptable condition number for the Gmap 
	# (as estimated from the linearized dynamics), when choosing 
	# the integration interval size. A larger max_condition number 
	# leads to fewer age-splits, thus faster computation but also 
	# lower accuracy. Hence, this number controls the trade-off 
	# between speed and accuracy. Typical values are 
	# 1e4 (slower, more accurate) up to 
	# 1e8 (faster, less accurate).
	
	# The condition number is kappa:
	# https://julia.quantecon.org/tools_and_techniques/iterative_methods_sparsity.html
	# The cond() operation can be expensive, so the inequality in Louca, Supp. Mat. Eq. 3 is useful
	# It looks *extremely* conservative, it blows up at e.g.
	# time = 0.00015
	#
	# upper_bound_condition_number: 18425.777249466213
	# 
	# At time t=0.00015, Condition number=2.175764, 1.658813, 1.70285
	# opnorms p=1 & p=2:
	# 32759.807937823098
	# 30274.76762619003
	# 30274.767626190034
	# 17479.145238634184
	#
	# vs.
	# 
	# upper_bound_condition_number: 1.0971658583069032e6
	# 
	# At time t=0.0002, Condition number=3.261165, 2.255482, 2.334215
	# opnorms p=1 & p=2:
	# 34805.07811454515
	# 32164.890247616975
	# 32164.89024761698
	# 18570.408042916435
	# 
	# Probably we could just use opnorm(A,1), or even opnorm(A,1)^2
	
	
	# Matrix norms
	# See: 
	# Lambers, Jim (2009). Vector Norms & Matrix Norms. pp. 1-16.
	# https://www.math.usm.edu/lambers/mat610/sum10/lecture2.pdf
	# These notes have the inequalities between norm forms
	
	# https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/
	# Note: operator norm = matrix norm
	# When p=1, ||A||1, much faster, seems to always be bigger
	# When p=2, ||A||2, the operator norm is the spectral norm, equal to the largest singular value of A
	# When p=Inf, ||A||Inf, this is just opnorm(t(A),1), and 
	#
	# ||A||2 <= sqrt(ncol(A))*||A||Inf
	# ||A||2 <= sqrt(nrow(A))*||A||Inf
	# 
	# (doesn't seem to be true, actually. But cond of opnorm(1) does seem conservative
	
	# Actually, instead of tracking the condition number kappa, they are tracking the 
	# *growth rate* of kappa:
	#
	# // determine linear dynamics (matrix form) of D at various representative 
	# ages, and keep track of the largest estimated exponential growth rate of 
	# the condition number of Gmap ("kappa_rate")
	# 
	# // calculate largest singular value (sigma1) of the dynamics at this age
	#	// then kappa_rate <= 2*sigma1, since 
	# // sigma1(exp(t*A))/sigma2(exp(t*A)) <= exp(t*2*sigma1(A))
	# [So and Thompson (2000) Singular values of Matrix Exponentials. Theorem 2.1]
	# 
	
	print("\nopnorms p=1 & p=2:\n")
	display(opnorm(A,1))
	display(opnorm(A,2))
	display(opnorm(transpose(A),2))
	display(1/(sqrt(size(A)[1]))*opnorm(transpose(A),2))
	print("\n")
	#display(opnorm(A,Inf))

	sigma1_of_A = opnorm(A,1)
	upper_bound_condition_number = exp(2*t*sigma1_of_A)
	upper_bound_kappa_growth_rate = 2*sigma1_of_A
	tmpstr = paste0(["\nupper_bound_condition_number=", string(upper_bound_condition_number), ", upper_bound_kappa_growth_rate=", string(upper_bound_kappa_growth_rate)])
	print(tmpstr)
	print("\n")
	


	
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

