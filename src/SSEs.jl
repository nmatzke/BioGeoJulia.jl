module SSEs
using DataFrames  # for e.g. DataFram()
export parameterized_ClaSSE, parameterized_ClaSSE_Es, parameterized_ClaSSE_Ds, parameterized_ClaSSE_v5, parameterized_ClaSSE_Es_v5, parameterized_ClaSSE_Ds_v5










#######################################################
# Various versions of ClaSSE equations
# trying for more and more efficient
#######################################################


# Load g3 -- 100 states
#include("/drives/GDrive/__GDrive_projects/2018-01-22_Marsden/software/Julia_EPIRK/ClaSSE_25states_v1.jl")
#include("/drives/GDrive/__GDrive_projects/2018-01-22_Marsden/software/Julia_EPIRK/ClaSSE_25states_v1.jl")

# Seems to be the same as
# function parameterized_ClaSSE!(du,u,p,t)
# (in-place function with "!")
parameterized_ClaSSE = (du,u,p,t) -> begin

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
	
	two = 1.0
  @inbounds for i in 1:n
		# Calculation of "E" (prob. of extinction)
		du[i] = mu[i] +                                         # case 1: lineage extinction
			-(sum(Cijk_vals[Carray_ivals .== i]) + sum(Qij_vals[Qarray_ivals .== i]) + mu[i])*u[i] +  # case 2: no event + eventual extinction
			(sum(Qij_vals[Qarray_ivals .== i] .* u[Qarray_jvals[Qarray_ivals .== i]])) + 			# case 3: change + eventual extinction
			(two * sum(Cijk_vals[Carray_ivals .== i] .* u[Carray_jvals[Carray_ivals .== i]] .* u[Carray_kvals[Carray_ivals .== i]])) 
			# case 4 & 5: speciation from i producing j,k or k,j, eventually both daughters go extinct
			# Because Carray contains all nonzero i->j,k rates, iterating once through on a particular
			# "i" does the double-summation required

		# Calculation of "D" (likelihood of tip data)
		du[n+i] = -(sum(Cijk_vals[Carray_ivals .== i]) + sum(Qij_vals[Qarray_ivals .== i]) + mu[i])*u[n+i] +  # case 1: no event
			(sum(Qij_vals[Qarray_ivals .== i] .* u[(n.+Qarray_jvals[Qarray_ivals .== i])])) + 	# case 2	
			(sum(Cijk_vals[Carray_ivals .== i] .*                                               # case 34: change + eventual extinction
				 (u[(n.+Carray_kvals[Carray_ivals .== i])].*u[Carray_jvals[Carray_ivals .== i]] 
			 .+ u[(n.+Carray_jvals[Carray_ivals .== i])].*u[Carray_kvals[Carray_ivals .== i]]) ))
  end
end


parameterized_ClaSSE_Es = (du,u,p,t) -> begin

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
	
	two = 1.0
  @inbounds for i in 1:n
		# Calculation of "E" (prob. of extinction)
		du[i] = mu[i] +                                         # case 1: lineage extinction
			-(sum(Cijk_vals[Carray_ivals .== i]) + sum(Qij_vals[Qarray_ivals .== i]) + mu[i])*u[i] +  # case 2: no event + eventual extinction
			(sum(Qij_vals[Qarray_ivals .== i] .* u[Qarray_jvals[Qarray_ivals .== i]])) + 			# case 3: change + eventual extinction
			(two * sum(Cijk_vals[Carray_ivals .== i] .* u[Carray_jvals[Carray_ivals .== i]] .* u[Carray_kvals[Carray_ivals .== i]])) 
			# case 4 & 5: speciation from i producing j,k or k,j, eventually both daughters go extinct
			# Because Carray contains all nonzero i->j,k rates, iterating once through on a particular
			# "i" does the double-summation required
  end
end

parameterized_ClaSSE_Ds = (du,u,p,t) -> begin

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
	sol_Es = p.sol_Es
	uE = p.uE
	uE = sol_Es(t)
	
	two = 1.0
  @inbounds for i in 1:n
		# Calculation of "D" (likelihood of tip data)
		du[i] = -(sum(Cijk_vals[Carray_ivals .== i]) + sum(Qij_vals[Qarray_ivals .== i]) + mu[i])*u[i] +  # case 1: no event
			(sum(Qij_vals[Qarray_ivals .== i] .* u[(Qarray_jvals[Qarray_ivals .== i])])) + 	# case 2	
			(sum(Cijk_vals[Carray_ivals .== i] .*                                               # case 34: change + eventual extinction
				 (u[(Carray_kvals[Carray_ivals .== i])].*uE[Carray_jvals[Carray_ivals .== i]] 
			 .+ u[(Carray_jvals[Carray_ivals .== i])].*uE[Carray_kvals[Carray_ivals .== i]]) ))
  end
end




 
 
# Modeled on:
# https://gitter.im/JuliaDiffEq/Lobby?at=5cbe49dbb4700e023db73d9e
# https://gist.github.com/nmatzke/b6332845747d7452b6e3b45564f460e8
# Pre-allocating the Carray_ivals .== i, Qarray_jvals[Qarray_ivals .== i
# Reduces GC (Garbage Collection) from 40% to ~5%
parameterized_ClaSSE_v5 = (du,u,p,t) -> begin

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
	
	two = 1.0
  @inbounds for i in 1:n
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]
		Qi_eq_i = p.p_TFs.Qi_eq_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]

		# Calculation of "E" (prob. of extinction)
		du[i] = mu[i] +                                         # case 1: lineage extinction
			-(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[i] +  # case 2: no event + eventual extinction
			(sum(Qij_vals[Qi_eq_i] .* u[Qj_sub_i])) + 			# case 3: change + eventual extinction
			(two * sum(Cijk_vals[Ci_eq_i] .* u[Cj_sub_i] .* u[Ck_sub_i])) 
			# case 4 & 5: speciation from i producing j,k or k,j, eventually both daughters go extinct
			# Because Carray contains all nonzero i->j,k rates, iterating once through on a particular
			# "i" does the double-summation required

		# Calculation of "D" (likelihood of tip data)
		du[n+i] = -(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[n+i] +  # case 1: no event
			(sum(Qij_vals[Qi_eq_i] .* u[(n.+Qj_sub_i)])) + 	# case 2	
			(sum(Cijk_vals[Ci_eq_i] .*                                               # case 34: change + eventual extinction
				 (u[(n.+Ck_sub_i)].*u[Cj_sub_i] 
			 .+ u[(n.+Cj_sub_i)].*u[Ck_sub_i]) ))
  end
end
 
 
parameterized_ClaSSE_Es_v5 = (du,u,p,t) -> begin

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
	
	two = 1.0
  @inbounds for i in 1:n
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]
		Qi_eq_i = p.p_TFs.Qi_eq_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]

		# Calculation of "E" (prob. of extinction)
		du[i] = mu[i] +                                         # case 1: lineage extinction
			-(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[i] +  # case 2: no event + eventual extinction
			(sum(Qij_vals[Qi_eq_i] .* u[Qj_sub_i])) + 			# case 3: change + eventual extinction
			(two * sum(Cijk_vals[Ci_eq_i] .* u[Cj_sub_i] .* u[Ck_sub_i])) 
			# case 4 & 5: speciation from i producing j,k or k,j, eventually both daughters go extinct
			# Because Carray contains all nonzero i->j,k rates, iterating once through on a particular
			# "i" does the double-summation required
  end
end

parameterized_ClaSSE_Ds_v5 = (du,u,p,t) -> begin

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
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]
		Qi_eq_i = p.p_TFs.Qi_eq_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]

		# Calculation of "D" (likelihood of tip data)
		du[i] = -(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[i] +  # case 1: no event
			(sum(Qij_vals[Qi_eq_i] .* u[Qj_sub_i])) + 	# case 2	
			(sum(Cijk_vals[Ci_eq_i] .*                                               # case 34: change + eventual extinction
				 (u[Ck_sub_i].*uE[Cj_sub_i] 
			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))
  end
end







end # end of module