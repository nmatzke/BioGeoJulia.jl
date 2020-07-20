
#######################################################
# Likelihood equation in the birthdeath function
# (derived by pulling apart the birthdeath() function from ape)
# This version stores all of the piece, for comparison
#######################################################
bd_liks <- function(tr, birthRate=1.0, deathRate=0.0)
	{
	ex='
	# Getting the birthRate and deathRate from
	# a = deathRate / birthRate	# relative death rate
	# r = birthRate - deathRate	# net diversification rate
	BD =  birthdeath(tr)
	BD
	names(BD)

	# Calculate the birthRate and deathRate from the outputs
	x1 = unname(BD$para["d/b"])
	x2 = unname(BD$para["b-d"])
	deathRate = (x2*x1) / (1-x1)
	birthRate = deathRate+x2
	c(birthRate, deathRate)
	'
	
	
	a = deathRate / birthRate	# relative death rate
	r = birthRate - deathRate	# net diversification rate

	N <- length(tr$tip.label)
	nb_node = tr$Nnode - 1
	sptimes <- c(NA, branching.times(tr)) # NA so the number of times equals number of tips?
	x = sptimes
	# a = "d/b"
	# r = "b-d"

	# dev <- function(a=0.1, r=0.2, N, x, return_deviance=FALSE)
	# 	{
	if (r < 0 || a > 1)
		{
		return(1e+100)
		}
	
	lnl_topology = lfactorial(tr$Nnode)
	lnl_numBirths = nb_node * log(r)
	lnl_Births_above_root = r * sum(sptimes[3:N])
	
	lnl_numtips_wOneMinusDeathRate = N * log(1 - a)
	# Interpretation: more tips are less likely, if relativeDeathRate is >0
	# If relativeDeathRate = 1, a=0, and lnl=-Inf... 
	#    CLEARLY WRONG EXCEPT IN A MAXIMUM LIKELIHOOD CONTEXT!!!
	# If relativeDeathRate = 0, a=0, and lnl=0, i.e. any number of tips is equiprobable
	
	lnl_branching_times = -2 * sum(log(exp(r * sptimes[2:N]) - a))
	# For each observed branching,
	# netDiversificationRate * timeOfEvent <- take exponent of that ; this means recorded events are less likely in the past
	# <- subtract "a", a constant (relativeDeathRate)
	#
	# This would be a straight likelihood as:
	# 1/
	# (exp(r*branching_time)-a)^2
	#
	# Sum the logs of these
	#
	# If deathRate = 0
	# lnl_branching_times = -2 * sum(log(exp(birthRate*sptimes[2:N]) - 0))
	# lnl_branching_times = -2 * sum(log( exp(birthRate*sptimes[2:N]) )
	# lnl_branching_times = -2 * sum( birthRate*sptimes[2:N] )
	#
	# Note: sum(X) = 9 = total branchlength of tr
	# In BD:
	# -2*sum(sptimes[2:N]) = -12
	# sum(sptimes[3:N]) = 3
	# So, lnl_branching_times + lnl_Births_above_root = yule's -lambda * X
	lnL = lnl_topology + lnl_numBirths + lnl_Births_above_root + lnl_numtips_wOneMinusDeathRate + lnl_branching_times
	dev =  -2 * lnL
	
	bd = NULL
	bd$tr = tr
	bd$birthRate = birthRate
	bd$deathRate = deathRate
	bd$relative_deathRate = a
	bd$net_diversification_rate = r
	bd$dev = dev
	bd$lnl_topology = lnl_topology
	bd$lnl_numBirths = lnl_numBirths
	bd$lnl_Births_above_root = lnl_Births_above_root
	bd$lnl_numtips_wOneMinusDeathRate = lnl_numtips_wOneMinusDeathRate
	bd$lnl_branching_times = lnl_branching_times
	bd$lnL = lnL
	
	return(bd)
	}






# Get total LnL and branch LnL from ClaSSE output
get_classe_LnLs <- function(classe_res)
	{
	
	branch_LnL = sum(attr(classe_res, "intermediates")$lq)
	ttl_LnL = classe_res
	attributes(ttl_LnL) = NULL
	
	return(c(ttl_LnL, branch_LnL))
	}



BGBres_into_classe_params <- function(res, classe_params, birthRate=0.2)
	{
	mats = get_Qmat_COOmat_from_res(res, numstates=ncol(res$ML_marginal_prob_each_state_at_branch_top_AT_node), include_null_range=TRUE, max_range_size=res$inputs$max_range_size, timeperiod_i=1)
	numstates = length(mats$states_list)
	include_null_range = res$inputs$include_null_range
	
	# BioGeoBEARS Cevent weights into DF
	Carray_df = get_Cevent_probs_df_from_mats(mats, include_null_range=include_null_range)
	
	# Diversitree params into a df
	lambda_ijk_df = classe_lambdas_to_df(classe_params, k=numstates)
	
	# Put the Cevent_weights * birthRate into diversitree table
	lambda_ijk_df = get_Cevent_lambdas_from_BGB_Carray(lambda_ijk_df, Carray_df, birthRate=birthRate)
	lambda_ijk_df[lambda_ijk_df$lambda != 0,]

	# Convert the diversitree table back to classe_params
	lambdas_to_put_in_params = rownames(lambda_ijk_df[lambda_ijk_df$lambda != 0,])
	indices_in_classe_params = match(x=lambdas_to_put_in_params, table=names(classe_params))
	classe_params[indices_in_classe_params] = lambda_ijk_df[lambdas_to_put_in_params,"lambda"]
	classe_params[153]
	classe_params["lambda020202"]
	classe_params["lambda060203"]
	
	
	# Now do the Qmat
	Qij_df = classe_Qs_to_df(classe_params, k=numstates)
	Qmat = mats$Qmat
	Qij_df = get_Qijs_from_BGB_Qarray(Qij_df, Qmat)
	Qs_to_put_in_params = rownames(Qij_df[Qij_df$q != 0,])
	indices_in_classe_params = match(x=Qs_to_put_in_params, table=names(classe_params))
	classe_params[indices_in_classe_params] = Qij_df[Qs_to_put_in_params,"q"]
	
	
	# And the mus	
	# Because mu (the lineage extinction rate) is always 0.0 in BioGeoBEARS models
	# (i.e., a pure-birth Yule process is being assumed)
	classe_params[grepl(pattern="mu", x=names(classe_params))] = 0.0 
	
	return(classe_params)
	} # BGBres_into_classe_params <- function(res, classe_params, birthRate=0.2)

# k= number of states
# Note that diversitree classe includes only e.g. q123, not q 132
classe_Qs_to_df <- function(classe_params, k=3)
	{
	ex='
	k = 3
	classe_params = c(lambda111 = 0.2, lambda112 = 0, lambda113 = 0, lambda122 = 0, 
lambda123 = 0, lambda133 = 0, lambda211 = 0, lambda212 = 0, lambda213 = 0, 
lambda222 = 0.2, lambda223 = 0, lambda233 = 0, lambda311 = 0, 
lambda312 = 0.0666666666666667, lambda313 = 0.0666666666666667, 
lambda322 = 0, lambda323 = 0.0666666666666667, lambda333 = 0, 
mu1 = 0.1, mu2 = 0.1, mu3 = 0.1, q12 = 0, q13 = 0, q21 = 0, q23 = 0, 
q31 = 0, q32 = 0)
	
	Qij_df = classe_qs_to_df(classe_params=classe_params, k=3)
	Qij_df
	' # end example
	
	# Number of qs per state (e.g., for 3 states, this is 6 qs per state
	#nsum <- k * (k + 1)/2

	# Convert the classe_params qs to a table
	Q_vals = classe_params[grepl(pattern="q", x=names(classe_params))]
	Q_names = names(Q_vals)
	
	# Get i, j, k indices
	Qij_txt = gsub(pattern="q", replacement="", x=Q_names)
	
	if (k <= 9)
		{
		ijs_vector = unlist(sapply(X=Qij_txt, FUN=strsplit, split="", USE.NAMES=FALSE))
		Qij_mat = matrix(as.numeric(ijs_vector), ncol=k, byrow=TRUE)
		Qij_mat
		} else {
		# Split by 2
		is_vec = unlist(sapply(X=Qij_txt, FUN=substr, start=1, stop=2, USE.NAMES=FALSE))
		js_vec = unlist(sapply(X=Qij_txt, FUN=substr, start=3, stop=4, USE.NAMES=FALSE))
		Qij_mat = cbind(is_vec, js_vec)
		}
	
	Qij_df = as.data.frame(cbind(Qij_mat, Q_vals), stringsAsFactors=FALSE)
	names(Qij_df) = c("i", "j", "q")
	# Convert qs to numeric
	Qij_df$q = as.numeric(as.character(Qij_df$q))
	return(Qij_df)
	} # END classe_Qs_to_df <- function(classe_params, k=3)


get_Qijs_from_BGB_Qarray <- function(Qij_df, Qmat)
	{
	for (r in 1:nrow(Qij_df))
		{
		# r = 153
		# Qij_df[153,]
		#               i  j  k lambda
		# lambda020202 02 02 02      0
		i = as.numeric(Qij_df$i[r])
		j = as.numeric(Qij_df$j[r])
	
		# Don't include the diagonals from Q
		if (i == j)
			{
			next()
			}

		Qij_df$q[r] = Qij_df$q[r] + Qmat[i,j]
		} # END for (r in nrow(Qij_df))

	Qij_df[Qij_df$q != 0,]
	return(Qij_df)
	} # END get_Qijs_from_BGB_Qarray <- function(Qij_df, Qmat, birthRate=1.0)


# k= number of states
# Note that diversitree classe includes only e.g. lambda123, not lambda 132
classe_lambdas_to_df <- function(classe_params, k=3)
	{
	ex='
	k = 3
	classe_params = c(lambda111 = 0.2, lambda112 = 0, lambda113 = 0, lambda122 = 0, 
lambda123 = 0, lambda133 = 0, lambda211 = 0, lambda212 = 0, lambda213 = 0, 
lambda222 = 0.2, lambda223 = 0, lambda233 = 0, lambda311 = 0, 
lambda312 = 0.0666666666666667, lambda313 = 0.0666666666666667, 
lambda322 = 0, lambda323 = 0.0666666666666667, lambda333 = 0, 
mu1 = 0.1, mu2 = 0.1, mu3 = 0.1, q12 = 0, q13 = 0, q21 = 0, q23 = 0, 
q31 = 0, q32 = 0)
	
	lambda_ijk_df = classe_lambdas_to_df(classe_params=classe_params, k=3)
	lambda_ijk_df
	' # end example
	
	# Number of lambdas per state (e.g., for 3 states, this is 6 lambdas per state
	nsum <- k * (k + 1)/2

	# Convert the classe_params lambdas to a table
	lambda_vals = classe_params[grepl(pattern="lambda", x=names(classe_params))]
	lambda_names = names(lambda_vals)
	
	# Get i, j, k indices
	lambda_ijk_txt = gsub(pattern="lambda", replacement="", x=lambda_names)
	
	if (k <= 9)
		{
		ijks_vector = unlist(sapply(X=lambda_ijk_txt, FUN=strsplit, split="", USE.NAMES=FALSE))
		lambda_ijk_mat = matrix(as.numeric(ijks_vector), ncol=k, byrow=TRUE)
		lambda_ijk_mat
		} else {
		# Split by 2
		is_vec = unlist(sapply(X=lambda_ijk_txt, FUN=substr, start=1, stop=2, USE.NAMES=FALSE))
		js_vec = unlist(sapply(X=lambda_ijk_txt, FUN=substr, start=3, stop=4, USE.NAMES=FALSE))
		ks_vec = unlist(sapply(X=lambda_ijk_txt, FUN=substr, start=5, stop=6, USE.NAMES=FALSE))
		lambda_ijk_mat = cbind(is_vec, js_vec, ks_vec)
		}
	
	lambda_ijk_df = as.data.frame(cbind(lambda_ijk_mat, lambda_vals), stringsAsFactors=FALSE)
	names(lambda_ijk_df) = c("i", "j", "k", "lambda")
	# Convert lambdas to numeric
	lambda_ijk_df$lambda = as.numeric(as.character(lambda_ijk_df$lambda))
	return(lambda_ijk_df)
	} # END classe_lambdas_to_df <- function(classe_params, k=3)


get_Cevent_lambdas_from_BGB_Carray <- function(lambda_ijk_df, Carray_df, birthRate=1.0)
	{
	for (r in 1:nrow(lambda_ijk_df))
		{
		# r = 153
		# lambda_ijk_df[153,]
		#               i  j  k lambda
		# lambda020202 02 02 02      0
		i = as.numeric(lambda_ijk_df$i[r])
		j = as.numeric(lambda_ijk_df$j[r])
		k = as.numeric(lambda_ijk_df$k[r])
	
		iTF = Carray_df$i == i
		jTF = Carray_df$j == j
		kTF = Carray_df$k == k
		TF = (iTF + jTF + kTF) == 3
		if (sum(TF) == 1)
			{
			lambda_ijk_df$lambda[r] = lambda_ijk_df$lambda[r] + (Carray_df$prob[TF] * birthRate)
			}
	
		# Don't repeat the search when j==k (left and right states are the same)
		if (j == k)
			{
			next()
			}
	
		# Repeat, switching j and k states
		jTF = Carray_df$j == k
		kTF = Carray_df$k == j
		TF = (iTF + jTF + kTF) == 3
		if (sum(TF) == 1)
			{
			lambda_ijk_df$lambda[r] = lambda_ijk_df$lambda[r] + (Carray_df$prob[TF] * birthRate)
			} # END if (sum(TF) == 1)
		} # END for (r in nrow(lambda_ijk_df))

	lambda_ijk_df[lambda_ijk_df$lambda != 0,]
	return(lambda_ijk_df)
	} # END get_Cevent_lambdas_from_BGB_Carray <- function(lambda_ijk_df, Carray_df)



# Calculate the sum of the log computed likelihoods at each node
# i.e., the likelihoods of the speciation events, assuming normalized
# likelihoods at the branch bottoms above each node
# "base" is "t(base)", actually
get_sum_log_computed_likes_at_each_node <- function(tr, base, lq, classe_params)
	{
	# Number of lambdas per state (e.g., for 3 states, this is 6 lambdas per state
	k = ncol(base) / 2
	nsum <- k * (k + 1)/2
	
	Ds_cols = (k+1):(2*k)
	base_likes = apply(X=base[,Ds_cols], MARGIN=2, FUN="*", exp(lq))
	base_normlikes = base_likes / rowSums(base_likes)
	
	# Get a data.frame tabulating the lambdas
	lambda_ijk_df = classe_lambdas_to_df(classe_params=classe_params, k=k)
	lambda_ijk_df$i = as.numeric(as.character(lambda_ijk_df$i))
	lambda_ijk_df$j = as.numeric(as.character(lambda_ijk_df$j))
	lambda_ijk_df$k = as.numeric(as.character(lambda_ijk_df$k))
	lambda_ijk_df$lambda = as.numeric(as.character(lambda_ijk_df$lambda))
	
	
	# Go through the ancestral states
	computed_likelihoods_at_each_node_just_before_speciation = matrix(0.0, nrow=nrow(base), ncol=k)
	
	# Reorder the edge matrix into pruningwise order
	# This is CRUCIAL!!
	tr2 <- reorder(tr, "pruningwise")
	num_internal_nodes = tr$Nnode
	
	# DEFINE DOWNPASS THROUGH THE BRANCHES	
	i = 1
	edges_to_visit = seq(from=1, by=2, length.out=num_internal_nodes)
	
	for (i in edges_to_visit)
		{
		# Get the node numbers at the tips of these two edges		
		j = i+1
		left_desc_nodenum <- tr2$edge[i, 2]
		right_desc_nodenum <- tr2$edge[j, 2]
		# And for the ancestor edge (i or j shouldn't matter, should produce the same result!!!)
		anc <- tr2$edge[i, 1]
		
		# For this anc node, go through the states and sum the likes
		tmp_likes_AT_node = rep(0.0, times=k)
		for (l in 1:k) # l = ancestral state number
			{
			i_TF = lambda_ijk_df$i == l
			j_ne_k_TF = lambda_ijk_df$j != lambda_ijk_df$k
			rows_use_lambda_div2_TF = (i_TF + j_ne_k_TF) == 2
			rows_use_lambda_div1_TF = (i_TF + rows_use_lambda_div2_TF) == 1
			lambda_ijk_df[rows_use_lambda_div1_TF,]
			lambda_ijk_df[rows_use_lambda_div2_TF,]
			
			# Skip e.g. null-range states
			if (sum(rows_use_lambda_div1_TF) > 0)
				{
				# Sum likes where the daughters the same
				ind = rows_use_lambda_div1_TF
				lcol = lambda_ijk_df$j[ind]
				rcol = lambda_ijk_df$k[ind]
				tmp_likes_AT_node[l] = sum(lambda_ijk_df$lambda[ind] * base_normlikes[left_desc_nodenum,lcol] * base_normlikes[right_desc_nodenum,rcol])
				} # END if (sum(rows_use_lambda_div1_TF) > 0)

			if (sum(rows_use_lambda_div2_TF) > 0)
				{
				# Sum likes where the daughters are NOT the same
				# (divide lambda by 2, but use twice)
				ind = rows_use_lambda_div2_TF
				lcol = lambda_ijk_df$j[ind]
				rcol = lambda_ijk_df$k[ind]
				# Left, then right
				tmp_likes_AT_node[l] = tmp_likes_AT_node[l] + sum(lambda_ijk_df$lambda[ind]/2 * base_normlikes[left_desc_nodenum,lcol] * base_normlikes[right_desc_nodenum,rcol])
				tmp_likes_AT_node[l] = tmp_likes_AT_node[l] + sum(lambda_ijk_df$lambda[ind]/2 * base_normlikes[left_desc_nodenum,rcol] * base_normlikes[right_desc_nodenum,lcol])
				} # END if (sum(rows_use_lambda_div1_TF) > 0)
			} # END for (l in 1:k)
		
		computed_likelihoods_at_each_node_just_before_speciation[anc,] = tmp_likes_AT_node
		} # END for (i in edges_to_visit) 
	
	return(computed_likelihoods_at_each_node_just_before_speciation)
	} # END get_sum_log_computed_likes_at_each_node <- function(tr, base, lq, classe_params)






# Convert a "res" object from diversitree 
# ClaSSE to a 
# BioGeoBEARS-like set of matrices
claSSE_res_to_prt <- function(res, tr, classe_params)
	{
	ex='
	res1t = structure(-6.22837508651605, intermediates = list(init = structure(c(0, 
0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0.153452947030521, 
0.153452947030521, 0.153452947030521, 0, 0, 0.0326226539925451, 
0.0868935659418745, 0.0868935659418745, 0.0868935659418745, 0, 
0, 0.0333321365204794), .Dim = 6:5), base = structure(c(0.0868935659418745, 
0.0868935659418745, 0.0868935659418745, 0, 0.994007973162888, 
0.00599202683711178, 0.0868935659418746, 0.0868935659418746, 
0.0868935659418745, 0.994007973162888, 0, 0.00599202683711177, 
0.153452946974297, 0.153452946974297, 0.153452946974297, 0, 0.978679619776353, 
0.0213203802236467, NA, NA, NA, NA, NA, NA, 0.153452947030521, 
0.153452947030521, 0.153452947030521, 0, 0, 1), .Dim = 6:5), 
    lq = c(-0.275795607421526, -0.275795607421526, -0.511628046092027, 
    0, -3.68502440109251), vals = c(0.153452947030521, 0.153452947030521, 
    0.153452947030521, 0, 0, 0.0326226539925451), branchLnL = -4.74824366202759, 
    root.p = c(0, 0, 1)), vals = c(0.153452947030521, 0.153452947030521, 
0.153452947030521, 0, 0, 0.0326226539925451))
	res = res1t
	
	trstr = "((chimp:1,human:1):1,gorilla:2);"
	tr = read.tree(file="", text=trstr)
	
	k = 3
	classe_params = c(lambda111 = 0.2, lambda112 = 0, lambda113 = 0, lambda122 = 0, 
lambda123 = 0, lambda133 = 0, lambda211 = 0, lambda212 = 0, lambda213 = 0, 
lambda222 = 0.2, lambda223 = 0, lambda233 = 0, lambda311 = 0, 
lambda312 = 0.0666666666666667, lambda313 = 0.0666666666666667, 
lambda322 = 0, lambda323 = 0.0666666666666667, lambda333 = 0, 
mu1 = 0.1, mu2 = 0.1, mu3 = 0.1, q12 = 0, q13 = 0, q21 = 0, q23 = 0, 
q31 = 0, q32 = 0)
	
	lambda_ijk_df = classe_lambdas_to_df(classe_params=classe_params, k=3)
	lambda_ijk_df
	
	classe_res_dfs = claSSE_res_to_prt(res, tr, classe_params)
	classe_res_dfs
	names(classe_res_dfs)
	' # END ex
	
	# Branch top values ("initial")
	init = t(attr(res,"intermediates")$init)
	
	# Branch bottom values ("base")
	base = t(attr(res,"intermediates")$base)
	
	# lqs = log-likelihoods at each branch bottom
	lq = attr(res,"intermediates")$lq
	q = exp(attr(res,"intermediates")$lq)
	
	vals = t(attr(res2, "intermediates")$vals)	# Es and Ds at the root	
	E_indices = 1:k
	d_root_orig = vals[-E_indices]							# Original D likelihoods at root

	# Assumes bifurcating tree
	numstates = ncol(init) / 2
	rootnode = length(tr$tip.label) + 1
	numnodes = nrow(init) # internal plus tip nodes
#	numTips = (nrow(init) + 1) / 2
	numTips = length(tr$tip.label)
	numInternal = numTips - 1
	
	# Intitializing
	Es_atNode_branchTop = matrix(data=0, ncol=numstates, nrow=numnodes)
	Es_atNode_branchBot = matrix(data=0, ncol=numstates, nrow=numnodes) 
	likes_at_each_nodeIndex_branchTop = matrix(data=0, ncol=numstates, nrow=numnodes)
	likes_at_each_nodeIndex_branchBot = matrix(data=0, ncol=numstates, nrow=numnodes) 
	normlikes_at_each_nodeIndex_branchTop = matrix(data=0, ncol=numstates, nrow=numnodes)
	normlikes_at_each_nodeIndex_branchBot = matrix(data=0, ncol=numstates, nrow=numnodes)
	
	# Filling in
	Ecols = 1:numstates
	Dcols = (numstates+1):(2*numstates)
	Es_atNode_branchTop = init[,Ecols]
	Es_atNode_branchBot = base[,Ecols]
	likes_at_each_nodeIndex_branchTop = init[,Dcols]
	normlikes_at_each_nodeIndex_branchBot = base[,Dcols]
	normlikes_at_each_nodeIndex_branchTop = likes_at_each_nodeIndex_branchTop / rowSums(likes_at_each_nodeIndex_branchTop)
	likes_at_each_nodeIndex_branchBot = normlikes_at_each_nodeIndex_branchBot * q
	
	# Likelihoods just after nodeOp, before normalization
	computed_likelihoods_at_each_node_just_before_speciation = get_sum_log_computed_likes_at_each_node(tr, base, lq, classe_params)
	log(rowSums(computed_likelihoods_at_each_node_just_before_speciation))
	TF = is.finite(log(rowSums(computed_likelihoods_at_each_node_just_before_speciation)))
	sum_log_nodelikes = sum(lq) + sum(log(rowSums(computed_likelihoods_at_each_node_just_before_speciation))[TF])
sum_log_nodelikes

	# LnLs1
	# If root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
	d_root_orig = likes_at_each_nodeIndex_branchTop[rootnode,]						# Original D likelihoods at root
	root.p = d_root_orig/sum(d_root_orig)
	loglik_s1 = log(sum(root.p * d_root_orig)) + sum(lq)
	loglik_s1

	# LnLs1t
	# If root=ROOT.OBS, root.p=NULL, condition.surv=TRUE
	root.p = d_root_orig/sum(d_root_orig)
	# MuSSE/ClaSSE
	pars = classe_params
	nsum <- k * (k + 1)/2
	lambda <- colSums(matrix(pars[1:(nsum * k)], nrow = nsum))
	i <- seq_len(k)
	e.root <- Es_atNode_branchTop[rootnode,]
	d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
	loglik_s1t = log(sum(root.p * d.root)) + sum(lq)
	loglik_s1t

	R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = sum_log_nodelikes


# R_result_branch_lnL = -4.748244
# R_result_total_LnLs1 = -8.170992
# R_result_total_LnLs1t = -6.228375
# R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -11.57223

	R_result_branch_lnL = sum(lq)
	R_result_total_LnLs1 = loglik_s1
	R_result_total_LnLs1t = loglik_s1t
	R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = sum_log_nodelikes
	
	
	classe_res_dfs = NULL
	classe_res_dfs$init = t(attr(res,"intermediates")$init)
	classe_res_dfs$base = t(attr(res,"intermediates")$init)
	classe_res_dfs$lq = lq
	classe_res_dfs$Es_atNode_branchTop = Es_atNode_branchTop
	classe_res_dfs$Es_atNode_branchBot = Es_atNode_branchBot
	classe_res_dfs$likes_at_each_nodeIndex_branchTop = likes_at_each_nodeIndex_branchTop
	classe_res_dfs$likes_at_each_nodeIndex_branchBot = likes_at_each_nodeIndex_branchBot
	classe_res_dfs$normlikes_at_each_nodeIndex_branchTop = normlikes_at_each_nodeIndex_branchTop
	classe_res_dfs$normlikes_at_each_nodeIndex_branchBot = normlikes_at_each_nodeIndex_branchBot
	classe_res_dfs$computed_likelihoods_at_each_node_just_before_speciation = computed_likelihoods_at_each_node_just_before_speciation
	classe_res_dfs$R_result_branch_lnL = R_result_branch_lnL
	classe_res_dfs$R_result_total_LnLs1 = R_result_total_LnLs1
	classe_res_dfs$R_result_total_LnLs1t = R_result_total_LnLs1t
	classe_res_dfs$R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
	
	names(classe_res_dfs)
	return(classe_res_dfs)
	} # END claSSE_res_to_prt <- function(res, tr, classe_params)





#diversitree:::all.branches.matrix
all.branches.matrix <- function(pars, cache, initial.conditions, branches, preset = NULL) 
{
    len <- cache$len
    depth <- cache$depth
    children <- cache$children
    #order <- cache$order[-length(cache$order)]
    order <- cache$order
    root <- cache$root
    n <- length(len)
    lq <- rep(0, n)
    n.tip <- cache$n.tip
    y <- cache$y
    branch.init <- branch.base <- matrix(NA, cache$info$ny, n)
    if (!is.null(preset)) {
        lq[preset$target] <- preset$lq
        branch.base[, preset$target] <- preset$base
    }
    if (is.null(names(y))) {
        for (x in y) {
            if (!is.null(x)) {
                idx <- x$target
                branch.init[, idx] <- x$y
                ans <- branches(x$y, x$t, pars, 0, idx)
                lq[idx] <- ans[[1]]
                branch.base[, idx] <- ans[[2]]
            }
        }
    }
    else {
        tip.t <- y$t
        tip.target <- y$target
        tip.y <- branch.init[, tip.target] <- y$y
        for (i in seq_along(tip.t)) {
            idx <- tip.target[i]
            ans <- branches(tip.y[, i], tip.t[i], pars, 0, idx)
            lq[idx] <- ans[[1]]
            branch.base[, idx] <- ans[[2]]
        }
    }
    
    # NJM ADD additional things to log
    ans_lists = NULL
    branches_function_list = NULL
    countval = 0
    for (i in order) {
    	# Count so you can treat the root specially
    	countval = countval + 1
        y.in <- initial.conditions(branch.base[, children[i, ]], pars, depth[i], i)
    	
    	# Normalize to sum to 1, if you are at the root state
    	# THIS CONVERTS LIKELIHOODS TO ROOT.OBS
    	#if (countval != length(cache$order))
    	#	{
    	#	y.in = y.in / sum(y.in)
	    #   }
   		#y.in = y.in / sum(y.in)
        if (!all(is.finite(y.in))) 
            stop("Bad initial conditions: calculation failure along branches?")
        branch.init[, i] <- y.in
        ans <- branches(y.in, len[i], pars, depth[i], i)
        lq[i] <- ans[[1]]
        branch.base[, i] <- ans[[2]]

	    # NJM ADD additional things to log
	    ans_lists[[i]] = ans
        branches_function_list[[i]] = body(branches)
    }
    y.in <- initial.conditions(branch.base[, children[root, ]], 
        pars, depth[root], root)
    branch.init[, root] <- y.in
    
    # This sets:
    # 
    # $init
    # $lq
    # $vals
    # $ans_lists
    # 
    list(init = branch.init, base = branch.base, lq = lq, vals = y.in, ans_lists=ans_lists, branches_function_list=branches_function_list)
}



