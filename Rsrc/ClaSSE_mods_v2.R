# Get total LnL and branch LnL from ClaSSE output
get_classe_LnLs <- function(classe_res)
	{
	
	branch_LnL = sum(attr(classe_res, "intermediates")$lq)
	ttl_LnL = classe_res
	attributes(ttl_LnL) = NULL
	
	return(c(ttl_LnL, branch_LnL))
	}







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
		lambda_ijk_mat = matrix(0.0, ncol=k, nrow=nsum)
		is_vec = unlist(sapply(X=lambda_ijk_txt, FUN=substr, start=1, stop=2, USE.NAMES=FALSE))
		js_vec = unlist(sapply(X=lambda_ijk_txt, FUN=substr, start=3, stop=4, USE.NAMES=FALSE))
		ks_vec = unlist(sapply(X=lambda_ijk_txt, FUN=substr, start=5, stop=6, USE.NAMES=FALSE))
		lambda_ijk_mat[,1] = as.numeric(is_vec)
		lambda_ijk_mat[,2] = as.numeric(js_vec)
		lambda_ijk_mat[,3] = as.numeric(ks_vec)
		}
	
	lambda_ijk_df = as.data.frame(cbind(lambda_ijk_mat, lambda_vals))
	names(lambda_ijk_df) = c("i", "j", "k", "lambda")
	return(lambda_ijk_df)
	} # END classe_lambdas_to_df <- function(classe_params, k=3)


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
	lambda_ijk_df = classe_lambdas_to_df(classe_params=classe_params, k=3)
	
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
		
		# For this node, go through the states and sum the likes
		tmp_likes_AT_node = rep(0.0, times=k)
		for (l in 1:k) # l = ancestral node number
			{
			i_TF = lambda_ijk_df$i == l
			j_ne_k_TF = lambda_ijk_df$j != lambda_ijk_df$k
			rows_use_lambda_div2_TF = (i_TF + j_ne_k_TF) == 2
			rows_use_lambda_div1_TF = (i_TF + rows_use_lambda_div2_TF) == 1
			
			lambda_ijk_df[rows_use_lambda_div1_TF,]
			lambda_ijk_df[rows_use_lambda_div2_TF,]
			
			# Sum likes where the daughters the same
			ind = rows_use_lambda_div1_TF
			lcol = lambda_ijk_df$k[ind]
			rcol = lambda_ijk_df$j[ind]
			tmp_likes_AT_node[l] = sum(lambda_ijk_df$lambda[ind] * base_normlikes[left_desc_nodenum,lcol] * base_normlikes[right_desc_nodenum,rcol])
			# Sum likes where the daughters are NOT the same
			# (divide lambda by 2, but use twice)
			ind = rows_use_lambda_div2_TF
			lcol = lambda_ijk_df$k[ind]
			rcol = lambda_ijk_df$j[ind]
			# Left, then right
			tmp_likes_AT_node[l] = tmp_likes_AT_node[l] + sum(lambda_ijk_df$lambda[ind]/2 * base_normlikes[left_desc_nodenum,lcol] * base_normlikes[right_desc_nodenum,rcol])
			tmp_likes_AT_node[l] = tmp_likes_AT_node[l] + sum(lambda_ijk_df$lambda[ind]/2 * base_normlikes[left_desc_nodenum,rcol] * base_normlikes[right_desc_nodenum,lcol])
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



