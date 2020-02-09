module TreePass
using BioGeoJulia.TrUtils 
using BioGeoJulia.SSEs 
using DataFrames
using PhyloNetworks		# for e.g. readTopology()
using Dates						# for e.g. DateTime, Dates.now()
using Distributed			# for e.g. @spawn
using Random					# for MersenneTwister()
using DifferentialEquations # for ODEProblem
export get_nodenumbers_above_node, get_postorder_nodenumbers_above_node, initialize_edgematrix, get_pruningwise_postorder_edgematrix, get_LR_uppass_edgematrix, get_LR_downpass_edgematrix, get_LR_uppass_nodeIndexes, get_LR_downpass_nodeIndexes, get_Rnodenums, get_nodeIndex_PNnumber, get_nodeIndex_from_PNnumber, prt, get_taxa_descending_from_each_node, isTip_TF, get_NodeIndexes_from_edge, get_NodeIndex_df_by_tree_edges, get_node_heights, get_node_ages, Res, construct_Res, count_nodes_finished, nodeOp, branchOp_example, branchOp_ClaSSE_Ds_v5, branchOp, setup_inputs_branchOp_ClaSSE_Ds_v5, countloop, iterative_downpass!, iterative_downpass_nonparallel!






# Get the 2 nodeIndexes descending from a node; iterates for uppass
# (assumes a binary, PhyloNetwork, rooted tree)
# Reverse for downpass

"""
using DataFrames
using PhyloNetworks
using PhyloPlots
include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr

indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

get_nodenumbers_above_node(tr, tr.root, indexNum_table=indexNum_table)
"""

function get_nodenumbers_above_node(tr, rootnodenum, nodeIndex_array, iterNum; indexNum_table)
  if (tr.node[rootnodenum].leaf != true)
  	# Left descendant edge
  	one_edge = tr.node[rootnodenum].edge[1]
  	anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
  	left_dec_PNnumber = anc_decPNnumbers[2]
		left_dec_nodeIndex = get_nodeIndex_from_PNnumber(left_dec_PNnumber, indexNum_table=indexNum_table)
		
  	one_edge = tr.node[rootnodenum].edge[2]
  	anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
  	right_dec_PNnumber = anc_decPNnumbers[2]
		right_dec_nodeIndex = get_nodeIndex_from_PNnumber(right_dec_PNnumber, indexNum_table=indexNum_table)
  	
  	# Then, iterate through left and right clades
  	#println(rootnodenum)
  	nodeIndex_array[iterNum] = rootnodenum
  	#print(nodeIndex_array)
  	iterNum = iterNum + 1
  	(nodeIndex_array, iterNum) = get_nodenumbers_above_node(tr, right_dec_nodeIndex, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
  	(nodeIndex_array, iterNum) = get_nodenumbers_above_node(tr, left_dec_nodeIndex, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
  	return (nodeIndex_array, iterNum)
  else
  	#println(rootnodenum)
  	nodeIndex_array[iterNum] = rootnodenum
  	#print(nodeIndex_array)
  	iterNum = iterNum + 1
  	return (nodeIndex_array, iterNum)
  end
end


function get_postorder_nodenumbers_above_node(tr, rootnodenum, nodeIndex_array, iterNum; indexNum_table)
  if (tr.node[rootnodenum].leaf != true)
  	# Left descendant edge
  	one_edge = tr.node[rootnodenum].edge[1]
  	anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
  	left_dec_PNnumber = anc_decPNnumbers[2]
		left_dec_nodeIndex = get_nodeIndex_from_PNnumber(left_dec_PNnumber, indexNum_table=indexNum_table)
		
  	one_edge = tr.node[rootnodenum].edge[2]
  	anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
  	right_dec_PNnumber = anc_decPNnumbers[2]
		right_dec_nodeIndex = get_nodeIndex_from_PNnumber(right_dec_PNnumber, indexNum_table=indexNum_table)
  	
  	# Then, iterate through left and right clades
  	#println(rootnodenum)
  	iterNum = iterNum + 1
  	nodeIndex_array[iterNum] = rootnodenum
  	(nodeIndex_array, iterNum) = get_postorder_nodenumbers_above_node(tr, right_dec_nodeIndex, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
  	(nodeIndex_array, iterNum) = get_postorder_nodenumbers_above_node(tr, left_dec_nodeIndex, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
  	#print(nodeIndex_array)

  	return (nodeIndex_array, iterNum)
  else
  	#println(rootnodenum)
  	iterNum = iterNum + 1
  	nodeIndex_array[iterNum] = rootnodenum
  	#print(nodeIndex_array)
  	return (nodeIndex_array, iterNum)
  end
end


function initialize_edgematrix(tr)
	ancNodeIndex = collect(repeat([0], 2*(tr.numNodes-tr.numTaxa)))
  decNodeIndex = collect(repeat([0], 2*(tr.numNodes-tr.numTaxa)))
  edgematrix = hcat(ancNodeIndex,decNodeIndex)
  return(edgematrix)
end


# Returns an edgematrix, with
# 1st column gives ancestral nodeIndex for the edge
# 2nd column gives descendant nodeIndex for the edge
# The edges are in order:
#   pair of edges descending from the root
#   pairs following up the right branch above root
#   pairs following up the left branch above root
#
# Iterate backwards for a postorder traversal in "pruningwise" order 
# (as in APE phylo.reorder "pruningwise")
#

"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
great_ape_newick_string = "(((human:6,(chimp1:0.5,bonobo:0.5):5.5):1,gorilla:7):5,orangutan:12);"
great_ape_newick_string = "(((human:6,(chimp1:0.5,bonobo:0.5):5.5):1,gorilla:7):5,(orang1:1,orang2:1):11);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root





rootnodenum = tr.root
#tipNodeIndex_array = collect(repeat([0], length(tr.numTaxa)))
NodeIndex_array = collect(repeat([0], 2*(tr.numNodes-tr.numTaxa)));
iterNum = 0;
indexNum_table = get_nodeIndex_PNnumber(tr);
uppass_edgematrix = initialize_edgematrix(tr)
res = get_pruningwise_postorder_edgematrix(tr, rootnodenum);
uppass_edgematrix = res[1]

get_LR_uppass_edgematrix(tr)
get_LR_downpass_edgematrix(tr)

"""

function get_pruningwise_postorder_edgematrix(tr, rootnodenum, iterNum=0; edgematrix=initialize_edgematrix(tr), indexNum_table=get_nodeIndex_PNnumber(tr))
  if (tr.node[rootnodenum].leaf != true)
  	# Left descendant edge
  	one_edge = tr.node[rootnodenum].edge[1]
  	anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
  	left_anc_PNnumber = anc_decPNnumbers[1]
  	left_dec_PNnumber = anc_decPNnumbers[2]
		left_anc_nodeIndex = get_nodeIndex_from_PNnumber(left_anc_PNnumber, indexNum_table=indexNum_table)
		left_dec_nodeIndex = get_nodeIndex_from_PNnumber(left_dec_PNnumber, indexNum_table=indexNum_table)
		
  	one_edge = tr.node[rootnodenum].edge[2]
  	anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
  	right_anc_PNnumber = anc_decPNnumbers[1]
  	right_dec_PNnumber = anc_decPNnumbers[2]
		right_anc_nodeIndex = get_nodeIndex_from_PNnumber(right_anc_PNnumber, indexNum_table=indexNum_table)
		right_dec_nodeIndex = get_nodeIndex_from_PNnumber(right_dec_PNnumber, indexNum_table=indexNum_table)
  	
  	# Then, iterate through left and right clades
  	#println(rootnodenum)
  	iterNum = iterNum + 1
  	edgematrix[iterNum,1] = right_anc_nodeIndex
  	edgematrix[iterNum,2] = right_dec_nodeIndex
  	iterNum = iterNum + 1
  	edgematrix[iterNum,1] = left_anc_nodeIndex
  	edgematrix[iterNum,2] = left_dec_nodeIndex
  	
  	(edgematrix, iterNum) = get_pruningwise_postorder_edgematrix(tr, right_dec_nodeIndex, iterNum, edgematrix=edgematrix, indexNum_table=indexNum_table)
  	(edgematrix, iterNum) = get_pruningwise_postorder_edgematrix(tr, left_dec_nodeIndex, iterNum, edgematrix=edgematrix, indexNum_table=indexNum_table)
  	#print(nodeIndex_array)

  	return (edgematrix, iterNum)
  else
  	#println(rootnodenum)
  	#iterNum = iterNum + 1
  	#nodeIndex_array[iterNum] = rootnodenum
  	#print(nodeIndex_array)
  	return (edgematrix, iterNum)
  end
	
	# Shouldn't get here
	return ("get_pruningwise_postorder_edgematrix(): Shouldn't get here!")
end




# Get the node indexes for an uppass from the root to the tips
# Reverse for downpass (postorder)
# Uses get_nodenumbers_above_node() recursively to iterate up tree

"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

uppass_edgematrix = get_LR_uppass_edgematrix(tr)
uppass_edgematrix
"""

function get_LR_uppass_edgematrix(tr)
	rootnodenum = tr.root
	iterNum = 0
	edgematrix = initialize_edgematrix(tr)
	indexNum_table = get_nodeIndex_PNnumber(tr)

	res = get_pruningwise_postorder_edgematrix(tr, rootnodenum, iterNum, edgematrix=edgematrix, indexNum_table=indexNum_table)
	uppass_edgematrix = res[1]
	return(uppass_edgematrix)
end


# Get the node indexes for an downpass from the root to the tips
# Reverse of get_LR_uppass_nodeIndexes()
# 

"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

uppass_edgematrix = get_LR_uppass_edgematrix(tr)
uppass_edgematrix

downpass_edgematrix = get_LR_downpass_edgematrix(tr)
downpass_edgematrix
"""
function get_LR_downpass_edgematrix(tr)
	rootnodenum = tr.root
	iterNum = 0
	edgematrix = initialize_edgematrix(tr)
	indexNum_table = get_nodeIndex_PNnumber(tr)

	res = get_pruningwise_postorder_edgematrix(tr, rootnodenum, iterNum, edgematrix=edgematrix, indexNum_table=indexNum_table)
	uppass_edgematrix = res[1]
	
	# Reverse
	numrows = size(uppass_edgematrix)[1]
	reverse_rownums = seq(numrows,1,-1)
	downpass_edgematrix = uppass_edgematrix[reverse_rownums,:]
	return(downpass_edgematrix)
end



# Preorder traversal (I think)
"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
#great_ape_newick_string = "(((human:6,(chimp1:0.5,bonobo:0.5):5.5):1,gorilla:7):5,orangutan:12);"
#great_ape_newick_string = "(((human:6,(chimp1:0.5,bonobo:0.5):5.5):1,gorilla:7):5,(orang1:1,orang2:1):11);"
tr = readTopology(great_ape_newick_string)

# Preorder traversal
uppass_nodeIndexes = get_LR_uppass_nodeIndexes(tr)

# Reverse preorder traversal (not the same as postorder!)
downpass_nodeIndexes = get_LR_downpass_nodeIndexes(tr)
"""
function get_LR_uppass_nodeIndexes(tr)
	rootnodenum = tr.root
	iterNum = 1
	nodeIndex_array = collect(repeat([0], tr.numNodes))
	indexNum_table = get_nodeIndex_PNnumber(tr)

	res = get_nodenumbers_above_node(tr, rootnodenum, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
	uppass_nodeIndexes = res[1]
	return(uppass_nodeIndexes)
end

# Reverse-Preorder traversal (I think)
"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
#great_ape_newick_string = "(((human:6,(chimp1:0.5,bonobo:0.5):5.5):1,gorilla:7):5,orangutan:12);"
#great_ape_newick_string = "(((human:6,(chimp1:0.5,bonobo:0.5):5.5):1,gorilla:7):5,(orang1:1,orang2:1):11);"
tr = readTopology(great_ape_newick_string)

# Preorder traversal
uppass_nodeIndexes = get_LR_uppass_nodeIndexes(tr)

# Reverse preorder traversal (not the same as postorder!)
downpass_nodeIndexes = get_LR_downpass_nodeIndexes(tr)
"""
function get_LR_downpass_nodeIndexes(tr)
	rootnodenum = tr.root
	iterNum = 1
	nodeIndex_array = collect(repeat([0], tr.numNodes))
	indexNum_table = get_nodeIndex_PNnumber(tr)

	res = get_nodenumbers_above_node(tr, rootnodenum, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
	downpass_nodeIndexes = reverse(res[1])
	return(downpass_nodeIndexes)
end


# Get Rnodenums
# They seem to be just 
# sort(tip_PNnumbers) + sort(abs(internal_PNnumbers))

"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root
trdf = prt(tr, rootnodenum)
df

# Get the edges, get the nodeIndexes corresponding to these
edge = tr.edge
edge

one_edge = edge[1]

edge_df = get_NodeIndex_df_by_tree_edges(tr, indexNum_table=indexNum_table)
edge_df

uppass_edgematrix = get_LR_uppass_edgematrix(tr)
uppass_edgematrix

downpass_edgematrix = get_LR_downpass_edgematrix(tr)
downpass_edgematrix


# Get the R node numbers, append to indexNum_table
Rnodenums_in_indexNum_table_order = get_Rnodenums(tr, indexNum_table)
indexNum_table2 = hcat(indexNum_table, Rnodenums_in_indexNum_table_order)
indexNum_table2

indexNum_table3 = indexNum_table2[sortperm(indexNum_table2[:,3]),:]
indexNum_table3

"""


function get_Rnodenums(tr, indexNum_table)
	numnodes = length(tr.node)
	Rnodenums = collect(1:numnodes)
	tipsTF = indexNum_table[:,2] .> 0
	internalTF = tipsTF .== false
	
	numtips = sum(tipsTF)
	numinternal = sum(internalTF)
	
	
	PNnumbers_in_Rnodenums_order = collect(repeat([0], numnodes))
	PNnumbers_in_Rnodenums_order[1:numtips] = sort(indexNum_table[:,2][tipsTF])
	
	PNnumbers_in_Rnodenums_order[(numtips+1):(numtips+numinternal)] = reverse(sort(indexNum_table[:,2][internalTF]))
	
	
	tmpmat = hcat(Rnodenums, PNnumbers_in_Rnodenums_order)
	
	# Match to indexNum_table
	indices = collect(1:numnodes)
	Rnodenums_in_indexNum_table_order = collect(repeat([0], numnodes))
	for i in 1:numnodes
		current_PNnumber = indexNum_table[i,2] 
		TF = tmpmat[:,2] .== current_PNnumber
		rownum_in_tmpmat = indices[TF][1]
		Rnodenums_in_indexNum_table_order[i] = tmpmat[rownum_in_tmpmat,1]
	end
	return(Rnodenums_in_indexNum_table_order)
end





"""
using DataFrames
using PhyloNetworks
using PhyloPlots
include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr

indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

"""
function get_nodeIndex_PNnumber(tr)
	# Get a numNode x 2 table
	# Index, then ".number"
	numnodes = length(tr.node)
	indexNum_table = Array{Int}(undef, numnodes, 2)
	for i in 1:numnodes
		indexNum_table[i,1] = i
		indexNum_table[i,2] = tr.node[i].number
	end
	return(indexNum_table)
end



# Go from a PhyloNetwork Node Number (PNnumber) to the node index 
# (i.e., the index of that node in the list of nodes)
function get_nodeIndex_from_PNnumber(PNnumber; indexNum_table)
	TF01 = indexNum_table[:,2] .== PNnumber
	# Adding the [1] returns a scalar
	nodeIndex = indexNum_table[:,1][TF01][1]
	return nodeIndex
end



"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root
trdf = prt(tr, rootnodenum)
trdf
"""
# Return a DataFrame with the edge numbers
function prt(tr, rootnodenum, get_taxa_by_node=true)
	#using DataFrames
	numnodes = length(tr.node)
	# number of digits for internal node numbers
	numdigits = length(digits(numnodes))
	if (numdigits < 2)
		numdigits = 2
	end
	
	# Initialize the dataframe
	trdf = DataFrames.DataFrame(nodeIndex=collect(1:numnodes), PNnumber=collect(repeat([0],numnodes)))
	
	# Fill in the PNnumber node numbers
	indexNum_table = get_nodeIndex_PNnumber(tr)
	trdf[!, :PNnumber] = indexNum_table[:,2]
	trdf
	
	# Add the R node numbers
	Rnodenums = get_Rnodenums(tr, indexNum_table)
	trdf[!,:Rnodenums] = Rnodenums
	trdf
	
	
	# Add the node ages
	node_age = get_node_ages(tr)
	trdf[!,:node_age] = node_age
	
	# Add branch lengths
	brlen = collect(repeat([0.0], numnodes))
	ancNodeIndex = collect(repeat([0], numnodes))
	leftNodeIndex = collect(repeat([0], numnodes))
	rightNodeIndex = collect(repeat([0], numnodes))
	nodeName = collect(repeat([""], numnodes))
	nodeType = collect(repeat([""], numnodes))
	
	edge_df = get_NodeIndex_df_by_tree_edges(tr, indexNum_table=indexNum_table)
	edge_df_rownums = collect(1:Rnrow(edge_df))
	for i in 1:Rnrow(trdf)
		nodeIndex = trdf[i,:nodeIndex]
		TF = edge_df[:,:edge_decNodeIndex] .== nodeIndex
		if (sum(TF) > 0)
			edge_df_rownum = edge_df_rownums[TF][1]
			brlen[i] = edge_df[:,:edge_length][edge_df_rownum]
			ancNodeIndex[i] = edge_df[:,:edge_ancNodeIndex][edge_df_rownum]
		else
			brlen[i] = -999.0
			ancNodeIndex[i] = -999
		end
		
		# Left and right descendant nodeIndexes
		TF_edgerows_descend_from_decNodeIndex = edge_df[:,:edge_ancNodeIndex] .== nodeIndex
		if (sum(TF_edgerows_descend_from_decNodeIndex) > 0)
			# Internal node
			leftNodeIndex[i] = edge_df[TF_edgerows_descend_from_decNodeIndex,:edge_decNodeIndex][1]
			rightNodeIndex[i] = edge_df[TF_edgerows_descend_from_decNodeIndex,:edge_decNodeIndex][2]
			nodeType[i] = "internal"
			if (nodeIndex == tr.root)
				nodeType[i] = "root"
			end
			
			# Node names
			tmp_nodeName = tr.node[nodeIndex].name
			if (tmp_nodeName == "")
				internal_nodeIndex_as_string = lpad(string(nodeIndex), numdigits, "0")
				nodeName[i] = join(["in", internal_nodeIndex_as_string], "")
			else
				nodeName[i] = tmp_nodeName
			end
		else
			# Tip node
			leftNodeIndex[i] = -999
			rightNodeIndex[i] = -999
			nodeType[i] = "tip"
			nodeName[i] = tr.node[nodeIndex].name
		end
	end # END for i in 1:Rnrow(trdf)
	
	# Add the fields
	trdf[!,:brlen] = brlen
	trdf[!,:ancNodeIndex] = ancNodeIndex
	trdf[!,:leftNodeIndex] = leftNodeIndex
	trdf[!,:rightNodeIndex] = rightNodeIndex
	trdf[!,:nodeName] = nodeName
	trdf[!,:nodeType] = nodeType
	
	if get_taxa_by_node == true
		downpass_edgematrix = get_LR_downpass_edgematrix(tr)
		taxa = get_taxa_descending_from_each_node(tr, trdf, downpass_edgematrix=downpass_edgematrix)
		trdf[!,:taxa] = taxa
	end
	
	return(trdf)
end



# Get tipnames descending from each node
# The list is comma-delimited and alphabetical
"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root

# Get trdf WITHOUT "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, false)

# Not all columns print now that trdf is big, let's print just the left and right
# Get trdf WITH "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, true);
headLR(trdf)

downpass_edgematrix = get_LR_downpass_edgematrix(tr)
taxa = get_taxa_descending_from_each_node(tr, trdf, downpass_edgematrix=get_LR_downpass_edgematrix(tr))
"""

function get_taxa_descending_from_each_node(tr, trdf; downpass_edgematrix=get_LR_downpass_edgematrix(tr))
	# Get sizes
	numnodes = length(tr.node)
	numIndices = length(sort(unique(flat2(downpass_edgematrix))))
	
	# Error check
	if (numnodes != numIndices)
		txt = ["ERROR in get_taxa_descending_from_each_node(): the number of nodes in tr (", string(numnodes), ") does not match the number of rows in downpass_edgematrix (", string(numIndices), ")."]
		error(join(txt, ""))
	end
	
	taxa = collect(repeat([""], numnodes))
	
	# Assumes a binary tree, and a nodeIndexes df in downpass order, from
	# downpass_edgematrix = get_LR_downpass_edgematrix(tr)
	
	# Step through the downpass_edgematrix, in pairs
	edgematrix_rows_to_visit = collect(1:2:Rnrow(downpass_edgematrix))
	for (iter,i) in enumerate(edgematrix_rows_to_visit)
		j = i+1
		nodeIndex_left = downpass_edgematrix[i,2]
		nodeIndex_right = downpass_edgematrix[j,2]
		nodeIndex_ancL = downpass_edgematrix[i,1]
		nodeIndex_ancR = downpass_edgematrix[j,1]
		if (nodeIndex_ancL == nodeIndex_ancR)
			nodeIndex_anc = nodeIndex_ancL
		else
			txt = ["Error in get_taxa_descending_from_each_node(): in downpass_edgematrix, at i=", string(i), ", j=", string(j), ", the nodeIndex_anc's don't match! Printing these rows of downpass_edgematrix."]
			msg = paste0(txt)
			println(msg)
			print(downpass_edgematrix[i:j,:])
			error(msg)
		end
		
		tmp_taxa = collect(repeat([""], 2))
		if (tr.node[nodeIndex_left].leaf == true)
			taxa[nodeIndex_left] = tr.node[nodeIndex_left].name
			tmp_taxa[1] = tr.node[nodeIndex_left].name
		else
			tmp_taxa[1] = taxa[nodeIndex_left]
		end
		if (tr.node[nodeIndex_right].leaf == true)
			taxa[nodeIndex_right] = tr.node[nodeIndex_right].name
			tmp_taxa[2] = tr.node[nodeIndex_right].name
		else
			tmp_taxa[2] = taxa[nodeIndex_right]
		end
		
		# Join the two lists and put in the ancestral node
		ancNodeIndex_of_Left = nodeIndex_ancL
		ancNodeIndex_of_Right = nodeIndex_ancR
		if (ancNodeIndex_of_Left != ancNodeIndex_of_Right)
			txt = ["ERROR in get_taxa_descending_from_each_node(): ancNodeIndex_of_Left must match ancNodeIndex_of_Right in trdf, but doesn't.\nAncestor of nodeIndex_left ", string(nodeIndex_left), " is ", string(ancNodeIndex_of_Left), ", but\nancestor of nodeIndex_right ", string(nodeIndex_right), " is ", string(ancNodeIndex_of_Right), "\n"]
			msg = join(txt, "")
			println(msg)
			error(msg)
		end
		# No error, so ancNodeIndex was found:
		ancNodeIndex = ancNodeIndex_of_Left
		taxa_unordered = flat2(split.(tmp_taxa, ","))
		taxa_ordered = sort(taxa_unordered)
		taxa_at_ancNode = join(taxa_ordered, ",")
		taxa[ancNodeIndex] = taxa_at_ancNode
	end
	
	return(taxa)
end









# Is the node a tip?
function isTip_TF(tr, nodeIndex)
	# Returns true if it's a tip (a "leaf"), false if not
	return(tr.node[nodeIndex].leaf)
end

# Get the nodeIndex ancestral to 2 nodes



# Elaborate function to MAKE SURE we are getting the ancestor & descendant PNnumbers correct

"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root
trdf = prt(tr, rootnodenum)
df

# Get the edges, get the nodeIndexes corresponding to these
edge = tr.edge
edge

one_edge = edge[1]

# Every edge has only 2 nodes. 
for i in 1:length(edge)
	one_edge = edge[i]
	anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
	print(anc_decPNnumbers)
end
"""

function get_NodeIndexes_from_edge(one_edge)
	# Make sure there are only 2 nodes for each edge.
	numnodec_on_this_edge = length(one_edge.node)
	if ( numnodec_on_this_edge != 2 )
		array_of_strings = ["\nError in get_NodeIndexes_from_edge(): this function assumes each edge has only 2 nodes, ancestor & descendant.\nInstead, this edge gives length(edge.node)=", string(numnodec_on_this_edge), ".\nPrinting edge to screen...\n"]
		txt = join(array_of_strings, "")
		println(txt)
		print(edge)
		error(txt)
	end
	
	# Initialize
	tmp_nodeIndices = collect(repeat([0], numnodec_on_this_edge))
	tmp_nodePNnumbers = collect(repeat([0], numnodec_on_this_edge))
	ancPNnumber = 0
	decPNnumber = 0
	
	# Record the PN (PhyloNetworks) node numbers
	tmp_nodePNnumbers[1] = one_edge.node[1].number
	tmp_nodePNnumbers[2] = one_edge.node[2].number
	tmp_nodeIndices_on_edge = collect(1:2)
	
	# Positive node numbers are tips; negative node numbers are internal

	# Declare error if both node numbers are the same
	PNnumber_eq_0_TF = tmp_nodePNnumbers .== 0
	if (sum(PNnumber_eq_0_TF) > 0)
		error("Error in get_NodeIndexes_from_edge(): both node numbers attached to this edge are the same. This should be impossible.")
	end

	
	# Declare error if both node numbers are positive (both are tips)
	positiveTF = tmp_nodePNnumbers .> 0
	if (sum(positiveTF) == 2)
		error("Error in get_NodeIndexes_from_edge(): both node numbers attached to this edge are tips (PhyloNetworks node number > 0). This should be impossible.")
	end
	
	# If one PNnumber is positive and one negative, then you have a descendant tip, and ancestor internal node
	if (sum(positiveTF) == 1)
		anc_nodeIndex_on_edge = tmp_nodeIndices_on_edge[positiveTF .== false]
		ancPNnumber = tmp_nodePNnumbers[anc_nodeIndex_on_edge] # internal nodenum
		dec_nodeIndex_on_edge = tmp_nodeIndices_on_edge[positiveTF .== true]
		decPNnumber = tmp_nodePNnumbers[dec_nodeIndex_on_edge]
	end
	
	# If both are negative, then you have 2 internal nodes. The one closest to 0 is the root
	if (sum(positiveTF) == 0)
		min_absval = -1 * minimum(abs.(tmp_nodePNnumbers))
		matches_min_TF = tmp_nodePNnumbers .== min_absval
		anc_nodeIndex_on_edge = tmp_nodeIndices_on_edge[matches_min_TF .== true]
		dec_nodeIndex_on_edge = tmp_nodeIndices_on_edge[matches_min_TF .== false]
		ancPNnumber = tmp_nodePNnumbers[anc_nodeIndex_on_edge] # ancestral internal nodenum
		decPNnumber = tmp_nodePNnumbers[dec_nodeIndex_on_edge] # ancestral internal nodenum
	end
	
	anc_decPNnumbers = [ancPNnumber[1], decPNnumber[1]]
	return(anc_decPNnumbers)
end


"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

edge_df = get_NodeIndex_df_by_tree_edges(tr, indexNum_table=indexNum_table)
"""
function get_NodeIndex_df_by_tree_edges(tr; indexNum_table)
	# Get the ancestral and descendant node numbers for each edge
	num_edges = length(tr.edge)
	edgeIndex = collect(1:num_edges)
	edge_ancNodeIndex = collect(repeat([0], num_edges))
	edge_decNodeIndex = collect(repeat([0], num_edges))
	edge_ancPNnumber = collect(repeat([0], num_edges))
	edge_decPNnumber = collect(repeat([0], num_edges))
	edge_length = collect(repeat([0.0], num_edges))

	# Initialize the dataframe
	edge_df = DataFrames.DataFrame(edgeIndex=edgeIndex, edge_ancNodeIndex=edge_ancNodeIndex, edge_decNodeIndex=edge_decNodeIndex, edge_ancPNnumber=edge_ancPNnumber, edge_decPNnumber=edge_decPNnumber, edge_length=edge_length)
	edge_df

	# Every edge has only 2 nodes. 
	for i in 1:length(tr.edge)
		one_edge = tr.edge[i]
		anc_dec_PNnumbers = get_NodeIndexes_from_edge(one_edge)
	
		# Use a named keyword argument, indexNum_table, to provide the translation from PNnumber to node index
		anc_dec_nodeIndices = get_nodeIndex_from_PNnumber.(anc_dec_PNnumbers, indexNum_table=indexNum_table)
		edge_df[i, :edge_ancNodeIndex] = anc_dec_nodeIndices[1]
		edge_df[i, :edge_decNodeIndex] = anc_dec_nodeIndices[2]
		edge_df[i, :edge_ancPNnumber] = anc_dec_PNnumbers[1]
		edge_df[i, :edge_decPNnumber] = anc_dec_PNnumbers[2]
		edge_df[i, :edge_length] = one_edge.length
	end

	return(edge_df)
end


"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

# Make a trdf, a DataFrame holding all of the key tree information
cumulative_height_at_each_node = get_node_heights(tr)
node_age = get_node_ages(tr)
"""
function get_node_heights(tr)
	indexNum_table = get_nodeIndex_PNnumber(tr)
	uppass_nodeIndexes = get_LR_uppass_nodeIndexes(tr)
	edge_df = get_NodeIndex_df_by_tree_edges(tr, indexNum_table=indexNum_table)
	cumulative_height_at_each_node = collect(repeat([0.0], length(tr.node)))
	
	# Iterate up through the nodes from the root
	for i in 1:length(uppass_nodeIndexes)
		nodeIndex = uppass_nodeIndexes[i]
		if (nodeIndex == tr.root)
			cumulative_height_at_each_node[nodeIndex] = 0.0
		else
			ancTF = edge_df[:,:edge_decNodeIndex] .== nodeIndex
			anc_nodeIndex = edge_df[:,:edge_ancNodeIndex][ancTF][1]
			# Get previous age
			previous_height_above_root = cumulative_height_at_each_node[anc_nodeIndex]
			length_of_this_branch = edge_df[:,:edge_length][ancTF][1]
			current_height = length_of_this_branch + previous_height_above_root
			cumulative_height_at_each_node[nodeIndex] = current_height
		end
	end
	return(cumulative_height_at_each_node)
end




"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

# Make a trdf, a DataFrame holding all of the key tree information
cumulative_height_at_each_node = get_node_heights(tr)
node_age = get_node_ages(tr)
"""
function get_node_ages(tr)
	cumulative_height_at_each_node = get_node_heights(tr)
	tree_height = maximum(cumulative_height_at_each_node)
	node_age = tree_height .- cumulative_height_at_each_node
	return node_age
end
















#######################################################
# Parallel operations on binary trees
#######################################################
# Threaded downpass that spawns new processes when the 2 nodes above are done.


# Results structure
struct Res
	# The states can be 
	# "not_ready" (value at branch top NOT available)
	# "ready_for_nodeOp" (values at branches above are ready)
	# "ready_for_branchOp" (value at branch top available)
	# "calculating" (value at branch bottom being calculated)
	# "done" (value at branch bottom available)
	node_state::Array{String}
	node_Lparent_state::Array{String}
	node_Rparent_state::Array{String}
	
	# Tree structure
	root_nodeIndex::Int64
	numNodes::Int64
	uppass_edgematrix::Array{Int64}
	likes_at_each_nodeIndex_branchTop::Array{Array{Float64,1},1}
	likes_at_each_nodeIndex_branchBot::Array{Array{Float64,1},1}
	thread_for_each_nodeOp::Array{Int64}
	thread_for_each_branchOp::Array{Int64}
	
	# Calculation stats
	calc_spawn_start::Array{DateTime}
	calc_start_time::Array{DateTime}
	calc_end_time::Array{DateTime}
	calc_duration::Array{Float64}
	calctime_iterations::Array{Float64}
end

# Construct a default, simple results structure
# (likes_at_each_nodeIndex_branchTop is an array)
function construct_Res_old()
	n = 1 # number of states
	node_state = ["ready_for_branchOp", "ready_for_branchOp", "not_ready", "ready_for_branchOp", "not_ready", "ready_for_branchOp", "not_ready"]
	node_Lparent_state = ["NA", "NA", "not_ready", "NA", "not_ready", "NA", "not_ready"]
	node_Rparent_state = ["NA", "NA", "not_ready", "NA", "not_ready", "NA", "not_ready"]
	root_nodeIndex = 7
	numNodes = 7
	uppass_edgematrix = [7 6; 7 5; 5 4; 5 3; 3 2; 3 1]
	likes_at_each_nodeIndex_branchTop = collect(repeat([1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0], n))
	likes_at_each_nodeIndex_branchBot = collect(repeat([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], n))
	thread_for_each_nodeOp = collect(repeat([0], 7))
	thread_for_each_branchOp = collect(repeat([0], 7))

	calc_spawn_start = collect(repeat([Dates.now()], 7))
	calc_start_time = collect(repeat([Dates.now()], 7))
	calc_end_time = collect(repeat([Dates.now()], 7))
	calc_duration = collect(repeat([0.0], 7))

	calctime_iterations = [0.0, 0.0]

	res = Res(node_state, node_Lparent_state, node_Rparent_state, root_nodeIndex, numNodes, uppass_edgematrix, likes_at_each_nodeIndex_branchTop, likes_at_each_nodeIndex_branchBot, thread_for_each_nodeOp, thread_for_each_branchOp, calc_spawn_start, calc_start_time, calc_end_time, calc_duration, calctime_iterations)
	return res
end


# Construct a default, simple results structure
# (likes_at_each_nodeIndex_branchTop is an array of arrays)
function construct_Res()
	n = 1 # number of states
	numNodes = 7  # number of nodes
	node_state = ["ready_for_branchOp", "ready_for_branchOp", "not_ready", "ready_for_branchOp", "not_ready", "ready_for_branchOp", "not_ready"]
	node_Lparent_state = ["NA", "NA", "not_ready", "NA", "not_ready", "NA", "not_ready"]
	node_Rparent_state = ["NA", "NA", "not_ready", "NA", "not_ready", "NA", "not_ready"]
	root_nodeIndex = 7
	uppass_edgematrix = [7 6; 7 5; 5 4; 5 3; 3 2; 3 1]
	likes_OneNode = collect(repeat([0.0], n))
	likes_at_each_nodeIndex_branchTop = repeat([likes_OneNode], numNodes)
	likes_at_each_nodeIndex_branchBot = repeat([likes_OneNode], numNodes)
	
	default_likes_at_each_nodeIndex_branchTop = [1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0]
	for i in 1:length(likes_at_each_nodeIndex_branchTop)
		likes_at_each_nodeIndex_branchTop[i][1] = default_likes_at_each_nodeIndex_branchTop[i]
	end
	
	typeof(likes_at_each_nodeIndex_branchTop)
	thread_for_each_nodeOp = collect(repeat([0], 7))
	thread_for_each_branchOp = collect(repeat([0], 7))

	calc_spawn_start = collect(repeat([Dates.now()], 7))
	calc_start_time = collect(repeat([Dates.now()], 7))
	calc_end_time = collect(repeat([Dates.now()], 7))
	calc_duration = collect(repeat([0.0], 7))

	calctime_iterations = [0.0, 0.0]

	res = Res(node_state, node_Lparent_state, node_Rparent_state, root_nodeIndex, numNodes, uppass_edgematrix, likes_at_each_nodeIndex_branchTop, likes_at_each_nodeIndex_branchBot, thread_for_each_nodeOp, thread_for_each_branchOp, calc_spawn_start, calc_start_time, calc_end_time, calc_duration, calctime_iterations)
	return res
end





"""
# Load a simple tree, see a simple list of nodeIndexes
using DataFrames
using PhyloNetworks
include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr.root
# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)

"""
function construct_Res(tr::HybridNetwork)
	root_nodeIndex = tr.root
	numNodes = tr.numNodes
	uppass_edgematrix = get_LR_uppass_edgematrix(tr)
	
	# Give tip nodeIndexes their nodeNodex as the "likelihood"
	indexNum_table = get_nodeIndex_PNnumber(tr)
	tipsTF = indexNum_table[:,2] .> 0
	tipLikes = indexNum_table[tipsTF,2] * 1.0
	likes_at_each_nodeIndex_branchTop = collect(repeat([0.0], numNodes))
	likes_at_each_nodeIndex_branchTop[tipsTF] = tipLikes
	
	# Fill in the node_states
	node_state = collect(repeat(["not_ready"], numNodes))
	node_state[tipsTF] .= "ready_for_branchOp"
	node_Lparent_state = collect(repeat(["not_ready"], numNodes))
	node_Rparent_state = collect(repeat(["not_ready"], numNodes))
	node_Lparent_state[tipsTF] .= "NA"
	node_Rparent_state[tipsTF] .= "NA"
	
	# Initialize with zeros for the other items
	likes_at_each_nodeIndex_branchBot = collect(repeat([0.0], numNodes))
	thread_for_each_nodeOp = collect(repeat([0], numNodes))
	thread_for_each_branchOp = collect(repeat([0], numNodes))
	
	calc_spawn_start = collect(repeat([Dates.now()], numNodes))
	calc_start_time = collect(repeat([Dates.now()], numNodes))
	calc_end_time = collect(repeat([Dates.now()], numNodes))
	calc_duration = collect(repeat([0.0], numNodes))

	calctime_iterations = [0.0, 0.0]
	number_of_whileLoop_iterations = [0]	

	# Initialize res object
	res = Res(node_state, node_Lparent_state, node_Rparent_state, root_nodeIndex, numNodes, uppass_edgematrix, likes_at_each_nodeIndex_branchTop, likes_at_each_nodeIndex_branchBot, thread_for_each_nodeOp, thread_for_each_branchOp, calc_spawn_start, calc_start_time, calc_end_time, calc_duration, calctime_iterations)
	return res
end





"""
# Load a simple tree, see a simple list of nodeIndexes
using DataFrames
using PhyloNetworks
include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr.root
# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
n=10 # number of states in the state space

"""
function construct_Res(tr::HybridNetwork, n)
	root_nodeIndex = tr.root
	numNodes = tr.numNodes
	uppass_edgematrix = get_LR_uppass_edgematrix(tr)
	
	
	# Set up an array of length nstates (n), to hold the likelihoods for each node
	blank_states = collect(repeat([0.0], n))
	likes_at_each_nodeIndex_branchTop = collect(repeat([blank_states], numNodes))

	# Give tip nodeIndexes a likelihood of 1 at all states
	indexNum_table = get_nodeIndex_PNnumber(tr)
	tipsTF = indexNum_table[:,2] .> 0
	tipnums = seq(1, length(tipsTF), 1)[tipsTF]
	
	for i in 1:length(tipnums)
		tipLikes = collect(repeat([1.0], n))
		likes_at_each_nodeIndex_branchTop[tipnums[i]] = tipLikes
	end

	
	# Fill in the node_states
	node_state = collect(repeat(["not_ready"], numNodes))
	node_state[tipsTF] .= "ready_for_branchOp"
	node_Lparent_state = collect(repeat(["not_ready"], numNodes))
	node_Rparent_state = collect(repeat(["not_ready"], numNodes))
	node_Lparent_state[tipsTF] .= "NA"
	node_Rparent_state[tipsTF] .= "NA"
	
	# Initialize with zeros for the other items
	likes_at_each_nodeIndex_branchBot = collect(repeat([blank_states], numNodes))
	thread_for_each_nodeOp = collect(repeat([0], numNodes))
	thread_for_each_branchOp = collect(repeat([0], numNodes))
	
	calc_spawn_start = collect(repeat([Dates.now()], numNodes))
	calc_start_time = collect(repeat([Dates.now()], numNodes))
	calc_end_time = collect(repeat([Dates.now()], numNodes))
	calc_duration = collect(repeat([0.0], numNodes))

	calctime_iterations = [0.0, 0.0]
	number_of_whileLoop_iterations = [0]	

	# Initialize res object
	res = Res(node_state, node_Lparent_state, node_Rparent_state, root_nodeIndex, numNodes, uppass_edgematrix, likes_at_each_nodeIndex_branchTop, likes_at_each_nodeIndex_branchBot, thread_for_each_nodeOp, thread_for_each_branchOp, calc_spawn_start, calc_start_time, calc_end_time, calc_duration, calctime_iterations)
	return res
end




# Convert this res object to a DataFrame
#function res_to_df(res)
#	
#end

function count_nodes_finished(node_state)
	sum(node_state .== "done")
end


# Combine likelihoods from above
function nodeOp(current_nodeIndex, res)
	res.node_state[current_nodeIndex] = "calculating_nodeOp"
	uppass_edgematrix = res.uppass_edgematrix
	
	# Record the thread, for kicks
	tmp_threadID = Threads.threadid()
	res.thread_for_each_nodeOp[current_nodeIndex] = tmp_threadID
	TF = uppass_edgematrix[:,1] .== current_nodeIndex
	if (sum(TF) == 2)
		# Get likelihoods from above (iterates up to tips)
		parent_nodeIndexes = uppass_edgematrix[TF,2]
		tmp1 = res.likes_at_each_nodeIndex_branchTop[parent_nodeIndexes[1]]
		tmp2 = res.likes_at_each_nodeIndex_branchTop[parent_nodeIndexes[2]]

		# Check that data are actually available
		if (sum(tmp1) == 0.0)
			txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): sum(tmp1) == 0.0, indicating data at parent nodes actually not available."], "")
			res.node_state[current_nodeIndex] = txt
			print("\n")
			print(txt)
			print("\n")
			return(error(txt))
		end

		if (sum(tmp2) == 0.0)
			txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): sum(tmp2) == 0.0, indicating data at parent nodes actually not available."], "")
			res.node_state[current_nodeIndex] = txt
			print("\n")
			print(txt)
			print("\n")
			return(error(txt))
		end

		nodeData_at_top = tmp1 + tmp2
		res.likes_at_each_nodeIndex_branchTop[current_nodeIndex] = nodeData_at_top
		
		# Check if it's the root node
		if (current_nodeIndex == res.root_nodeIndex)
			res.node_state[current_nodeIndex] = "done"
		else
			res.node_state[current_nodeIndex] = "ready_for_branchOp"
		end
		return()
	elseif (sum(TF) == 0)
	  # If a tip
	  txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): shouldn't be run on a tip node."], "")
	  print("\n")
	  print(txt)
	  print("\n")
		return(error(txt))
	else
	  txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): sum(TF) should be 0 or 2"], "")
	  print("\n")
	  print(txt)
	  print("\n")
		return(error(txt))
	end
	txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): shouldn't get here."], "")
	print("\n")
	print(txt)
	print("\n")
	return(error(txt))
end

# Calculate down a branch
# This function can read from res, but writing to res is VERY BAD as 
# it created conflicts apparently when there were more @spawns than cores
# Do all the writing to res in the while() loop
function branchOp_example(current_nodeIndex, res; num_iterations=10000000)
	calc_start_time = Dates.now()
	spawned_nodeIndex = current_nodeIndex
	tmp_threadID = Threads.threadid()
	
	# Example slow operation
	y = countloop(num_iterations, current_nodeIndex)

	nodeData_at_top = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
	nodeData_at_bottom = nodeData_at_top / 2.0

	return(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time)
end




# Calculate Ds down a branch
#
# Modifies branchOp to do Ds calculation down a branch
#
# This function can read from res, but writing to res is VERY BAD as 
# it created conflicts apparently when there were more @spawns than cores
# Do all the writing to res in the while() loop
function branchOp_ClaSSE_Ds_v5(current_nodeIndex, res; u0, tspan, p_Ds_v5)
	calc_start_time = Dates.now()
	spawned_nodeIndex = current_nodeIndex
	tmp_threadID = Threads.threadid()
	
	# Example slow operation
	#y = countloop(num_iterations, current_nodeIndex)
	prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p_Ds_v5)
	sol_Ds = solve(prob_Ds_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)

	nodeData_at_top = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
	#nodeData_at_bottom = nodeData_at_top / 2.0
	nodeData_at_bottom = sol_Ds
	
	return(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time)
end



"""
	Set up inputs
"""

function setup_inputs_branchOp_ClaSSE_Ds_v5(u0, tspan, p_Ds_v5; solver="Tsit5()", 
				 save_everystep="false", abstol="1e-9", reltol="1e-9")
	
	prob_str = "prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p_Ds_v5)"
	solve_str = join(["sol_Ds = solve(prob_Ds_v5, ", solver, ", save_everystep=", save_everystep, ", abstol=", abstol, ", reltol=", reltol, ")"])
	store_str = "nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)]"
	
	# Assemble the NamedTuple of inputs
	inputs = (u0=u0, tspan=tspan, p_Ds_v5=p_Ds_v5, solver=solver, save_everystep=save_everystep, abstol=abstol, reltol=reltol, prob_str=prob_str, solve_str=solve_str, store_str=store_str)
	return inputs
end


# Calculate Ds down a branch
#
# Modifies branchOp to do Ds calculation down a branch
#
# This function can read from res, but writing to res is VERY BAD as 
# it created conflicts apparently when there were more @spawns than cores
# Do all the writing to res in the while() loop
function branchOp(current_nodeIndex, res, inputs)
	calc_start_time = Dates.now()
	spawned_nodeIndex = current_nodeIndex
	tmp_threadID = Threads.threadid()
	
	
	# The old practice input was an Int64
	if (typeof(inputs) != Int64)
		u0 = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
		tspan = inputs.tspan
		p_Ds_v5 = inputs.p_Ds_v5
		#solver = inputs.solver
		#save_everystep = inputs.save_everystep
		#abstol = inputs.abstol
		#reltol = inputs.reltol


		# Example slow operation
		#y = countloop(num_iterations, current_nodeIndex)
		#	prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p_Ds_v5)
		#	sol_Ds = solve(prob_Ds_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)

		eval(Meta.parse(inputs.prob_str))
		eval(Meta.parse(inputs.solve_str))
		eval(Meta.parse(inputs.store_str))
	else
		nodeData_at_top = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
		nodeData_at_bottom = nodeData_at_top / 2.0
		#nodeData_at_bottom = sol_Ds
	end

	
	return(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time)
end




function countloop(num_iterations, current_nodeIndex)
	x = 0.0
	random_number_generator = MersenneTwister(current_nodeIndex);

	for i in 1:num_iterations
	   x = x + (randn(random_number_generator, 1)[1] / num_iterations)
	end
	return(x)
end



"""
Iterate through the "res" object many times to complete the downpass, spawning jobs along the way
"""
function iterative_downpass!(res; max_iterations=10^10, num_iterations=10000000)

	# Check number of threads
	numthreads = Threads.nthreads()
	parallel_TF = numthreads > 1
	if (parallel_TF == false)
		txt = "Error in iterative_downpass!(): This function probably requires multiple threads operating to consistently compile. Try starting julia with e.g. 'JULIA_NUM_THREADS=8 julia'."
		error(txt)
	end

	diagnostics = collect(repeat([Dates.now()], 3))
	diagnostics[1] = Dates.now()
	
	# Setup
	current_nodeIndex = res.root_nodeIndex
	tasks = Any[]
	tasks_fetched_TF = Bool[]
	are_we_done = false

	iteration_number = 0
	while(are_we_done == false)
		iteration_number = iteration_number+1
		# As long as all the nodes are not done,
		# check for "ready" nodes
		# When they finish, change to "done"
		indexes_ready = findall(res.node_state .== "ready_for_branchOp")
		for current_nodeIndex in indexes_ready
			# Before spawning, do some checks
			res.node_state[current_nodeIndex] = "calculating_branchOp"
			# Check for root; no calculation on root branch for now
			if current_nodeIndex == res.root_nodeIndex
				res.node_state[current_nodeIndex] = "done"
				return()
			end

			# Spawn a branch operation, and a true-false of whether they are fetched
			res.calc_spawn_start[current_nodeIndex] = Dates.now()
			push!(tasks, @spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
			push!(tasks_fetched_TF, false)
		end
	
		# Check which jobs are done, fetch them, and update status of that node
		num_tasks = length(tasks)
		for i in 1:num_tasks
			if (tasks_fetched_TF[i] == false)
				if (istaskdone(tasks[i]) == true)
					# Get the results
					calc_end_time = Dates.now()
					(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time) = fetch(tasks[i])
					
					# Store run information
					res.calc_start_time[spawned_nodeIndex] = calc_start_time
					res.calc_end_time[spawned_nodeIndex] = calc_end_time
					res.calc_duration[spawned_nodeIndex] = (calc_end_time - calc_start_time).value / 1000.0
					tasks_fetched_TF[i] = true
					
					# Record information
					res.thread_for_each_branchOp[spawned_nodeIndex] = tmp_threadID
					res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = nodeData_at_bottom
					# Get the ancestor nodeIndex
					uppass_edgematrix = res.uppass_edgematrix
					TF = uppass_edgematrix[:,2] .== spawned_nodeIndex
					parent_nodeIndex = uppass_edgematrix[TF,1][1]

					# Get the left daughter nodeIndex (1st in the uppass_edgematrix)
					edge_rows_TF = uppass_edgematrix[:,1] .== parent_nodeIndex
					left_nodeIndex = uppass_edgematrix[edge_rows_TF,2][1]
					right_nodeIndex = uppass_edgematrix[edge_rows_TF,2][2]

					# Update the state of the parent_node's daughters
					if (spawned_nodeIndex == left_nodeIndex)
						res.node_Lparent_state[parent_nodeIndex] = "ready"
					end
					if (spawned_nodeIndex == right_nodeIndex)
						res.node_Rparent_state[parent_nodeIndex] = "ready"
					end

					# Update the state of the current node
					res.node_state[spawned_nodeIndex] = "done"
				end
			end
		end
	
		# Update which nodes have had both parents complete
		TF1 = res.node_state .== "not_ready"
		TF2 = res.node_Lparent_state .== "ready"
		TF3 = res.node_Rparent_state .== "ready"
		TF = (TF1 + TF2 + TF3) .== 3
		res.node_state[TF] .= "ready_for_nodeOp"
	
		# Update nodes when the branches above finish
		indexes_ready = findall(res.node_state .== "ready_for_nodeOp")
		for current_nodeIndex in indexes_ready
			# Spawn a node operation
			#push!(tasks, @spawn nodeOp(current_nodeIndex, res))
			nodeOp(current_nodeIndex, res)
		end
	
		# Check if we are done?
		are_we_done = count_nodes_finished(res.node_state) >= res.numNodes
		
		# Error trap
		if (iteration_number >= max_iterations)
			txt = join(["Error in iterative_downpass(): iteration_number ", string(iteration_number), " exceeded max_iterations. Probably your loop is not concluding, or you have a massively huge tree or slow calculation, and need to set max_iterations=Inf."], "")
			error(txt)
		end
		
		# Test for concluding the while loop
		are_we_done && break
	end
	
	# This breaks it for some reason:
	# ERROR: setfield! immutable struct of type Res cannot be changed
	#global res.number_of_whileLoop_iterations = iteration_number

	print_num_iterations = false
	if print_num_iterations
		txt = join(["\nFinished at iteration_number ", string(iteration_number), "."], "")
		print(txt)
		print("\n")
	end
	
	# Final run diagnostics
	diagnostics[2] = Dates.now()
	diagnostics[3] = diagnostics[2]-diagnostics[1]
	total_calctime_in_sec = (diagnostics[2]-diagnostics[1]).value / 1000
	
	res.calctime_iterations[1] = total_calctime_in_sec
	res.calctime_iterations[2] = iteration_number / 1.0
	
	return(total_calctime_in_sec, iteration_number)
end # END iterative_downpass!




"""
Iterate through the "res" object many times to complete the downpass, spawning jobs along the way
Non-parallel version (no istaskdone, etc.)
"""
function iterative_downpass_nonparallel!(res; max_iterations=10^10, num_iterations=10000000)
	diagnostics = collect(repeat([Dates.now()], 3))
	diagnostics[1] = Dates.now()
	
	# Setup
	current_nodeIndex = res.root_nodeIndex

	# Check number of threads
	numthreads = Threads.nthreads()
	parallel_TF = numthreads > 1
	tasks = Any[]
	tasks_fetched_TF = Bool[]
	are_we_done = false

	iteration_number = 0
	while(are_we_done == false)
		iteration_number = iteration_number+1
		# As long as all the nodes are not done,
		# check for "ready" nodes
		# When they finish, change to "done"
		indexes_ready = findall(res.node_state .== "ready_for_branchOp")
		for current_nodeIndex in indexes_ready
			# Before spawning, do some checks
			res.node_state[current_nodeIndex] = "calculating_branchOp"
			# Check for root; no calculation on root branch for now
			if current_nodeIndex == res.root_nodeIndex
				res.node_state[current_nodeIndex] = "done"
				return()
			end

			# Spawn a branch operation, and a true-false of whether they are fetched
			res.calc_spawn_start[current_nodeIndex] = Dates.now()
			print(join(["\nbranchOp on current_nodeIndex=", string(current_nodeIndex)], ""))
# 			if (parallel_TF == true)
# 				push!(tasks, @spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
# 			else
				tmp_results = branchOp(current_nodeIndex, res, num_iterations)
				push!(tasks, tmp_results)
# 			end
			push!(tasks_fetched_TF, false)
		end
	
		# Check which jobs are done, fetch them, and update status of that node
		num_tasks = length(tasks)
		for i in 1:num_tasks
			if (tasks_fetched_TF[i] == false)
				#if (istaskdone(tasks[i]) == true)
					# Get the results
					calc_end_time = Dates.now()
# 					if (parallel_TF == true)
# 						(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time) = fetch(tasks[i])
# 					else
						(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time) = tasks[i]
# 					end
					# Store run information
					res.calc_start_time[spawned_nodeIndex] = calc_start_time
					res.calc_end_time[spawned_nodeIndex] = calc_end_time
					res.calc_duration[spawned_nodeIndex] = (calc_end_time - calc_start_time).value / 1000.0
					tasks_fetched_TF[i] = true
					
					# Record information
					res.thread_for_each_branchOp[spawned_nodeIndex] = tmp_threadID
					res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = nodeData_at_bottom
					# Get the ancestor nodeIndex
					uppass_edgematrix = res.uppass_edgematrix
					TF = uppass_edgematrix[:,2] .== spawned_nodeIndex
					parent_nodeIndex = uppass_edgematrix[TF,1][1]

					# Get the left daughter nodeIndex (1st in the uppass_edgematrix)
					edge_rows_TF = uppass_edgematrix[:,1] .== parent_nodeIndex
					left_nodeIndex = uppass_edgematrix[edge_rows_TF,2][1]
					right_nodeIndex = uppass_edgematrix[edge_rows_TF,2][2]

					# Update the state of the parent_node's daughters
					if (spawned_nodeIndex == left_nodeIndex)
						res.node_Lparent_state[parent_nodeIndex] = "ready"
					end
					if (spawned_nodeIndex == right_nodeIndex)
						res.node_Rparent_state[parent_nodeIndex] = "ready"
					end

					# Update the state of the current node
					res.node_state[spawned_nodeIndex] = "done"
				#end
			end
		end
	
		# Update which nodes have had both parents complete
		TF1 = res.node_state .== "not_ready"
		TF2 = res.node_Lparent_state .== "ready"
		TF3 = res.node_Rparent_state .== "ready"
		TF = (TF1 + TF2 + TF3) .== 3
		res.node_state[TF] .= "ready_for_nodeOp"
	
		# Update nodes when the branches above finish
		indexes_ready = findall(res.node_state .== "ready_for_nodeOp")
		for current_nodeIndex in indexes_ready
			# Spawn a node operation
			#push!(tasks, @spawn nodeOp(current_nodeIndex, res))
			nodeOp(current_nodeIndex, res)
		end
	
		# Check if we are done?
		are_we_done = count_nodes_finished(res.node_state) >= res.numNodes
		
		# Error trap
		if (iteration_number >= max_iterations)
			txt = join(["Error in iterative_downpass_nonparallel(): iteration_number ", string(iteration_number), " exceeded max_iterations. Probably your loop is not concluding, or you have a massively huge tree or slow calculation, and need to set max_iterations=Inf."], "")
			error(txt)
		end
		
		# Test for concluding the while loop
		are_we_done && break
	end
	
	# This breaks it for some reason:
	# ERROR: setfield! immutable struct of type Res cannot be changed
	#global res.number_of_whileLoop_iterations = iteration_number

	print_num_iterations = false
	if print_num_iterations
		txt = join(["\nFinished at iteration_number ", string(iteration_number), "."], "")
		print(txt)
		print("\n")
	end
	
	# Final run diagnostics
	diagnostics[2] = Dates.now()
	diagnostics[3] = diagnostics[2]-diagnostics[1]
	total_calctime_in_sec = (diagnostics[2]-diagnostics[1]).value / 1000
	
	res.calctime_iterations[1] = total_calctime_in_sec
	res.calctime_iterations[2] = iteration_number / 1.0
	
	return(total_calctime_in_sec, iteration_number)
end # END iterative_downpass_nonparallel!





















end # end of module