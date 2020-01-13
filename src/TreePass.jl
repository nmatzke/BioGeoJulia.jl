module TreePass
using DataFrames
using PhyloNetworks  # for e.g. readTopology()
export get_nodenumbers_above_node, get_postorder_nodenumbers_above_node, initialize_edgematrix, get_pruningwise_postorder_edgematrix, get_LR_uppass_edgematrix, get_LR_downpass_edgematrix, get_LR_uppass_nodeIndexes, get_LR_downpass_nodeIndexes, get_Rnodenums, get_nodeIndex_PNnumber, get_nodeIndex_from_PNnumber, prt, get_taxa_descending_from_each_node, isTip_TF, get_NodeIndexes_from_edge, get_NodeIndex_df_by_tree_edges, get_node_heights, get_node_ages






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








end # end of module