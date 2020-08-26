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

include("/GitHub/BioGeoJulia.jl/src/TreePass.jl")
import .TreePass

# Repeat calculation in Julia
include("/GitHub/BioGeoJulia.jl/notes/ModelLikes.jl")
import .ModelLikes

########################



readTopology("(((chimp:1,human:1):1,gorilla:2):1,orang:3);")
length(tr.node)
tr.edge[2]

newnode, newedge = PhyloNetworks.breakedge!(tr.edge[2], tr)

# getting and issue here. 
# Error reads: ERROR: BoundsError: attempt to access Int64
#  at index [2]

# Also cannot convert tr.edge[2] to an float. is it trying to look for an Int but only finding a floater? hmm

"""
PhyloNetwork's practice code reads as: 

But even in their practice, it breaks at newnode, newedge

net = readTopology("(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));")
length(net.node)
net.edge[4]
newnode, newedge = PhyloNetworks.breakedge!(net.edge[4], net)
length(net.node)
newedge # new edge 21 goes from node -8 and 11 (new)
net.edge[4] # original edge 4 now goes from node 11 (new) to 3
writeTopology(net)

"""