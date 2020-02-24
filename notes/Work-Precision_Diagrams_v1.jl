
#######################################################
# Work-Precision Diagrams on ClaSSE Es & Ds calculations
#######################################################

# Modified from:
# 100 Independent Linear Work-Precision Diagrams
# Chris Rackauckas
# https://benchmarks.juliadiffeq.org/html/NonStiffODE/linear_wpd.html

using Pkg
using OrdinaryDiffEq, Sundials, DiffEqDevTools, Plots, ODEInterfaceDiffEq, ODE, LSODA
Pkg.add(PackageSpec(url="https://github.com/JuliaDiffEq/deSolveDiffEq.jl"))
using deSolveDiffEq


using Random
Random.seed!(123)
gr()

abstols = 1.0 ./ 10.0 .^ (3:13)
reltols = 1.0 ./ 10.0 .^ (0:10);