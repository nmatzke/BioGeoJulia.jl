#######################################################
# Parsers
#######################################################

module Parsers

print("\n\nStarting module 'Parsers'...loading dependencies...\n")
using BenchmarkTools # for @time
using InvertedIndices # for Not
using LSODA           # for lsoda()
using Sundials        # for CVODE_BDF(linear_solver=:GMRES)
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

# (1) List all function names here:
export say_hello3, getranges_from_LagrangePHYLIP

#######################################################
# Temporary file to store functions under development
#
# Start with:
# 
# Setup:

"""
cd("/GitHub/BioGeoJulia.jl/notes/")
include("Parsers.jl")
"""
#######################################################


#######################################################
# (2) write the functions here
#######################################################

say_hello3() = println("Hello dude3!")
say_hello3()

"""
# Goal: Read in various text files,
# starting with a C++Lagrange / BioGeoBEARS-format
# geography file, e.g.:
#
# =================
# 19	4 (K O M H)
# P_mariniana_Kokee2	1000
# P_mariniana_Oahu	0100
# P_mariniana_MauiNui	0010
# ...
# =================
#
# This is basically a "PHYLIP"-formatted file.


include("/GitHub/BioGeoJulia.jl/notes/Parsers.jl")
import .Parsers
lgdata_fn = "/GitHub/BioGeoJulia.jl/Rsrc/Psychotria_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)
"""


"""
Function to read in Lagrange/BioGeoBEARS-type PHYLIP-formatted
 geography input files.

lgdata_fn = lagrange-style geography input file
block_allQs = give error if a row has all-question-marks
"""

function getranges_from_LagrangePHYLIP(lgdata_fn; block_allQs=true)
	#lgdata_fn = "/GitHub/BioGeoJulia.jl/Rsrc/Psychotria_geog.data"

	# Obtain a file handle
	fhandle = open(lgdata_fn)

	# Read the lines to an array
	lines = readlines(fhandle)
	close(fhandle)
	
	# Count the species read in
	spnum = 0
	
	# We have to assume that the first line is a header
	for i in 1:length(lines)
		# If the line is blank, skip it
		line = lines[i]
		if (length(line) == 0)
			continue
		end
		#println(i, ": ", line)

		# Parse the header line
		if (i == 1)
			# Split on "("
			parts = split(line, "(")
			# remove whitespace, convert to Integer using parse.()
			parts1 = strip(parts[1])
			global (numtaxa, numareas) = parse.(Int, split(strip(parts1)))

			# Part 2: remove trailing ")", parse the rest to areas list
			parts2 = strip(strip(strip(parts[2]), [')']))
			# Split on whitespace to get the list of areas
			global areas_list = split(parts2)

			# Check that the given number matches the actual number of areas
			if (length(areas_list) != numareas)
				txt = paste0(["STOP ERROR in getranges_from_LagrangePHYLIP(). In your input file,\n'", lgdata_fn, "',\nthe given number of areas in line 1 (numareas=", numareas, ") does not equal the length of the areas_list (", parts2, ").\nPlease correct your input geography file and retry. Have a nice day."])
				error(txt)
			end # end error check
	 
			# Set up the output matrix
			global sp_names = collect(repeat([""], numtaxa))
			tmp_areas = collect(repeat(["0"], numareas))
			global sp_areas = collect(repeat([tmp_areas], numtaxa))

			continue # Go to next round of for-loop without finishing this one
		end # END if (i == 1)

		# All of the other lines, fill in sp_names and sp_areas
		words = strip.(split(line))
		spnum = spnum + 1
	
		# Error check
		if (spnum > numtaxa)
			txt = paste0(["STOP ERROR in getranges_from_LagrangePHYLIP(). While reading in species #", spnum, ", you exceeded the limit declared in the first line of your input geography file, where numtaxa=", numtaxa, " species. Please correct your input geography file and retry. Have a nice day."])
			error(txt)
		end

		
		sp_names[spnum] = words[1]
		# Split areanums to "1", "0", etc., parse to Integers
		#areanums_for_this_species = parse.(Int, split(words[2], ""))
		areanums_for_this_species = split(words[2], "")
		#print(areanums_for_this_species)
		# [:] avoids creating a linked reference
		sp_areas[spnum] = areanums_for_this_species[:]

	end

	# Error check
	if (spnum != numtaxa)
		txt = paste0(["STOP ERROR in getranges_from_LagrangePHYLIP(). While reading in species, only ", spnum, " taxa were read in, however the first line of your input geography file declared there would be numtaxa=", numtaxa, " species. Please correct your input geography file and retry. Have a nice day."])
		error(txt)
	end

	# DataFrame of the area presence/absences
	# Flatten to vector, reshape to array, transpose to get
	# back original inputs in matrix form
	sp_areas2 = permutedims(reshape(flat2(sp_areas), numareas, numtaxa))
	sp_areas2_Matrix = convert(Matrix, sp_areas2)
	sp_areas3_Matrix = Rcbind(sp_names, sp_areas2_Matrix)
	geog_df = convert(DataFrame, sp_areas3_Matrix)
	new_names = Rrbind(["tipnames"], areas_list)
	geog_df = rename(geog_df, new_names)

	geog_df[!,:K]
	geog_df[!,:tipnames]

	return geog_df
end # END getranges_from_LagrangePHYLIP




end # ENDING Parsers