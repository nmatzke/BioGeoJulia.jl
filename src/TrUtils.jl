module TrUtils
using DataFrames
#using RCall

export getwd
export setwd
export recursive_find
export include_jls
export source
export seq
export Rchoose
export Rcbind
export Rrbind
export paste
export paste0
export type
export class
export slashslash
export addslash
export df_to_Rdata
export Reval
export Rdput
export Rnames
export saveopen
export Rnrow
export Rncol
export Rsize
export Rorder
export headLR
export flat2
export single_element_array_to_scalar
export headf
export moref
export scr2str

# R-like utilities, and other short functions

# Handy aliases

# getwd
function getwd()
	pwd()
end

# setwd
function setwd(path=expanduser("~"))
	cd(path)
end


# Find all code *.jl files in a package

"""
package_path = "/GitHub/BioGeoJulia.jl"
recursive_find(package_path)
include_jls(package_path)
"""

function recursive_find(package_path)
	# Look in the source directory, "src"
	srcpath = joinpath(package_path, "src")
	srcfiles = readdir(srcpath)
	jl_files = []  # empty array
	for fn in srcfiles
		TF = endswith(fn, ".jl")
		if (TF == true)
			tmpfn = joinpath(srcpath, fn)
			push!(jl_files, tmpfn)
		end
	end
	
	# Look in the tests directory, "test"
# 	srcpath = joinpath(package_path, "test")
# 	srcfiles = readdir(srcpath)
# 	for fn in srcfiles
# 		TF = endswith(fn, ".jl")
# 		if (TF == true)
# 			tmpfn = joinpath(srcpath, fn)
# 			push!(jl_files, tmpfn)
# 		end
# 	end
	
	return jl_files
end

"""
package_path = "/GitHub/BioGeoJulia.jl"
recursive_find(package_path)
include_jls(package_path)
"""
function include_jls(package_path)
	srcfiles = recursive_find(package_path)
	for fn in srcfiles
		include(fn)
	end
end

# source
function source(str)
	include(str)
end

# seq
function seq(from, to, by=1)
	return(collect(from:by:to))
end


# choose
# n choose k
function Rchoose(n,k)
	return(binomial(n,k))
end


# cbind()
function Rcbind(A...)
	hcat(A...)
end

# rbind
# c()
# concatenate
function Rrbind(A...)
	vcat(A...)
end

# paste
function paste(array_of_strings, delim)
	newtxt = join(array_of_strings, delim)
	return(newtxt)
end

# paste0
function paste0(array_of_strings, delim=0)
	newtxt = join(array_of_strings, delim)
	return(newtxt)
end

# type
function type(obj)
	typeof(obj)
end

# class
function class(obj)
	typeof(obj)
end

# Convert any multiple slashes to single slashes
function slashslash(txt)
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	return(txt)
end

# Add a slash to the end of a string, if it is not there
# (Unless the txt string is "", then return "")
function addslash(txt)
	if (txt == "")
		return("")
	end
	if (endswith(txt, "/") == false)
		txt = join([txt, "/"])
	end
	return(txt)
end


# Save a julia DataFrame to an R data.frame
# as an Rdata file that can be easily loaded.
# Source: https://stackoverflow.com/questions/28084403/saving-julia-dataframe-to-read-in-r-using-hdf5/57903489#57903489
using DataFrames
#using RCall
function df_to_Rdata(df; fn="dfjulia.RData", path=expanduser("~"))
	# Create the output fn
	pathfn = slashslash(join([addslash(path), fn], ""))
		
	# R environment in a session started from Julia
	g = globalEnv  # This requires "using RCall" previously, or you get:
								 # ERROR: UndefVarError: globalEnv not defined
	reval(rparse("dfls <- NULL"))

	# add columns one at a time converting Julia vectors to R-types via   RCall.sexp
	#  https://github.com/JuliaStats/RCall.jl/blob/master/src/sexp.jl
	for cnm in DataFrames._names(df)
		g[:colcnm] = sexp(convert(Array, df[!,cnm]))
		reval(rparse("dfls\$$cnm <- colcnm"))
	end
	reval(rparse("df <- data.frame(dfls)"))
	
	# Make and run the command to save the .Rdata file
	# (the .Rdata file will load to object df in R)
	txt = join(["save(file='", pathfn, "', df)"], "")
	reval(rparse(txt))
	return(pathfn)
end


# Input an object from a string representation from repr()
# dget, load, eval()
# eval
"""
tmpmatrix = [3 1; 3 2; 5 3; 5 4; 7 5; 7 6]
tmpstr = repr(tmpmatrix)
tmpstr2 = eval(Meta.parse(tmpstr))
tmpstr2
"""
function Reval(tmpstr)
	eval(Meta.parse(tmpstr))
end

# Output an object to a string representation with repr()
# dput, dump, str()
"""
tmpstr = "Array{Any,1}[[1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]"
states_list = Reval(tmpstr)
tmpstr2 = Rdput(states_list)
"""
function Rdput(item)
	tmpstr = repr(item)
end


# fields / "names" of an object
# https://stackoverflow.com/questions/41687418/how-to-get-fields-of-a-julia-object
function Rnames(obj)
	fieldnames(typeof(obj))
end

# Send the just-done plot to PDF, and open
function saveopen(fn)
	print("Saving ", """'""", fn, """'""")
	savefig(fn)
	cmdtxt = join(["`", "open ", fn, "`"], "")
	print("""Running """, cmdtxt, "\n")
	run(`open $fn`)
end


function Rnrow(obj)
	return size(obj)[1]
end

function Rncol(obj)
	return size(obj)[2]
end

function Rsize(obj)
	return size(obj)
end

function Rorder(obj)
	return sortperm(obj)
end

# Print the left and rightmost columns of a table
function headLR(df, num_startcols=4, num_endcols=4)
	ncols = Rncol(df)
	startcols = collect(1:num_startcols)
	endcols = collect((ncols-(num_endcols-1)):ncols)
	colnums = flat2(collect([startcols, endcols]))
	colnums
	print(df[:,colnums])
end


# Flattens an array of arrays into a vector
# Similar to R's unlist()
function flat2(arr)
    rst = Any[]
    grep(v) = for x in v
        if isa(x, Array) grep(x) else push!(rst, x) end
    end
    grep(arr)
    rst
end


# Convert a single-element array to scalar
# Julia often produces single-element arrays. 
# To convert to scalar, just take item [1]
# https://stackoverflow.com/questions/39079428/1-element-array-to-scalar-in-julia
function single_element_array_to_scalar(tmparray)
	if length(tmparray) != 1
		txt = ["STOP ERROR in single_element_array_to_scalar().\nThe input 'tmparray' has to be of length 1, but length(tmparray)=", string(length(tmparray)), ".\nPrinting input tmparray...\n"]
		errortxt = join(txt, "")
		println(errortxt)
		print(tmparray)
		error(errortxt)
	end
	
	# If check passed, go ahead.
	tmpscalar = tmparray[1]
	return tmpscalar
end



# Print the file to screen, with line numbers
function headf(fn; numlines=5)
	open(fn, "r") do f
		for (i,ln) in enumerate(eachline(f))
			if i > numlines
				break
			end
			println("$i $ln")
		end
	end
end

# Print the file to screen, with line numbers
function moref(fn)
	open(fn, "r") do f
	 for (i,ln) in enumerate(eachline(f))
		 println("$i $ln")
	 end
	end
end


function scr2str(obj)
	io = IOBuffer()
	show(io, "text/plain", obj)
	str = String(take!(io))
	return str
end




end