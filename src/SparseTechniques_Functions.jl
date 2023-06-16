# SparseTechniquesInPowerSystems.jl

using SparseArrays
using DataFrames

include("Helper_Functions.jl")

function nnzRowConstructor(compElem::DataFrameRow)
	nnzElem = DataFrame(ID = -1, Val = compElem.Val, NROW = compElem.i, NCOL = compElem.j, NIR = -1, NIC = -1)
	return first(nnzElem) # Otherwise, Julia will interpret it as a DataFrame and NOT a DataFrameRow
	# The reason we want a DataFrameRow is that when we invoke nnzElem.NROW, it extracts the inner element, instead of a 1-vector
end

function sparmat(compMatrix::DataFrame;
	verbose::Bool = false)
	N = maximum([compMatrix.i compMatrix.j])
	(firs, fics) = (repeat([-1], N), repeat([-1], N))

	NVec = DataFrame(FIR = firs, FIC = fics)
	nnzVec = DataFrame(ID = Int64[], Val = ComplexF64[], NROW = Int64[], NCOL = Int64[], NIR = Int64[], NIC = Int64[])

	numElems = size(compMatrix, 1)
	
	for elemNum in 1:numElems
		compElem = compMatrix[elemNum, :]

		nnzElem = nnzRowConstructor(compElem)

		NVec, nnzVec = updateSparse(NVec, nnzVec, nnzElem, type="replace", verbose=verbose)
	end

	return NVec, nnzVec
end



function updateSparse(NVec::DataFrame, nnzVec::DataFrame, nnzElem::DataFrameRow;
	type::String = "replace",
	verbose::Bool = false)

	updateFlag = false #Check if the value is actually updated instead of added the sparse matrix
	numExistingElems = size(nnzVec, 1)
	myprintln(verbose, "Currently the Sparse Matrix has $numExistingElems elements.")

	elID = numExistingElems + 1
	FIR = NVec.FIR
	FIC = NVec.FIC
	row = nnzElem.NROW
	col = nnzElem.NCOL
	# elID = nnzElem.ID
	elVal = nnzElem.Val

	# Check for row placement
	if FIR[row] == -1 # First element to be inserted into that row
		myprintln(verbose, "This elem is the first in row $row !")
		FIR[row] = elID
		nnzElem.NIR = -1 #auto, unneeded
	else # There exists at least one element already inserted into that row
		incumbentID = FIR[row]
		myprintln(verbose, "Row $row already has an element, with ID: $incumbentID")
		
		incumbent = nnzVec[incumbentID, :]

		if col < incumbent.NCOL # New elem comes before incumbent in that row => replace it in FIR and shift incumbent in the nnzVec
			myprintln(verbose, "Our elem comes before element $incumbentID, can topple it!")
			FIR[row] = elID
			el.NIR = incumbentID
		elseif col == incumbent.NCOL # Same-same? Replace or add?
			myprintln(verbose, "Stand down. It's a draw between our elem and element $incumbentID.")
			if type == "replace"
				nnzVec.Val[incumbentID] = elVal
				myprintln(verbose, "Replacing/Updating the element's previous value with the new elem's value.")
				updateFlag = true
			elseif type == "add"
				nnzVec.Val[incumbentID] += elVal
				myprintln(verbose, "Added the the new elem's value to the incumbent's value.")
			else
				error("Not prepared for this scenario.")
			end
		elseif col > incumbent.NCOL # elem comes only after the incumbent, check for next element in line
			myprintln(verbose, "Our elem will come after the incumbent element $incumbentID. Continue Searching.")
			prevIncumbentID = incumbentID
			incumbentID = incumbent.NIR # going to the next element in the row

			StillFindingAPlace = true # Initalizer for a loop, which will keep checking the next element in the row until elem gets its own place
			while StillFindingAPlace
				if incumbentID == -1 # elem is the last in the row
					myprintln(verbose, "Our elem will sit right after the FIR element $prevIncumbentID !")
					nnzVec.NIR[prevIncumbentID] = elID
					nnzElem.NIR = -1 # auto
					StillFindingAPlace = false
				else # more elements to check in the row
					incumbent = nnzVec[incumbentID, :]
					myprintln(verbose, "Not the second element to be added to this row either. Keep searching.")
					if col < incumbent.NCOL
						myprintln(verbose, "Our element can topple the incumbent element $incumbentID.")
						nnzElem.NIR = incumbentID
						nnzVec[prevIncumbentID] = elID
						StillFindingAPlace = false

					elseif col == incumbent.NCOL # Same-same? Either replace or add
						myprintln(verbose, "Stand down. It's a draw between our elem and element $incumbentID.")
						if type == "replace"
							nnzVec.Val[incumbentID] = elVal
							myprintln(verbose, "Replacing/Updating the element's previous value with the new elem's value.")
							updateFlag = true
						elseif type == "add"
							nnzVec.Val[incumbentID] += elVal
							myprintln(verbose, "Added the the new elem's value to the incumbent's value.")
						else
							error("Not prepared for this scenario")
						end
						StillFindingAPlace = false

					elseif col > incumbent.NCOL
						myprintln(verbose, "Not coming before incumbent element $incumbentID. Keep searching.")
						prevIncumbentID = incumbentID
						incumbentID = incumbent.NIR

					else
						error("Shouldn't be possible.")
					end
				end
			end
		
		else
			error("Shouldn't be possible!")
		end
	end

	if updateFlag == false
		myprintln(verbose, "Another element to be added to the sparse matrix.")
		myprintln(verbose, "The element:")
		myprintln(verbose, nnzElem)
		myprintln(verbose, "nnzVec before:")
		myprintln(verbose, nnzVec)
		push!(nnzVec, nnzElem)
		myprintln(verbose, "nnzVec after:")
		myprintln(verbose, nnzVec)
	end

	return NVec, nnzVec
end


values = vec([-1, -2, 2, 8, 1, 3, -2, -3, 2, 1, 2, -4]);
rows = vec([1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 5]);
cols = vec([1, 3, 1, 2, 4, 3, 5, 2, 3, 1, 2, 5]);
compMatrix = DataFrame(Val = values, i = rows, j = cols);

NVec, nnzVec = sparmat(compMatrix, verbose = false)

"""
    compressed2Full(compMatrix::DataFrame)

Converts a compressed matrix representation stored in a DataFrame into a full 
matrix.

## Arguments
- `compMatrix::DataFrame`: A DataFrame representing the compressed matrix. 
It should have three columns: `i`, `j`, and `val`. Column `i` contains the 
row indices, column `j` contains the column indices, and column `val` contains 
the corresponding values.

## Returns
- `fullMatrix::Matrix`: A full matrix representation of the input compressed 
matrix.

## Dependencies
This function requires the following packages to be imported:
- `SparseArrays`: Provides support for sparse matrix operations.
- `DataFrames`: Provides support for working with tabular data in a 
DataFrame format.

## Example
```julia
using SparseArrays
using DataFrames

# Define the compressed matrix
values = vec([-1, -2, 2, 8, 1, 3, -2, -3, 2, 1, 2, -4])
rows = vec([1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 5])
cols = vec([1, 3, 1, 2, 4, 3, 5, 2, 3, 1, 2, 5])
compMatrix = DataFrame(Val = values, i = rows, j = cols)

# Convert the compressed matrix to a full matrix
matFull = compressed2Full(compMatrix)
"""

function compressed2Full(compMatrix::DataFrame)
	# Use SparseArrays's sparse function to conveniently convert the compressed 
	# matrix (i, j, Val) into Compressed Storage Column CSC format
	# I don't care how it does it. 
	# This function is only to be called for testing purposes anyway.
    sparseMatrix = sparse(compMatrix.i, compMatrix.j, compMatrix.Val)
	# Convert the sparse matrix into the full matrix.
    fullMatrix = Matrix(sparseMatrix)
	return fullMatrix
end