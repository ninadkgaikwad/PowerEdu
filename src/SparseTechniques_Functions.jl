# SparseTechniquesInPowerSystems.jl

using SparseArrays
using DataFrames

function nnzRowConstructor(compElem::DataFrameRow, elemNum::Int64)
	nnzElem = DataFrame(ID = elemNum, Val = compElem.Val, NROW = compElem.i, NCOL = compElem.j, NIR = -1, NIC = -1)
	@show typeof(nnzElem)
	return first(nnzElem)
end

function sparmat(compMatrix::DataFrame)
	N = maximum([compMatrix.i compMatrix.j])
	(firs, fics) = (repeat([-1], N), repeat([-1], N))

	NVec = DataFrame(FIR = firs, FIC = fics)
	nnzVec = DataFrame(ID = Int64[], Val = ComplexF64[], NROW = Int64[], NCOL = Int64[], NIR = Int64[], NIC = Int64[])

	numElems = size(compMatrix, 1)
	
	for elemNum in 1:numElems
		compElem = compMatrix[elemNum, :]

		nnzElem = nnzRowConstructor(compElem, elemNum)
		NVec, nnzVec = updateSparse(NVec, nnzVec, nnzElem, type="replace")
	end

	return NVec, nnzVec
end



function updateSparse(NVec::DataFrame, nnzVec::DataFrame, nnzElem::DataFrameRow;
	type::String = "replace")

	updateFlag = false #Check if the value is actually updated instead of added the sparse matrix
	numExistingElems = size(nnzElem, 1)
	FIR = NVec.FIR
	FIC = NVec.FIC
	row = nnzElem.NROW
	col = nnzElem.NCOL
	elID = nnzElem.ID
	elVal = nnzElem.Val

	# Check for row placement
	if FIR[row] == -1 # First element to be inserted into that row
		println("First Element!")
		FIR[row] = elID
		nnzElem.NIR = -1 #auto, unneeded
	else # There exists at least one element already inserted into that row
		x = FIR[row]
		println("NOOO! $x is NOT equal to $(-1)")
		@show incumbentID = FIR[row]
		incumbent = nnzVec[incumbentID, :]
		@show typeof(incumbent)

		if col < incumbent.NCOL # New elem comes before incumbent in that row => replace it in FIR and shift incumbent in the nnzVec
			FIR[row] = elID
			el.NIR = incumbentID
		elseif col == incumbent.NCOL # Same-same? Replace or add?
			if type == "replace"
				nnzVec.Val[incumbentID] = elVal
				updateFlag = true
			elseif type == "add"
				nnzVec.Val[incumbentID] += elVal
			else
				error("Not prepared for this scenario.")
			end
		elseif col > incumbent.NCOL # elem comes only after the incumbent, check for next element in line
			prevIncumbentID = incumbentID
			incumbentID = incumbent.NIR # going to the next element in the row
			incumbent = nnzVec[incumbentID, :]

			StillFindingAPlace = true # Initalizer for a loop, which will keep checking the next element in the row until elem gets its own place
			while StillFindingAPlace
				if incumbentID == -1 # elem is the last in the row
					nnzVec.NIR[prevIncumbentID] = elID
					nnzElem.NIR[prevIncumbentID] = -1 # auto
					StillFindingAPlace = false
				else # more elements to check in the row
					incumbent = nnzVec[incumbentID, :]

					if col < incumbent.NCOL
						nnzElem.NIR = incumbentID
						nnzVec[prevIncumbentID] = elID
					elseif col == incumbent.NCOL # Same-same? Either replace or add
						if type == "replace"
							nnzVec.Val[incumbentID] = elVal
							updateFlag = true
							StillFindingAPlace = false
						elseif type == "add"
							nnzVec.Val[incumbentID] += elVal
						else
							error("Not prepared for this scenario")
						end
					elseif col < incumbent.NCOL
						prevIncumbentID = incumbentID
						incumbentID = incumbent.NIR
						incumbent = nnzVec[incumbentID, :]
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
		push!(nnzVec, nnzElem)
		println("Another element added to the sparse matrix.")
	end

	return NVec, nnzVec
end


values = vec([-1, -2, 2, 8, 1, 3, -2, -3, 2, 1, 2, -4]);
rows = vec([1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 5]);
cols = vec([1, 3, 1, 2, 4, 3, 5, 2, 3, 1, 2, 5]);
compMatrix = DataFrame(Val = values, i = rows, j = cols);

NVec, nnzVec = sparmat(compMatrix)

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