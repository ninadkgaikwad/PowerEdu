# SparseTechniquesInPowerSystems.jl
using Revise
using SparseArrays
using DataFrames
using Test
using CSV

include("Helper_Functions.jl")

const SparseMatrix = NamedTuple{(:NVec, :MVec, :nnzVec), Tuple{DataFrame, DataFrame, DataFrame}}

""""
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


"""
    spar2Full(sparMat::SparseMatrix;
    readMethod::String = "row-wise",
    verbose::Bool = false)

Convert a sparse matrix represented by sparMat to a regular (full) matrix. 
It is literally the inverse function of `sparmat`.

## Arguments
- `sparMat::SparseMatrix`: Named tuple representing the sparse matrix, with fields `NVec`, `MVec`, and `nnzVec`.
- `readMethod::String`: Optional. The method to read the elements of the sparse matrix. Default is "row-wise". Possible values are "row-wise", "col-wise", and "nnzVecElems-wise".
- `verbose::Bool`: Optional. If set to `true`, print verbose output. Default is `false`.

## Returns
- `mat::Matrix{ComplexF64}`: The converted regular matrix.

## Details
The `spar2Full` function converts the sparse matrix represented by `sparMat` to a regular matrix. The sparse matrix is defined by the non-zero elements stored in `nnzVec` DataFrame, with the row and column information stored in `NVec` and `MVec` DataFrames, respectively.

The conversion can be done in different ways depending on the `readMethod` argument:
- `"row-wise"`: The elements are read row by row, and the corresponding positions in the regular matrix are filled.
- `"col-wise"`: The elements are read column by column, and the corresponding positions in the regular matrix are filled.
- `"nnzVecElems-wise"`: The elements are read from `nnzVec` DataFrame in the order they appear, and the corresponding positions in the regular matrix are filled.

## Example
```julia
sparMat = (NVec=DataFrame(FIR=[1, 2, 4], NIR=[2, -1, -1]),
MVec=DataFrame(FIC=[1, 3, 2, -1]), 
nnzVec=DataFrame(ID=[1, 2, 3],
NROW=[1, 2, 3],
NCOL=[1, 2, 3],
Val=[1.0, 2.0, 3.0],
NIR=[2, -1, -1],
NIC=[-1, -1, -1]))
mat = spar2Full(sparMat, readMethod="row-wise", verbose=true)
```
```julia
3×3 Matrix{ComplexF64}:
1.0+0.0im  0.0+2.0im  0.0+0.0im
0.0+0.0im  0.0+0.0im  2.0+0.0im
0.0+0.0im  0.0+0.0im  0.0+3.0im
```
```julia
mat = spar2Full(sparMat, readMethod="col-wise", verbose=true)
```
```julia
3×3 Matrix{ComplexF64}:
1.0+0.0im  0.0+2.0im  0.0+0.0im
0.0+0.0im  0.0+0.0im  2.0+0.0im
0.0+0.0im  0.0+0.0im  0.0+3.0im
```
```julia
mat = spar2Full(sparMat, readMethod="nnzVecElems-wise", verbose=true)
```
```julia
3×3 Matrix{ComplexF64}:
1.0+0.0im  0.0+2.0im  0.0+0.0im
0.0+0.0im  0.0+0.0im  2.0+0.0im
0.0+0.0im  0.0+0.0im  0.0+3.0im
```
"""
function spar2Full(sparMat::SparseMatrix;
    readMethod::String = "row-wise",
    verbose::Bool = false)

    NVec = sparMat.NVec
    MVec = sparMat.MVec
    nnzVec = sparMat.nnzVec
    T = eltype(nnzVec.Val)

    FIR = NVec.FIR
    FIC = MVec.FIC
    N = length(FIR)
    M = length(FIC)
    nnz = length(nnzVec.ID)
    mat = zeros(T, N, M)

    if readMethod == "row-wise"
        for row = 1:N
            elemIdx = FIR[row]
            while elemIdx != -1
                nnzElem = nnzVec[elemIdx, :]
                i = nnzElem.NROW
                j = nnzElem.NCOL
                mat_ij = nnzElem.Val
                mat[i, j] = mat_ij
                elemIdx = nnzElem.NIR
            end
        end
    elseif readMethod == "col-wise"
        for col = 1:M
            elemIdx = FIC[col]
            while elemIdx != -1
                nnzElem = nnzVec[elemIdx, :]
                i = nnzElem.NROW
                j = nnzElem.NCOL
                mat_ij = nnzElem.Val
                mat[i, j] = mat_ij
                elemIdx = nnzElem.NIC
            end
        end
    elseif readMethod == "nnzVecElems-wise"
        for elemNum = 1:nnz
            nnzElem = nnzVec[elemNum, :]
            i = nnzElem.NROW
            j = nnzElem.NCOL
            mat_ij = nnzElem.Val
            mat[i, j] = mat_ij
        end
    else
        error("Unknown read method for conversion of sparse matrix
        to full")
    end
        
    return mat
end

"""
nnzRowConstructor(compElem::DataFrameRow)

The `nnzRowConstructor` function constructs a single row of the `nnz` DataFrame. It initializes the row with information extracted from the `compElem` DataFrameRow.

Arguments:
- `compElem::DataFrameRow`: A DataFrameRow containing the information to initialize the `nnz` DataFrame row.

Returns:
- A DataFrameRow representing a single row of the `nnz` DataFrame, initialized with the appropriate values extracted from `compElem`.

Example:
```julia
compElem = DataFrame(ID = 1, Val = 0.5, i = 2, j = 3)
nnzRow = nnzRowConstructor(compElem)
println(nnzRow.ID)    # -1
println(nnzRow.Val)   # 0.5
println(nnzRow.NROW)  # 2
println(nnzRow.NCOL)  # 3
println(nnzRow.NIR)   # -1
println(nnzRow.NIC)   # -1
"""
function nnzRowConstructor(compElem::DataFrameRow)
	nnzElem = DataFrame(ID = -1, Val = compElem.Val, NROW = compElem.i, NCOL = compElem.j, NIR = -1, NIC = -1)
	return first(nnzElem) # Otherwise, Julia will interpret it as a DataFrame and NOT a DataFrameRow
	# The reason we want a DataFrameRow is that when we invoke nnzElem.NROW, it extracts the inner element, instead of a 1-vector
end

"""
    sparmat(input::T where T<:Union{DataFrame, Matrix}; verbose::Bool = false)

The `sparmat` function creates a named tuple `sparMat` that contains three DataFrames, `NVec`, `MVec`, and `nnzVec`, which together represent a sparse matrix. It takes either a compressed matrix `compMatrix` (a DataFrame) or a regular (full) matrix as input. If the input is a regular matrix, it first converts it into `compMatrix` using the `full2comp` function.

Arguments:
- `input::T`: The input matrix, which can be either a compressed matrix represented as a DataFrame (`compMatrix`) or a regular (full) matrix.
- `verbose::Bool`: (optional) A Boolean value indicating whether to display verbose output during the construction process. Default is `false`.

Returns:
- A named tuple `sparMat = (NVec = NVec, MVec = MVec, nnzVec = nnzVec)` representing the three DataFrames `NVec`, `MVec`, and `nnzVec` that together represent the sparse matrix.

Example:
```julia
compMatrix = DataFrame(ID = [1, 2, 3], Val = [0.5, 0.3, 0.8], i = [2, 1, 3], j = [3, 2, 1])
sparMat = sparmat(compMatrix, verbose = true)
"""
function sparmat(input::T where T<:Union{DataFrame, Matrix};
    verbose::Bool=false)
    if isa(input, DataFrame)
        # Input is a compMat, so use it directly
        compMatrix = input
    elseif isa(input, Matrix)
        # Input is a regular (full) matrix, convert it to compMat using full2comp
        compMatrix = full2comp(input)
    else
        throw(ArgumentError("Unsupported input type. Must be either DataFrame (compMat) or Matrix"))
    end

    N = maximum(compMatrix.i)
    M = maximum(compMatrix.j)
    (firs, fics) = (repeat([-1], N), repeat([-1], M))

    dataType = eltype(compMatrix.Val)

    NVec = DataFrame(FIR = firs)
    MVec = DataFrame(FIC = fics)
	nnzVec = DataFrame(ID = Int64[], Val = dataType[], NROW = Int64[], NCOL = Int64[], NIR = Int64[], NIC = Int64[])

	numElems = size(compMatrix, 1)
	
	for elemNum in 1:numElems
		compElem = compMatrix[elemNum, :]

		nnzElem = nnzRowConstructor(compElem)

        sparMat = (NVec=NVec, MVec=MVec, nnzVec=nnzVec)
		NVec, MVec, nnzVec = updateSparse(sparMat, nnzElem, type="replace", verbose=verbose)
	end

    sparMat = (NVec=NVec, MVec=MVec, nnzVec=nnzVec) 
	return sparMat
end



"""
    resolveTie(nnzVec::DataFrame, incumbentID::Int64, elVal::ComplexF64; type="replace", verbose=false)

Resolves a tie between an incumbent element and a new element by either replacing the incumbent's value or adding the new element's value to it.

## Note
If `type` is set to `"maintain"`, the `nnzElem` value will not be updated. 
Only the other parameter values of `NVec` and `nnzVec` vectors will be changed.

## Arguments
- `nnzVec::DataFrame`: The DataFrame representing the sparse matrix.
- `incumbentID::Int64`: The ID of the incumbent element.
- `elVal::ComplexF64`: The value of the new element.
- `type::String`: (optional) The tie resolution strategy. Valid options are "replace" (default), "add", and "maintain".
- `verbose::Bool`: (optional) Whether to print verbose output. Default is `false`.

## Returns
The updated `nnzVec` DataFrame after resolving the tie.

# Examples
```julia
nnzVec = DataFrame(ID = 1:5, NROW = [1, 1, 2, 2, 3], NCOL = [1, 2, 1, 2, 1], Val = [1.0, 2.0, 3.0, 4.0, 5.0], NIR = [2, -1, 4, -1, -1], NIC = [-1, -1, -1, -1, -1])
nnzVec = resolveTie(nnzVec, 3, 2.0 + 1.5im, type="add", verbose=true)


```julia 
nnzVec = resolveTie(nnzVec, 3, 2.0 + 1.5im, type="maintain", verbose=true)
"""
function resolveTie(nnzVec::DataFrame, incumbentID::Int64, elVal::ComplexF64;
	type::String = "replace",
	verbose::Bool = false)
	
    if type == "replace"
        nnzVec.Val[incumbentID] = elVal
        myprintln(verbose, "Replacing/Updating the element's previous value with the new elem's value.")
    elseif type == "add"
        nnzVec.Val[incumbentID] += elVal
        myprintln(verbose, "Added the new elem's value to the incumbent's value.")
    else
        if type != "maintain"
            error("Not prepared for this scenario.")
        end
    end
    
    return nnzVec
end

"""
    checkElementIntoRow(sparMat::SparseMatrix,
    nnzElem::DataFrameRow;
	type::String="replace", verbose::Bool=false)

Check and insert an element into a row of a sparse matrix.

This function checks if the given `nnzElem` can be inserted into the row specified by `nnzElem.NROW`
in the `nnzVec` DataFrame. If the row is empty, the element is inserted as the first element. If there
is already an element present, the function determines the correct position for the element based on
the `nnzElem.NCOL` value. The position can be before the incumbent element, after it, or it may result
in a tie, which can be resolved using the `resolveTie` function.

## Arguments
- `sparMat::SparseMatrix`: A named tuple representing the sparse matrix with elements `NVec`, `MVec`, and `nnzVec`. `NVec` is a DataFrame containing additional information about the sparse matrix, `MVec` is a DataFrame containing additional information about the columns, and `nnzVec` is a DataFrame representing the non-zero elements of the sparse matrix.
- `nnzElem::DataFrameRow`: DataFrameRow representing the new element to be inserted.
- `type::String`: Optional. The type of tie resolution. Default is "replace". Possible values are "replace" and "add".
- `verbose::Bool`: Optional. If set to `true`, print verbose output. Default is `false`.

## Returns
- `sparMat::SparseMatrix`: Updated named tuple `sparMat` with the modified `NVec`, `MVec`, and `nnzVec`.
- `nnzElem::DataFrameRow`: Updated DataFrameRow `nnzElem` with the updated element information.
- `updateFlag::Bool`: A flag indicating if the sparse matrix was updated.

## Example
```julia
nnzVec = DataFrame(ID = 1:5, NROW = [1, 1, 2, 2, 3], NCOL = [1, 2, 1, 2, 1], Val = [1.0, 2.0, 3.0, 4.0, 5.0], NIR = [2, -1, 4, -1, -1], NIC = [-1, -1, -1, -1, -1])
nnzElem = DataFrameRow([0, 3, 2, 3, 1], ["ID", "NROW", "NCOL", "Val", "NIR", "NIC"])
sparMat = (NVec = DataFrame(FIR = [1, -1, 3]), MVec = DataFrame(FIC = [-1, -1, -1]), nnzVec = nnzVec)
sparMat, nnzElem, updateFlag = checkElementIntoRow(sparMat, nnzElem, type="replace", verbose=true)
```
"""
function checkElementIntoRow(sparMat::SparseMatrix,
    nnzElem::DataFrameRow;
	type::String="replace", verbose::Bool=false)
	
    NVec = sparMat.NVec
    MVec = sparMat.MVec
    nnzVec = sparMat.nnzVec
	numExistingElems = size(nnzVec, 1)
	myprintln(verbose, "Currently the Sparse Matrix has $numExistingElems elements.")
	nnzElem.ID = numExistingElems + 1
	elID = nnzElem.ID
	FIR = NVec.FIR
	row = nnzElem.NROW
	col = nnzElem.NCOL
	elVal = nnzElem.Val

	updateFlag = false
    stillFindingAPlaceInRow = true

    # Check if our elem is the very first element to be inserted for that row
    if FIR[row] == -1 # First element to be inserted into that row
        myprintln(verbose, "This elem is the first in row $row !")
        FIR[row] = elID
        stillFindingAPlaceInRow = false
    else
        incumbentID = FIR[row]
        myprintln(verbose, "Row $row already has an element, with ID: $incumbentID")
        incumbent = nnzVec[incumbentID, :]
    end

    # If there exists at least one element (the incumbent) already present in the row,
    # check if our elem can come before it.
    if stillFindingAPlaceInRow
        if col < incumbent.NCOL # New elem comes before incumbent in that row => replace it in FIR and shift incumbent in the nnzVec
            myprintln(verbose, "Our elem comes before element $incumbentID, can topple it!")
            FIR[row] = elID
            nnzElem.NIR = incumbentID
            stillFindingAPlaceInRow = false
        elseif col == incumbent.NCOL
            myprintln(verbose, "Stand down. It's a draw between our elem and element $incumbentID.")
            nnzVec = resolveTie(nnzVec, incumbentID, elVal, type=type, verbose=verbose)
            updateFlag = true
            stillFindingAPlaceInRow = false
        else
            myprintln(verbose, "Our elem will come after the incumbent element $incumbentID. Continue Searching.")
            prevIncumbentID = incumbentID
            incumbentID = incumbent.NIR # going to the next element in the row
        end
    end

    # If there is no change to the FIR,
    # Keep checking until our elem finds its place in that row.
    while stillFindingAPlaceInRow
        if incumbentID == -1 # elem is the last in the row
            myprintln(verbose, "Our elem will sit right after the FIR element $prevIncumbentID !")
            nnzVec.NIR[prevIncumbentID] = elID
            stillFindingAPlaceInRow = false
        else # more elements to check in the row
            incumbent = nnzVec[incumbentID, :]
            myprintln(verbose, "Not the second element to be added to this row either. Keep searching.")
        end

        if stillFindingAPlaceInRow
            if col < incumbent.NCOL
                myprintln(verbose, "Our element can topple the incumbent element $incumbentID.")
                nnzElem.NIR = incumbentID
                nnzVec.NIR[prevIncumbentID] = elID
                stillFindingAPlaceInRow = false
            elseif col == incumbent.NCOL # Same-same? Either replace or add
                myprintln(verbose, "Stand down. It's a draw between our elem and element $incumbentID.")
                nnzVec = resolveTie(nnzVec, incumbentID, elVal, type=type, verbose=verbose)
                updateFlag = true
                stillFindingAPlaceInRow = false
            else # col > incumbent.NCOL
                myprintln(verbose, "Not coming before incumbent element $incumbentID. Keep searching.")
                prevIncumbentID = incumbentID
                incumbentID = incumbent.NIR
            end
        end
    end

	NVec.FIR = FIR
    sparMat = (NVec=NVec, MVec=MVec, nnzVec=nnzVec)
	return sparMat, nnzElem, updateFlag
end

"""
    checkElementIntoColumn(sparMat::SparseMatrix,
    nnzElem::DataFrameRow;
	type::String="replace", verbose::Bool=false)

Check and insert an element into a column of a sparse matrix.

This function is analogous to the `checkElementIntoRow` function and serves the same purpose,
but it operates on columns instead of rows. It checks if the given `nnzElem` can be inserted
into the column specified by `nnzElem.NCOL` in the `nnzVec` DataFrame. If the column is empty,
the element is inserted as the first element. If there is already an element present, the function
determines the correct position for the element based on the `nnzElem.NROW` value. The position
can be before the incumbent element, after it, or it may result in a tie, which can be resolved
using the `resolveTie` function.

## Arguments
- `sparMat::SparseMatrix`: A named tuple representing the sparse matrix with elements `NVec`, `MVec`, and `nnzVec`. `NVec` is a DataFrame containing additional information about the sparse matrix, `MVec` is a DataFrame containing additional information about the columns, and `nnzVec` is a DataFrame representing the non-zero elements of the sparse matrix.
- `nnzElem::DataFrameRow`: DataFrameRow representing the new element to be inserted.
- `type::String`: Optional. The type of tie resolution. Default is "replace". Possible values are "replace" and "add".
- `verbose::Bool`: Optional. If set to `true`, print verbose output. Default is `false`.

## Returns
- `sparMat::SparseMatrix`: Updated named tuple `sparMat` with the modified `NVec`, `MVec`, and `nnzVec`.
- `nnzElem::DataFrameRow`: Updated DataFrameRow `nnzElem` with the updated element information.
- `updateFlag::Bool`: A flag indicating if the sparse matrix was updated.

## Example
```julia
nnzVec = DataFrame(ID = 1:5, NROW = [1, 1, 2, 2, 3], NCOL = [1, 2, 1, 2, 1], Val = [1.0, 2.0, 3.0, 4.0, 5.0], NIR = [2, -1, 4, -1, -1], NIC = [-1, -1, -1, -1, -1])
nnzElem = DataFrameRow([0, 3, 2, 3, 1], ["ID", "NROW", "NCOL", "Val", "NIR", "NIC"])
sparMat = (NVec = DataFrame(FIR = [1, -1, 3]), MVec = DataFrame(FIC = [-1, -1, -1]), nnzVec = nnzVec)
sparMat, nnzElem, updateFlag = checkElementIntoColumn(sparMat, nnzElem, type="replace", verbose=true)
```
"""
function checkElementIntoColumn(sparMat::SparseMatrix,
    nnzElem::DataFrameRow;
	type::String="replace", verbose::Bool=false)
	
    NVec = sparMat.NVec
    MVec = sparMat.MVec
    nnzVec = sparMat.nnzVec

	numExistingElems = size(nnzVec, 1)
	myprintln(verbose, "Currently the Sparse Matrix has $numExistingElems elements.")
	nnzElem.ID = numExistingElems + 1
	elID = nnzElem.ID
	FIC = MVec.FIC
	row = nnzElem.NROW
	col = nnzElem.NCOL
	elVal = nnzElem.Val

	updateFlag = false
    stillFindingAPlaceInColumn = true

    # Check if our elem is the very first element to be inserted for that column
    if FIC[col] == -1 # First element to be inserted into that column
        myprintln(verbose, "This elem is the first in column $col !")
        FIC[col] = elID
        stillFindingAPlaceInColumn = false
    else
        incumbentID = FIC[col]
        myprintln(verbose, "Column $col already has an element, with ID: $incumbentID")
        incumbent = nnzVec[incumbentID, :]
    end

    # If there exists at least one element (the incumbent) already present in the column,
    # check if our elem can come before it.
    if stillFindingAPlaceInColumn
        if row < incumbent.NROW # New elem comes before incumbent in that column => replace it in FIC and shift incumbent in the nnzVec
            myprintln(verbose, "Our elem comes before element $incumbentID, can topple it!")
            FIC[col] = elID
            nnzElem.NIC = incumbentID
            stillFindingAPlaceInColumn = false
        elseif row == incumbent.NROW
            myprintln(verbose, "Stand down. It's a draw between our elem and element $incumbentID.")
            nnzVec = resolveTie(nnzVec, incumbentID, elVal, type=type, verbose=verbose)
            updateFlag = true
            stillFindingAPlaceInColumn = false
        else
            myprintln(verbose, "Our elem will come after the incumbent element $incumbentID. Continue Searching.")
            prevIncumbentID = incumbentID
            incumbentID = incumbent.NIC # going to the next element in the column
        end
    end

    # If there is no change to the FIC,
    # Keep checking until our elem finds its place in that column.
    while stillFindingAPlaceInColumn
        if incumbentID == -1 # elem is the last in the column
            myprintln(verbose, "Our elem will sit right after the FIC element $prevIncumbentID !")
            nnzVec.NIC[prevIncumbentID] = elID
            stillFindingAPlaceInColumn = false
        else # more elements to check in the column
            incumbent = nnzVec[incumbentID, :]
            myprintln(verbose, "Not the second element to be added to this column either. Keep searching.")
        end

        if stillFindingAPlaceInColumn
            if row < incumbent.NROW
                myprintln(verbose, "Our element can topple the incumbent element $incumbentID.")
                nnzElem.NIC = incumbentID
                nnzVec.NIC[prevIncumbentID] = elID
                stillFindingAPlaceInColumn = false
            elseif row == incumbent.NROW # Same-same? Either replace or add
                myprintln(verbose, "Stand down. It's a draw between our elem and element $incumbentID.")
                nnzVec = resolveTie(nnzVec, incumbentID, elVal, type=type, verbose=verbose)
                updateFlag = true
                stillFindingAPlaceInColumn = false


            else # row > incumbent.NROW
                myprintln(verbose, "Not coming before incumbent element $incumbentID. Keep searching.")
                prevIncumbentID = incumbentID
                incumbentID = incumbent.NIC
            end
        end
    end

	MVec.FIC = FIC
    sparMat = (NVec=NVec, MVec=MVec, nnzVec=nnzVec)
	return sparMat, nnzElem, updateFlag
end


"""
    updateSparse(sparMat::SparseMatrix,
    nnzElem::DataFrameRow;
	type::String = "replace", verbose::Bool = false)

Update a sparse matrix by inserting a new element.

This function updates the sparse matrix represented by `nnzVec` and additional information
in `NVec` by inserting a new element specified by the `nnzElem` DataFrameRow. The element is
inserted by calling the `checkElementIntoRow` and `checkElementIntoColumn` functions, which
determine the correct position for the element based on its row and column indices.

If the specified position for the new element is already occupied by a previous "incumbent"
element, the `type` argument determines how the tie is resolved. If `type` is set to "replace",
the new element's value replaces the incumbent element's value. If `type` is set to "add",
the new element's value is added to the incumbent element's value.

## Arguments
- `sparMat::SparseMatrix`: A named tuple representing the sparse matrix with elements `NVec`, `MVec`, and `nnzVec`. `NVec` is a DataFrame containing additional information about the sparse matrix, `MVec` is a DataFrame containing additional information about the columns, and `nnzVec` is a DataFrame representing the non-zero elements of the sparse matrix.
- `nnzElem::DataFrameRow`: DataFrameRow representing the new element to be inserted.
- `type::String`: Optional. The type of tie resolution. Default is "replace". Possible values are "replace" and "add".
- `verbose::Bool`: Optional. If set to `true`, print verbose output. Default is `false`.

## Returns
- `sparMat::SparseMatrix`: Updated named tuple `sparMat` with the modified `NVec`, `MVec`, and `nnzVec`.

## Example
```julia
NVec = DataFrame(FIR = [-1, -1, -1], FIC = [-1, -1, -1])
nnzVec = DataFrame(ID = 1:5, NROW = [1, 1, 2, 2, 3], NCOL = [1, 2, 1, 2, 1], Val = [1.0, 2.0, 3.0, 4.0, 5.0], NIR = [2, -1, 4, -1, -1], NIC = [-1, -1, -1, -1, -1])
nnzElem = DataFrameRow([0, 3, 2, 3, 1], ["ID", "NROW", "NCOL", "Val", "NIR", "NIC"])
sparMat = (NVec = NVec, MVec = DataFrame(FIC = [-1, -1, -1]), nnzVec = nnzVec)
sparMat = updateSparse(sparMat, nnzElem, type="replace", verbose=true)
```
"""
function updateSparse(sparMat::SparseMatrix,
    nnzElem::DataFrameRow;
	type::String = "replace",
	verbose::Bool = false)

	sparMat, nnzElem, updateFlag = checkElementIntoRow(sparMat, nnzElem, type=type, verbose=verbose)
    sparMat, nnzElem, updateFlag = checkElementIntoColumn(sparMat, nnzElem, type="maintain", verbose=verbose)

	if updateFlag == false
		myprintln(verbose, "Another element to be added to the sparse matrix.")
		myprintln(verbose, "The element:")
		myprintln(verbose, nnzElem)
		myprintln(verbose, "nnzVec before:")
		myprintln(verbose, sparMat.nnzVec)
		push!(sparMat.nnzVec, nnzElem)
		myprintln(verbose, "nnzVec after:")
		myprintln(verbose, sparMat.nnzVec)
	else
		myprintln(verbose, "An element was updated, but not added to the sparse matrix.")
	end

	return sparMat
end

function sparseMatrixConstructor(N::Int64, M::Int64; T::DataType=ComplexF64)
    (firs, fics) = (repeat([-1], N), repeat([-1], M))

    NVec = DataFrame(FIR = firs)
    MVec = DataFrame(FIC = fics)
    nnzVec = DataFrame(ID = Int64[], Val = T[], NROW = Int64[], NCOL = Int64[], NIR = Int64[], NIC = Int64[])
    sparMat = (NVec=NVec, MVec=MVec, nnzVec=nnzVec)

    return sparMat
end

"""
    constructSparseYBus(CDF_DF_List_pu::Vector{DataFrame};
    disableTaps::Bool = false,
    sortBy::String = "busNumbers",
    verbose::Bool = false,
    saveTables::Bool = false,
    saveLocation::String = "processedData/")

This function constructs a sparse YBus representation from the provided `busData` and `branchData` extracted from `CDF_DF_List_pu`.

## Arguments:
- `CDF_DF_List_pu::Vector{DataFrame}`: A vector of DataFrames containing `busData` and `branchData`.
- `disableTaps::Bool`: A boolean flag indicating whether to disable taps. Default is `false`.
- `sortBy::String`: A string indicating the sorting order. Default is `"busNumbers"`.
- `verbose::Bool`: A boolean flag indicating whether to print verbose output. Default is `false`.
- `saveTables::Bool`: A boolean flag indicating whether to save tables. Default is `false`.
- `saveLocation::String`: A string indicating the location to save the tables. Default is `"processedData/".

## Returns:
- `NVec::DataFrame`: A `DataFrame` representing the N vector of the YBus matrix.
- `nnzVec::DataFrame`: A `DataFrame` representing the non-zero values of the YBus matrix.
"""
function constructSparseYBus(CDF_DF_List_pu::Vector{DataFrame};
    disableTaps::Bool = false,
    sortBy::String = "busNumbers",
    verbose::Bool = false,
    saveTables::Bool = false,
    saveLocation::String = "processedData/")

    systemName = extractSystemName(CDF_DF_List_pu)
    busData_pu = CDF_DF_List_pu[2]
    branchData_pu = CDF_DF_List_pu[3]
    N = size(busData_pu, 1)
    numBranch = size(branchData_pu, 1)

    sparMat = sparseMatrixConstructor(N, N)

    for branch = 1:numBranch
        currentBranch = branchData_pu[branch, :]
        i = currentBranch.Tap_Bus_Num
        k = currentBranch.Z_Bus_Num

        if disableTaps
            a = 1
        elseif currentBranch.Transformer_t != 0
            a = currentBranch.Transformer_t
        else
            a = 1
        end

        y_ik = 1/(currentBranch.R_pu + im*currentBranch.X_pu)

        ybus_ii = y_ik/(a^2) + im*currentBranch.B_pu / 2
        ybus_kk = y_ik + im*currentBranch.B_pu / 2
        ybus_ik = -y_ik/a
        ybus_ki = -y_ik/a

        cases = ["ii", "kk", "ik", "ki"]

        for case in cases
            if case == "ii"
                compElem = DataFrame(ID = -1, Val = ybus_ii, i = i, j = i)
            elseif case == "kk"
                compElem = DataFrame(ID = -1, Val = ybus_kk, i = k, j = k)
            elseif case == "ik"
                compElem = DataFrame(ID = -1, Val = ybus_ik, i = i, j = k)
            elseif case == "ki"
                compElem = DataFrame(ID = -1, Val = ybus_ki, i = k, j = i)
            end

            compElem = compElem[1, :]
            nnzElem = nnzRowConstructor(compElem)
            sparMat = updateSparse(sparMat, nnzElem, type="add", verbose=false)
        end

    end

    for bus = 1:N
        ybus_ii = busData_pu.G_pu[bus] + im*busData_pu.B_pu[bus]
        compElem = DataFrame(ID = -1, Val = ybus_ii, i = bus, j = bus)
        compElem = compElem[1, :]
        nnzElem = nnzRowConstructor(compElem)
        sparMat = updateSparse(sparMat, nnzElem, type="add", verbose=false)
    end

    return sparMat
end

"""
    computeMismatchesViaSparseYBus(PSpecified::Vector{Float64},
    QSpecified::Vector{Float64},
    V::Vector{Float64},
    delta::Vector{Float64},
    sparMat::SparseMatrix) -> Tuple{Vector{Float64}, Vector{Float64}}

Calculate mismatches between specified and computed power values using the sparse YBus matrix.

## Arguments
- `PSpecified::Vector{Float64}`: Vector of specified active power values for each bus.
- `QSpecified::Vector{Float64}`: Vector of specified reactive power values for each bus.
- `V::Vector{Float64}`: Vector of voltage magnitudes for each bus.
- `delta::Vector{Float64}`: Vector of voltage phase angles for each bus.
- `sparMat::SparseMatrix`: Named tuple representing the sparse YBus matrix, with fields `NVec`, `MVec`, and `nnzVec`.

## Returns
- `deltaP::Vector{Float64}`: Vector representing the mismatches between specified and computed active power values.
- `deltaQ::Vector{Float64}`: Vector representing the mismatches between specified and computed reactive power values.

## Details
The `computeMismatchesViaSparseYBus` function calculates the mismatches between the specified active and reactive power values (`PSpecified` and `QSpecified`) and the computed power values using the sparse YBus matrix represented by `sparMat`. It takes the voltage magnitudes (`V`) and voltage phase angles (`delta`) as inputs.

The function iterates over each bus in the system and calculates the active and reactive power mismatches using the specified power values, voltage magnitudes, voltage phase angles, and the non-zero elements of the YBus matrix.

The active power mismatch `deltaP` is computed as the difference between the specified active power and the computed active power for each bus.

The reactive power mismatch `deltaQ` is computed as the difference between the specified reactive power and the computed reactive power for each bus.

Note that the `sparMat` named tuple should contain the necessary dataframes representing the non-zero elements of the YBus matrix.

```julia
deltaP, deltaQ = computeMismatchesViaSparseYBus(PSpecified, QSpecified, V, delta, sparMat)
```
"""
function computeMismatchesViaSparseYBus(PSpecified::Vector{Float64}, 
    QSpecified::Vector{Float64}, 
    V::Vector{Float64}, 
    delta::Vector{Float64}, 
    sparMat::SparseMatrix)

    NYBus = sparMat.NVec
    MYBus = sparMat.MVec
    nnzYBus = sparMat.nnzVec
    
    N = length(NYBus.FIR)
    (P, Q) = (repeat([0.0000], N), repeat([0.0000], N))

    for bus = 1:N
        i = bus
        elemkNum = NYBus.FIR[bus]
        while elemkNum != -1
            elemk = nnzYBus[elemkNum, :]
            k = elemk.NCOL
            YBus_ik = elemk.Val
            P[i] += abs( YBus_ik*V[i]*V[k] ) * cos( angle(YBus_ik) + delta[k] - delta[i])
            Q[i] += -abs( YBus_ik*V[i]*V[k]) * sin( angle(YBus_ik) + delta[k] - delta[i])
            elemkNum = elemk.NIR
        end
    end

    deltaP = PSpecified - P
    deltaQ = QSpecified - Q
    return deltaP, deltaQ
end

function constructSparseJacobian(CDF_DF_List_pu::Vector{DataFrame},
    P::Vector{Float64},
    Q::Vector{Float64},
    V::Vector{Float64},
    delta::Vector{Float64},
    YBus::SparseMatrix;
    verbose::Bool=false,
    combinationOrder::String="hcat-then-vcat",
    saveTable::Bool=false,
    processedDataFolder::String="processedData/")

    powSysData = initializeVectors_pu(CDF_DF_List_pu)
    nPV = powSysData.nPV
    nPQ = powSysData.nPQ
    lPV = powSysData.listOfPVBuses
    lPQ = powSysData.listOfPQBuses

    J11 = constructSparseJacobianSubMatrix(CDF_DF_List_pu, P, Q, V, delta, YBus, type="J11")
    J12 = constructSparseJacobianSubMatrix(CDF_DF_List_pu, P, Q, V, delta, YBus, type="J12")
    J21 = constructSparseJacobianSubMatrix(CDF_DF_List_pu, P, Q, V, delta, YBus, type="J21")
    J22 = constructSparseJacobianSubMatrix(CDF_DF_List_pu, P, Q, V, delta, YBus, type="J22")

    if combinationOrder == "hcat-then-vcat"
        JTop = hcatSparse(J11, J12)
        JBottom = hcatSparse(J21, J22)
        J = vcatSparse(JTop, JBottom)
    elseif combinationOrder == "vcat-then-hcat"
        JLeft = vcatSparse(J11, J21)
        JRight = vcatSparse(J12, J22)
        J = hcatSparse(JLeft, JRight)
    else
        error("Unknown combination order.")
    end

    return J
end

function constructSparseJacobianSubMatrix(CDF_DF_List_pu::Vector{DataFrame},
    P::Vector{Float64},
    Q::Vector{Float64},
    V::Vector{Float64},
    delta::Vector{Float64},
    sparYBus::SparseMatrix;
    type::String="J11",
    verbose::Bool=false)

    busData_pu = CDF_DF_List_pu[2]
    N = size(busData_pu, 1)

    busPositions = busPositionsForJacobian(CDF_DF_List_pu)

    NYBus = sparYBus.NVec
    MYBus = sparYBus.MVec
    nnzYBus = sparYBus.nnzVec

    powSysData = initializeVectors_pu(CDF_DF_List_pu)
    nPV = powSysData.nPV
    nPQ = powSysData.nPQ
    lPV = powSysData.listOfPVBuses
    lPQ = powSysData.listOfPQBuses
    lNonSlack = powSysData.listOfNonSlackBuses

    if type == "J11"
        JSub = sparseMatrixConstructor(N-1, N-1, T=Float64)
    elseif type == "J12"
        JSub = sparseMatrixConstructor(N-1, nPQ, T=Float64)
    elseif type == "J21"
        JSub = sparseMatrixConstructor(nPQ, N-1, T=Float64)
    elseif type == "J22"
        JSub = sparseMatrixConstructor(nPQ, nPQ, T=Float64)
    else
        error("Unknown Jacobian Sub-matrix.")
    end

    for i = 1:N
        elemNum = NYBus.FIR[i]
        while elemNum != -1
            elem = nnzYBus[elemNum, :]
            k = elem.NCOL

            if type == "J11"
                row = busPositions.nonSlackBusIdx[i]
                col = busPositions.nonSlackBusIdx[k]
            elseif type == "J12"
                row = busPositions.nonSlackBusIdx[i]
                col = busPositions.PQBusIdx[k]
            elseif type == "J21"
                row = busPositions.PQBusIdx[i]
                col = busPositions.nonSlackBusIdx[k]
            elseif type == "J22"
                row = busPositions.PQBusIdx[i]
                col = busPositions.PQBusIdx[k]
            else
                error("Unknown Jacobian Submatrix type")
            end

            if row != -1 && col != -1

                Y_ik = elem.Val
                if i == k
                    B_ii = imag(Y_ik)
                    G_ii = real(Y_ik)
                    Vi_squared = V[i]^2
                    j11 = -Q[i] - B_ii*Vi_squared
                    j12 = P[i] + G_ii*Vi_squared
                    j21 = P[i] - G_ii*Vi_squared
                    j22 = Q[i] - B_ii*Vi_squared
                else
                    YVkVi = abs(Y_ik)*V[k]*V[i]
                    gamma_deltak_delta_i = angle(Y_ik) + delta[k] - delta[i]
                    j11 = -YVkVi*sin(gamma_deltak_delta_i)
                    j12 = YVkVi*cos(gamma_deltak_delta_i)
                    j21 = -j12
                    j22 = j11
                end
                
                if type == "J11"
                    compElem = DataFrame(i = row, j = col, Val = j11)
                elseif type == "J12"
                    compElem = DataFrame(i = row, j = col, Val = j12)
                elseif type == "J21"
                    compElem = DataFrame(i = row, j = col, Val = j21)
                elseif type == "J22"
                    compElem = DataFrame(i = row, j = col, Val = j22)
                else
                    error("Unknown Jacobian Submatrix type.")
                end

                nnzElem = nnzRowConstructor(compElem[1, :])
                JSub = updateSparse(JSub, nnzElem)
            end
            
            elemNum = elem.NIR;
        end
    end

    return JSub
end

function busPositionsForJacobian(CDF_DF_List_pu::Vector{DataFrame};
    verbose::Bool=true)

    busPositions = DataFrame(busNum = Int64[], nonSlackBusIdx = Int64[], PQBusIdx = Int64[])
    counterPQ = 0
    counterPV = 0

    busData = CDF_DF_List_pu[2]
    N = size(busData, 1)

    for busNum = 1:N
        busType = busData.Type[busNum]
        if busType == 0
            counterPQ = counterPQ + 1
            push!(busPositions.busNum, busNum)
            push!(busPositions.nonSlackBusIdx, counterPQ + counterPV)
            push!(busPositions.PQBusIdx, counterPQ)
        elseif busType == 2
            counterPV = counterPV + 1
            push!(busPositions.busNum, busNum)
            push!(busPositions.nonSlackBusIdx, counterPQ + counterPV)
            push!(busPositions.PQBusIdx, -1)
        elseif busType == 3
            push!(busPositions.busNum, busNum)
            push!(busPositions.nonSlackBusIdx, -1)
            push!(busPositions.PQBusIdx, -1)
        else
            error("Not prepared for this bus type.")
        end
    end

    return busPositions
end

"""
    hcatSparse(matLeft::SparseMatrix, 
        matRight::SparseMatrix;
        verbose::Bool = false)

The `hcatSparse` function performs horizontal concatenation of two sparse matrices represented by named tuples `matLeft` and `matRight`. It concatenates the `matLeft` and `matRight` DataFrames horizontally and returns a new sparse matrix represented by the named tuple `matHorz`.

Arguments:
- `matLeft::SparseMatrix`: A named tuple containing three DataFrames representing the left sparse matrix with columns `NVec`, `MVec`, and `nnzVec`.
- `matRight::SparseMatrix`: A named tuple containing three DataFrames representing the right sparse matrix with columns `NVec`, `MVec`, and `nnzVec`.
- `verbose::Bool`: (optional) A Boolean value indicating whether to display verbose output during the concatenation process. Default is `false`.

Returns:
- A named tuple `matHorz` containing three DataFrames `NVec`, `MVec`, and `nnzVec` representing the horizontally concatenated sparse matrix.

Example:
```julia
NLeft1 = DataFrame(FIR = [1, 2, 4])
MLeft1 = DataFrame(FIC = [2, 1, 3])
nnzLeft1 = DataFrame(ID = [1, 2, 3, 4, 5], NROW = [1, 2, 2, 3, 3], NCOL = [2, 1, 3, 2, 3], NIR = [-1, 3, -1, 5, -1], NIC = [4, -1, 5, -1, -1], Val = [1, 1, 1, 1, 1])
matLeft = (NVec=NLeft1, MVec=MLeft1, nnzVec=nnzLeft1)

NRight1 = DataFrame(FIR = [1, 2, 3])
MRight1 = DataFrame(FIC = [1, 2])
nnzRight1 = DataFrame(ID = [1, 2, 3], NROW = [1, 2, 3], NCOL = [1, 2, 1], NIR = [-1, -1, -1], NIC = [3, -1, -1], Val = [1, 1, 1])
matRight = (NVec=NRight1, MVec=MRight1, nnzVec=nnzRight1)

matHorz = hcatSparse(matLeft, matRight)
"""
function hcatSparse(matLeft::SparseMatrix, 
    matRight::SparseMatrix;
    verbose::Bool = false)

    if length(matRight.NVec.FIR) != length(matLeft.NVec.FIR)
        error("Horizontal Concatenation NOT possible for mismatching number of rows.")
    end

    matRightMod = modify_matRight(matLeft, matRight)
    NRightVec, MRightVec, nnzRightVec = matRightMod.NVec, matRightMod.MVec, matRightMod.nnzVec

    matLeftMod = modify_matLeft(matLeft, matRightMod)
    NLeftVec, MLeftVec, nnzLeftVec = matLeftMod.NVec, matLeftMod.MVec, matLeftMod.nnzVec


    NVec = NLeftVec
    MVec = DataFrame(FIC = vcat(MLeftVec.FIC, MRightVec.FIC))
    nnzVec = vcat(nnzLeftVec, nnzRightVec)
    matHorz = (NVec=NVec, MVec=MVec, nnzVec=nnzVec)

    return matHorz
end

"""
    vcatSparse(matTop::SparseMatrix, 
        matBottom::SparseMatrix;
        verbose::Bool = false)

The `vcatSparse` function performs vertical concatenation of two sparse matrices represented by named tuples `matTop` and `matBottom`. It concatenates the `matTop` and `matBottom` DataFrames vertically and returns a new sparse matrix represented by the named tuple `matVert`.

Arguments:
- `matTop::SparseMatrix`: A named tuple containing three DataFrames representing the top sparse matrix with columns `NVec`, `MVec`, and `nnzVec`.
- `matBottom::SparseMatrix`: A named tuple containing three DataFrames representing the bottom sparse matrix with columns `NVec`, `MVec`, and `nnzVec`.
- `verbose::Bool`: (optional) A Boolean value indicating whether to display verbose output during the concatenation process. Default is `false`.

Returns:
- A named tuple `matVert` containing three DataFrames `NVec`, `MVec`, and `nnzVec` representing the vertically concatenated sparse matrix.

Example:
```julia
matTop = (NVec=DataFrame(FIR=[1, 3, 4]), MVec=DataFrame(FIC=[2, 1, 3]), nnzVec=DataFrame(ID=[1, 2, 3, 4, 5], NROW=[1, 1, 2, 3, 3], NCOL=[1, 3, 2, 2, 3], NIR=[-1, 3, -1, 5, -1], NIC=[4, -1, 5, -1, -1], Val=[11, 22, 33, 44, 55]))
matBottom = (NVec=DataFrame(FIR=[1]), MVec=DataFrame(FIC=[1, 2, 3]), nnzVec=DataFrame(ID=[1, 2, 3], NROW=[1, 2, 3], NCOL=[1, 2, 1], NIR=[-1, -1, -1], NIC=[3, -1, -1], Val=[66, 77, 88]))

matVert = vcatSparse(matTop, matBottom, verbose=true)
"""
function vcatSparse(matTop::SparseMatrix, 
    matBottom::SparseMatrix;
    verbose::Bool = false)

    if length(matTop.MVec.FIC) != length(matBottom.MVec.FIC)
        error("Vertical Concatenation NOT possible for mismatching number of columns.")
    end

    matBottomMod = modify_matBottom(matTop, matBottom)
    NBottomVec, MBottomVec, nnzBottomVec = matBottomMod.NVec, matBottomMod.MVec, matBottomMod.nnzVec

    matTopMod = modify_matTop(matTop, matBottomMod)
    NTopVec, MTopVec, nnzTopVec = matTopMod.NVec, matTopMod.MVec, matTopMod.nnzVec

    NVec = DataFrame(FIR = vcat(NTopVec.FIR, NBottomVec.FIR))
    MVec = MTopVec
    nnzVec = vcat(nnzTopVec, nnzBottomVec)
    matVert = (NVec=NVec, MVec=MVec, nnzVec=nnzVec)

    return matVert
end

"""
    modify_matRight(matLeft::SparseMatrix,
        matRight::SparseMatrix)

The `modify_matRight` function modifies the `matRight` named tuple to prepare it for horizontal concatenation with `matLeft`. It updates the column indices of the elements in the right sparse matrix and ensures that the matrix can be properly concatenated with the left sparse matrix.

Arguments:
- `matLeft::SparseMatrix`: A named tuple containing three DataFrames representing the left sparse matrix with columns `NVec`, `MVec`, and `nnzVec`.
- `matRight::SparseMatrix`: A named tuple containing three DataFrames representing the right sparse matrix with columns `NVec`, `MVec`, and `nnzVec`.

Returns:
- A modified named tuple `matRightCopy` with updated column indices to prepare it for horizontal concatenation with `matLeft`.

Example:
```julia
matLeft = (NVec=DataFrame(FIR=[1, 3, 4]), MVec=DataFrame(FIC=[2, 1, 3]), nnzVec=DataFrame(ID=[1, 2, 3, 4, 5], NROW=[1, 1, 2, 3, 3], NCOL=[1, 3, 2, 2, 3], NIR=[-1, 3, -1, 5, -1], NIC=[4, -1, 5, -1, -1], Val=[11, 22, 33, 44, 55]))
matRight = (NVec=DataFrame(FIR=[1, 2, 3]), MVec=DataFrame(FIC=[1, 2]), nnzVec=DataFrame(ID=[1, 2, 3], NROW=[1, 2, 3], NCOL=[1, 2, 1], NIR=[-1, -1, -1], NIC=[3, -1, -1], Val=[66, 77, 88]))

matRightCopy = modify_matRight(matLeft, matRight)
"""
function modify_matRight(matLeft::SparseMatrix,
    matRight::SparseMatrix)
    
    NLeftVec, MLeftVec, nnzLeftVec = matLeft.NVec, matLeft.MVec, matLeft.nnzVec
    N = length(NLeftVec.FIR)
    MLeft = length(MLeftVec.FIC)
    nnzLeft = length(nnzLeftVec.ID)

    matRightCopy = deepcopy(matRight)
    NRightVec, MRightVec, nnzRightVec = matRightCopy.NVec, matRightCopy.MVec, matRightCopy.nnzVec

    for row in 1:N
        if NRightVec.FIR[row] != -1
            NRightVec.FIR[row] += nnzLeft
        end
    end

    MRight = length(MRightVec.FIC)
    for col in 1:MRight
        if MRightVec.FIC[col] != -1
            MRightVec.FIC[col] += nnzLeft
        end
    end
    
    nnzRightVec.ID .+= nnzLeft
    nnzRightVec.NCOL .+= MLeft
    nnzRight = length(nnzRightVec.ID)
    for elemNum in 1:nnzRight
        if nnzRightVec.NIR[elemNum] != -1
            nnzRightVec.NIR[elemNum] += nnzLeft
        end
        if nnzRightVec.NIC[elemNum] != -1
            nnzRightVec.NIC[elemNum] += nnzLeft
        end
    end

    matRightCopy = (NVec=NRightVec, MVec=MRightVec, nnzVec=nnzRightVec)
    return matRightCopy
end

"""
    modify_matBottom(matTop::SparseMatrix,
        matBottom::SparseMatrix)

The `modify_matBottom` function modifies the `matBottom` named tuple to prepare it for vertical concatenation with `matTop`. It updates the row indices of the elements in the bottom sparse matrix to properly align them with the top sparse matrix.

Arguments:
- `matTop::SparseMatrix`: A named tuple containing three DataFrames representing the top sparse matrix with columns `NVec`, `MVec`, and `nnzVec`.
- `matBottom::SparseMatrix`: A named tuple containing three DataFrames representing the bottom sparse matrix with columns `NVec`, `MVec`, and `nnzVec`.

Returns:
- A modified named tuple `matBottomCopy` with updated row indices to prepare it for vertical concatenation with `matTop`.

Example:
```julia
matTop = (NVec=DataFrame(FIR=[1, 3, 4]), MVec=DataFrame(FIC=[2, 1, 3]), nnzVec=DataFrame(ID=[1, 2, 3, 4, 5], NROW=[1, 1, 2, 3, 3], NCOL=[1, 3, 2, 2, 3], NIR=[-1, 3, -1, 5, -1], NIC=[4, -1, 5, -1, -1], Val=[11, 22, 33, 44, 55]))
matBottom = (NVec=DataFrame(FIR=[1, 2, 3]), MVec=DataFrame(FIC=[1, 2]), nnzVec=DataFrame(ID=[1, 2, 3], NROW=[1, 2, 3], NCOL=[1, 2, 1], NIR=[-1, -1, -1], NIC=[3, -1, -1], Val=[66, 77, 88]))

matBottomCopy = modify_matBottom(matTop, matBottom)
"""
function modify_matBottom(matTop::SparseMatrix,
    matBottom::SparseMatrix)
    
    NTopVec, MTopVec, nnzTopVec = matTop.NVec, matTop.MVec, matTop.nnzVec
    NTop = length(NTopVec.FIR)
    M = length(MTopVec.FIC)
    nnzTop = length(nnzTopVec.ID)

    matBottomCopy = deepcopy(matBottom)
    NBottomVec, MBottomVec, nnzBottomVec = matBottomCopy.NVec, matBottomCopy.MVec, matBottomCopy.nnzVec

    for col in 1:M
        if MBottomVec.FIC[col] != -1
            MBottomVec.FIC[col] += nnzTop
        end
    end

    NBottom = length(NBottomVec.FIR)
    for row in 1:NBottom
        if NBottomVec.FIR[row] != -1
            NBottomVec.FIR[row] += nnzTop
        end
    end
    
    nnzBottomVec.ID .+= nnzTop
    nnzBottomVec.NROW .+= NTop
    nnzBottom = length(nnzBottomVec.ID)
    for elemNum in 1:nnzBottom
        if nnzBottomVec.NIC[elemNum] != -1
            nnzBottomVec.NIC[elemNum] += nnzTop
        end
        if nnzBottomVec.NIR[elemNum] != -1
            nnzBottomVec.NIR[elemNum] += nnzTop
        end
    end

    matBottomCopy = (NVec=NBottomVec, MVec=MBottomVec, nnzVec=nnzBottomVec)
    return matBottomCopy
end

"""
    modify_matLeft(matLeft::SparseMatrix,
        matRight::SparseMatrix)

The `modify_matLeft` function modifies the `matLeft` named tuple to prepare it for vertical concatenation with `matRight`. It updates the row indices of the elements in the left sparse matrix to properly align them with the right sparse matrix.

Arguments:
- `matLeft::SparseMatrix`: A named tuple containing three DataFrames representing the left sparse matrix with columns `NVec`, `MVec`, and `nnzVec`.
- `matRight::SparseMatrix`: A named tuple containing three DataFrames representing the right sparse matrix with columns `NVec`, `MVec`, and `nnzVec`.

Returns:
- A modified named tuple `matLeftCopy` with updated row indices to prepare it for vertical concatenation with `matRight`.

Example:
```julia
matLeft = (NVec=DataFrame(FIR=[1, 3, 4]), MVec=DataFrame(FIC=[2, 1, 3]), nnzVec=DataFrame(ID=[1, 2, 3, 4, 5], NROW=[1, 1, 2, 3, 3], NCOL=[1, 3, 2, 2, 3], NIR=[-1, 3, -1, 5, -1], NIC=[4, -1, 5, -1, -1], Val=[11, 22, 33, 44, 55]))
matRight = (NVec=DataFrame(FIR=[1, 2, 3]), MVec=DataFrame(FIC=[1, 2]), nnzVec=DataFrame(ID=[1, 2, 3], NROW=[1, 2, 3], NCOL=[1, 2, 1], NIR=[-1, -1, -1], NIC=[3, -1, -1], Val=[66, 77, 88]))

matLeftCopy = modify_matLeft(matLeft, matRight)
"""
function modify_matLeft(matLeft::SparseMatrix,
    matRight::SparseMatrix)

    matLeftCopy = deepcopy(matLeft)
    NLeftVec, MLeftVec, nnzLeftVec = matLeftCopy.NVec, matLeftCopy.MVec, matLeftCopy.nnzVec
    NRightVec, MRightVec, nnzRightVec = matRight.NVec, matRight.MVec, matRight.nnzVec

    N = length(NLeftVec.FIR)
    for row = 1:N
        lastKnownElemInRow = NLeftVec.FIR[row]
        if lastKnownElemInRow == -1
            NLeftVec.FIR[row] = NRightVec.ID[row]
        else
            nextElemInRow = lastKnownElemInRow
            while nextElemInRow != -1
                lastKnownElemInRow = nextElemInRow
                nextElemInRow = nnzLeftVec.NIR[lastKnownElemInRow]
            end
            nnzLeftVec.NIR[lastKnownElemInRow] = NRightVec.FIR[row]
        end
    end

    matLeftCopy = (NVec=NLeftVec, MVec=MLeftVec, nnzVec=nnzLeftVec)
    return matLeftCopy
end

"""
    modify_matTop(matTop::SparseMatrix,
        matBottom::SparseMatrix;
        verbose::Bool=false)

The `modify_matTop` function modifies the `matTop` named tuple to prepare it for vertical concatenation with `matBottom`. It updates the column indices of the elements in the top sparse matrix to properly align them with the bottom sparse matrix.

Arguments:
- `matTop::SparseMatrix`: A named tuple containing three DataFrames representing the top sparse matrix with columns `NVec`, `MVec`, and `nnzVec`.
- `matBottom::SparseMatrix`: A named tuple containing three DataFrames representing the bottom sparse matrix with columns `NVec`, `MVec`, and `nnzVec`.
- `verbose::Bool`: (optional) A Boolean value indicating whether to display verbose output during the modification process. Default is `false`.

Returns:
- A modified named tuple `matTopCopy` with updated column indices to prepare it for vertical concatenation with `matBottom`.

Example:
```julia
matTop = (NVec=DataFrame(FIR=[1, 3, 4]), MVec=DataFrame(FIC=[2, 1, 3]), nnzVec=DataFrame(ID=[1, 2, 3, 4, 5], NROW=[1, 1, 2, 3, 3], NCOL=[1, 3, 2, 2, 3], NIR=[-1, 3, -1, 5, -1], NIC=[4, -1, 5, -1, -1], Val=[11, 22, 33, 44, 55]))
matBottom = (NVec=DataFrame(FIR=[1, 2, 3]), MVec=DataFrame(FIC=[1, 2]), nnzVec=DataFrame(ID=[1, 2, 3], NROW=[1, 2, 3], NCOL=[1, 2, 1], NIR=[-1, -1, -1], NIC=[3, -1, -1], Val=[66, 77, 88]))

matTopCopy = modify_matTop(matTop, matBottom)
"""
function modify_matTop(matTop::SparseMatrix,
    matBottom::SparseMatrix;
    verbose::Bool=false)

    matTopCopy = deepcopy(matTop)
    NTopVec, MTopVec, nnzTopVec = matTopCopy.NVec, matTopCopy.MVec, matTopCopy.nnzVec
    NBottomVec, MBottomVec, nnzBottomVec = matBottom.NVec, matBottom.MVec, matBottom.nnzVec

    M = length(MTopVec.FIC)
    for col = 1:M
        lastKnownElemInCol = MTopVec.FIC[col]
        if lastKnownElemInCol == -1
            myprintln(verbose, "Column $col in top matrix was empty!?")
            myprintln(verbose, "So, the next element to be connected would be
            $(MBottom.ID[col])")
            MTopVec.FIC[col] = MBottomVec.FIC[col]
        else
            incumbentID = MTopVec.FIC[col]
            myprintln(verbose, "Column $col already has an element, with ID: $incumbentID")
            myprintln(verbose, "So we find the next element.")
            nextElemInCol = lastKnownElemInCol
            while nextElemInCol != -1
                lastKnownElemInCol = nextElemInCol
                nextElemInCol = nnzTopVec.NIC[lastKnownElemInCol]
            end
            nnzTopVec.NIC[lastKnownElemInCol] = MBottomVec.FIC[col]
        end
    end

    matTopCopy = (NVec=NTopVec, MVec=MTopVec, nnzVec=nnzTopVec)
    return matTopCopy
end

"""
    full2comp(matFull::Matrix{T}) where T

The `full2comp` function converts a full matrix `matFull` into a compressed matrix `compMat`. It iterates over the elements of the full matrix, constructs `compMat` by storing the non-zero elements with their row and column indices.

Arguments:
- `matFull::Matrix{T}`: The full matrix to be converted into a compressed matrix, where `T` represents the element type.

Returns:
- A DataFrame `compMat` representing the compressed matrix with columns `i`, `j`, and `Val` containing the row indices, column indices, and values of the non-zero elements, respectively.

Example:
```julia
matFull = [11 0 22; 0 33 0; 0 44 55]
compMat = full2comp(matFull)
"""
function full2comp(matFull::Matrix{T}) where T
    m, n = size(matFull)
    i, j, Val = Int[], Int[], T[]

    for row in 1:m
        for col in 1:n
            if matFull[row, col] != zero(T) # Use zero(T) to handle different numeric types
                push!(i, row)
                push!(j, col)
                push!(Val, matFull[row, col])
            end
        end
    end

    compMat = DataFrame(i=i, j=j, Val=Val)
    return compMat
end

function computeSparsity(A::SparseMatrix;
    returnValue::String="print only",
    verbose::Bool = false)

    N, M, nnz = length(A.NVec.FIR), length(A.MVec.FIC), length(A.nnzVec.ID)

    if N == 0 || M == 0
        error("Divide by zero error. This is a non-empty 2D array, right?")
    else
        sparsity = 100*nnz/(N*M)
    end

    if returnValue == "print only"
        println("Sparsity = $(sparsity) %")
        return
    elseif returnValue == "return value and print"
        println("Sparsity = $(sparsity) %")
        return sparsity
    elseif returnValue == "return value only"
        return sparsity
    elseif returnValue == "nothing"
        return
    else
        error("Unkown type of output request.")
    end
    
    return sparsity
end

function getValueFromSparMat(A::SparseMatrix,
    row::Int64,
    col::Int64;
    verbose::Bool=false)

    FIR, FIC, nnzVec = A.NVec.FIR, A.MVec.FIC, A.nnzVec
    N, M = length(FIR), length(FIC)

    val = 0

    if row > N || col > M
        error("Requested element is outside the bounds of the matrix.")
    end

    elemFound = false
    keepSearching = true

    if N ≤ M

        myprintln(verbose, "Iterating over row $row of the sparse matrix as "*
        "number of rows $N ≤ number of columns $M.")
        rowPtr = FIR[row]
        while rowPtr != -1 && !elemFound && keepSearching 
            nnzElem = nnzVec[rowPtr, :]
            colNum = nnzElem.NCOL
            if colNum == col
                elemFound = true
                keepSearching = false
                val = nnzElem.Val
                return val
            elseif colNum > col
                keepSearching = false
                myprintln(verbose, "Element not found in this row, "* 
                "Are you sure it is a non-sparse element?")
                return val
            elseif colNum < col
                myprintln(verbose, "Not there yet. Go to next element in row.")
                rowPtr = nnzElem.NIR
            else
                error("How did this even happen? Wasn't supposed to happen.")
            end
        end

        return val

    else

        myprintln(verbose, "Iterating over column $col of the sparse matrix as "*
        "number of rows $N ≥ number of columns $M.")        
        colPtr = FIC[col]
        while colPtr != -1 && !elemFound && keepSearching
            nnzElem = nnzVec[:, colPtr]
            rowNum = nnzElem.NROW
            if rowNum == row
                elemFound = true
                keepSearching = false
                val = nnzElem.Val
                return val
            elseif rowNum > row
                keepSearching = false
                myprintln(verbose, "Element not found in this column, " *
                "Are you sure it is a non-sparse element?")
                return val
            elseif rowNum < row
                myprintln(verbose, "Not there yet. Go to next element in column.")
                colPtr = nnzElem.NIC
            else
                error("How did this even happen? Wasn't supposed to happen.")
            end
        end

        return val
    end
    
    myprintln(true, "This line should never have been reached. Review function.")
end

function dotProductSparLU(A::SparseMatrix,
    row::Int64,
    col::Int64;
    returnType::String="productAndAlpha",
    verbose::Bool=false)

    j = min(row, col)

    prod, α = 0, 0

    productComputed = false

    if j==1
        myprintln(verbose, "First row/column element, no multiplication required.")
        # prod = getValueFromSparMat(A, row, col, verbose=false)
        prod = 0
        α = 0
        productComputed = true    
    end

    if ~productComputed
        myprintln(verbose, "Since the element is NOT a first row/col element,"*
        " multiplication is required.")
        FIR, FIC, nnzVec = A.NVec.FIR, A.MVec.FIC, A.nnzVec

        rowPtr = FIR[row]
        colPtr = FIC[col]

        while rowPtr != -1 && colPtr != -1 && ~productComputed
            nnzElemRow = nnzVec[rowPtr, :]
            irow = nnzElemRow.NCOL    
            nnzElemCol = nnzVec[colPtr, :]
            icol = nnzElemCol.NROW

            if irow ≥ j || icol ≥ j
                myprintln(verbose, "Stopping product computation as at least "*
                "one of the pointers has gone above required bounds.")
                productComputed = true
            elseif irow == icol
                myprintln(verbose, "A good match for multiplication found.")
                myprintln(verbose, "The indices are ($(row), $(irow)) and ($(icol), $(col))")
                myprintln(verbose, " Multiplying them.")
                # prod += getValueFromSparMat(A, row, irow, verbose=false) * getValueFromSparMat(A, icol, col, verbose=false)
                prod += nnzElemRow.Val*nnzElemCol.Val
                α += 1
                rowPtr = nnzElemRow.NIR
                colPtr = nnzElemCol.NIC
                
            elseif irow < icol
                myprintln(verbose, "Row pointer is lagging behind column pointer."*
                " Will increment row pointer by one.")
                rowPtr = nnzElemRow.NIR
            elseif icol < irow
                myprintln(verbose, "Column pointer is lagging behind row pointer."*
                " Will increment column pointer by one.")
                colPtr = nnzElemCol.NIC
            else
                error("Unprepared for this.")
            end
        end
    end


    if returnType == "productOnly"
        return prod
    elseif returnType == "productAndAlpha"
        return prod, α
    elseif returnType == "alphaOnly"
        return α
    else
        error("Unknown return type.")
    end

end

function solveUsingSparseLU(A::SparseMatrix,
    b::Vector;
    method::String="Crout",
    returnType::String="xAndAlphaBeta",
    verbose::Bool = false)
    
    qlu = sparLU(A, verbose=verbose)
    Q = qlu.Q
    α = qlu.α
    y, β_1 = sparForwardSolve(Q, b, verbose=verbose)
    x, β_2 = sparBackwardSolve(Q, y, verbose=verbose)
    β = β_1 + β_2
    numOperations = α + β
    
    if returnType == "xAndAlphaBeta"
        return (x=x, numOperations=numOperations, α=α, β=β)
    elseif returnType == "xOnly"
        return (x=x)
    elseif returnType == "alphaBetaOnly"
        return (numOperations=numOperations, α=α, β=β)
    else
        error("Unknown return type")
    end
end

function sparLU(A::SparseMatrix;
    method::String="Crout",
    verbose::Bool=false)

    FIR, FIC, nnzVec = A.NVec.FIR, A.MVec.FIC, A.nnzVec
    N, M= length(FIR), length(FIC)
    α, fills = 0, 0

    Tinput = eltype(nnzVec.Val)
    if Tinput == Int64
        T = Float64
    else
        T = Tinput
    end

    Q, L, U = [sparseMatrixConstructor(N, M, T=T) for _ in 1:3]

    for j = 1:M
        myprintln(verbose, "Completing column $(j) of the Q Matrix.")
        for k = j:N
            myprintln(verbose, "About to compute Q matrix elements for row $k and column $j .")
            prod, αₖⱼ = dotProductSparLU(Q, k, j, verbose=verbose)
            aₖⱼ = getValueFromSparMat(A, k, j)

            qₖⱼ = aₖⱼ - prod
            if qₖⱼ != 0
                if aₖⱼ == 0
                    myprintln(verbose, "Ahah! Fill-in detected at (row, col) = ($(k), $(k)).")
                    fills += 1
                end
                compElemQ = DataFrame(i = k, j = j, Val = qₖⱼ)
                nnzElemQ = nnzRowConstructor(compElemQ[1, :])
                if method == "Doolittle" && k == j
                    compElemL = DataFrame(i = k, j =j, Val = 1)
                    nnzElemL = nnzRowConstructor(compElemL[1, :])
                    nnzElemU = nnzElemQ
                    updateSparse(U, nnzElemU)
                    myprintln(verbose, "U[$(nnzElemU.NROW), $(nnzElemU.NCOL)] = $(nnzElemU.Val).")
                elseif method == "Crout" && k == j
                    nnzElemL = nnzElemQ
                    # compElemU = DataFrame(i = j, j = k, Val = 1)
                    compElemU = DataFrame(i = k, j = j, Val = 1)
                    nnzElemU = nnzRowConstructor(compElemU[1, :])
                    updateSparse(U, nnzElemU)
                    myprintln(verbose, "U[$(nnzElemU.NROW), $(nnzElemU.NCOL)] = $(nnzElemU.Val).")
                else
                    nnzElemL = nnzElemQ
                end
                updateSparse(Q, nnzElemQ)
                myprintln(verbose, "Q[$(nnzElemQ.NROW), $(nnzElemQ.NCOL)] = $(nnzElemQ.Val).")
                updateSparse(L, nnzElemL)
                myprintln(verbose, "L[$(nnzElemL.NROW), $(nnzElemL.NCOL)] = $(nnzElemL.Val).")
                α += αₖⱼ 
            elseif qₖⱼ == 0
                myprintln(verbose, "Everything is zero. Therefore nothing should be added anywhere.")
            else
                error("Shouldn't be reached.")
            end
        end
        
        if j != M
            myprintln(verbose, "Completing row $(j) of the Q Matrix.")
            qⱼⱼ = getValueFromSparMat(Q, j, j)

            for k = j+1:M
                myprintln(verbose, "About to compute Q matrix elements for row $j and column $k.")
                prod, αⱼₖ = dotProductSparLU(Q, j, k, verbose=verbose)
                aⱼₖ = getValueFromSparMat(A, j, k)
                qⱼₖ = (aⱼₖ - prod)/qⱼⱼ
                if qⱼₖ != 0
                    if aⱼₖ == 0
                        myprintln(verbose, "Ahah! Fill-in detected at (row, col) = ($(k), $(k)).")
                        fills += 1
                    end
                    compElemQ = DataFrame(i = j, j = k, Val = qⱼₖ)
                    nnzElemQ = nnzRowConstructor(compElemQ[1, :])
                    updateSparse(Q, nnzElemQ)
                    myprintln(verbose, "Q[$(nnzElemQ.NROW), $(nnzElemQ.NCOL)] = $(nnzElemQ.Val).")
                    updateSparse(U, nnzElemQ)
                    myprintln(verbose, "U[$(nnzElemQ.NROW), $(nnzElemQ.NCOL)] = $(nnzElemQ.Val).")
                    α += (αⱼₖ + 1)
                elseif qⱼₖ == 0
                    myprintln(verbose, "Everything is zero. Therefore nothing should be added anywhere.")
                else
                    error("Shouldn't be reached.")
                end
            end

        end
    end

    qluMatrices = (Q = Q, L = L, U = U, α = α, fills = fills)
    return qluMatrices
end

function fill_ins(A::Matrix, Q::Matrix)::Int64
    count = 0
    N = size(A, 1)
    M = size(A, 2)
    if size(Q) != (N, M)
        error("Incompatible sizes of matrices. Cannot be compared.")
    end

    for i in N
        for j in M
            if A[i, j] == 0 && Q[i, j] != 0
                count += 1
            end
        end
    end
    return count
end

function sparForwardSolve(Q::SparseMatrix, 
    b::Vector{<:Union{Int64, Float64, ComplexF64}};
    returnType::String="yAndBeta",
    verbose::Bool=false)

    NVec, MVec, nnzVec = Q.NVec, Q.MVec, Q.nnzVec
    FIR = NVec.FIR
    FIC = MVec.FIC
    Vals = nnzVec.Val
    TA = eltype(Vals)
    Tb = eltype(b)
    if TA != Tb
        println("Hmmm.. DataTypes of SparseMatrix and the Vector are NOT the same.")
    end

    N = length(b)
    if N != length(FIR)
        println("Hmmm.. Number of rows in SparseMatrix and Vector are NOT the same.")
    end
    M = length(FIC)
    y = zeros(Tb, N)
    β = 0

    for i = 1:N
        qᵢᵢ = getValueFromSparMat(Q, i, i)
        myprintln(verbose, "Q[$(i), $(i)] = $(qᵢᵢ)")
        dotProductValues = dotProductSparFwd(Q, y, i, verbose=verbose)
        product, β₁ = dotProductValues.product, dotProductValues.β
        y[i] = (b[i] - product)/qᵢᵢ
        β += (β₁ + 1)
        myprintln(verbose, "y[$(i)] = $(y[i])")
        myprintln(verbose, "Multiplication/Division counter is now $(β).")
    end

    if returnType == "yAndBeta"
        return (y=y, β=β)
    elseif returnType == "yOnly"
        return (y=y)
    elseif returnType == "betaOnly"
        return (β=β)
    else
        error("Unknown return type")
    end
end

function dotProductSparFwd(A::SparseMatrix, 
    y::Vector{<:Union{Int64, Float64, ComplexF64}},
    row::Int64;
    returnType::String="productAndBeta",
    verbose::Bool=false)

    product, β = 0, 0

    i = row
    if i == 1
        myprintln(verbose, "First row, don't perform any computation. "*
        "returning zero.")
        if returnType == "productAndBeta"
            return (product=product, β=β)
        elseif returnType == "yOnly"
            return (product=product)
        elseif returnType == "betaOnly"
            return (β=β)
        else
            error("Unknown return type")
        end
    end
        
    NVec, MVec, nnzVec = A.NVec, A.MVec, A.nnzVec
    FIR = NVec.FIR

    rowEnded = false
    diagonalNotTouching = false
    currentElemID = FIR[i]
    if currentElemID == -1 
        rowEnded = true
        myprintln(verbose, "What? Row never starts?")
    else
        nnzElem = nnzVec[currentElemID, :]
        myprintln(verbose, "Starting iterating over row $(i) with elemID = $(currentElemID).")
    end

    while ~rowEnded && ~diagonalNotTouching
        j = nnzElem.NCOL
        myprintln(verbose, "Current element with ID $(currentElemID) is at ($(i), $(j)) "*
        "and has value $(nnzElem.Val)")
        if j == i - 1
            diagonalNotTouching = true
            myprintln(verbose, "Looks like this will be the final " *
            "computation for this row, as diagonal for Q has been reached.")
        end
        qᵢⱼ = nnzElem.Val
        product += qᵢⱼ*y[j]
        myprintln(verbose, "Product is now $(product) after an increment of " * 
        "$(qᵢⱼ)*$(y[j]) = $(qᵢⱼ*y[j]).")
        β += 1
        myprintln(verbose, "Beta counter increased by one to become $(β).")
        currentElemID = nnzElem.NIR

        if currentElemID == -1
            rowEnded = true
            myprintln(verbose, "Looks like the computation will end prematurely, "*
            "as the row is already ended. (This shouldn't really happen in this system.)")
        else
            nnzElem =  nnzVec[currentElemID, :]
        end
    end

    if returnType == "productAndBeta"
        return (product=product, β=β)
    elseif returnType == "productOnly"
        return (product=product)
    elseif returnType == "betaOnly"
        return (β=β)
    else
        error("Unknown return type")
    end
    
    error("Invalid code zone reached.")
end

function sparBackwardSolve(U::SparseMatrix, 
    b::Vector{<:Union{Int64, Float64, ComplexF64}};
    returnType::String="xAndBeta",
    verbose::Bool=false)

    NVec, MVec, nnzVec = U.NVec, U.MVec, U.nnzVec
    FIR = NVec.FIR
    FIC = MVec.FIC
    Vals = nnzVec.Val
    TA = eltype(Vals)
    Tb = eltype(b)
    
    if TA != Tb
        println("Hmmm.. DataTypes of SparseMatrix and the Vector are NOT the same.")
    end

    N = length(b)
    if N != length(FIR)
        println("Hmmm.. Number of rows in SparseMatrix and Vector are NOT the same.")
    end
    
    x = zeros(Tb, N)
    β = 0

    for i = N:-1:1
        myprintln(verbose, "U[$(i), $(i)] is 1 by Crout's method.")
        dotProductValues = dotProductSparBckwd(U, x, i, verbose=verbose)
        product, β₁ = dotProductValues.product, dotProductValues.β
        x[i] = b[i] - product
        β += β₁
        myprintln(verbose, "x[$(i)] = $(x[i])")
        myprintln(verbose, "Multiplication/Division counter is now $(β).")
    end

    if returnType == "xAndBeta"
        return (x=x, β=β)
    elseif returnType == "xOnly"
        return (x=x)
    elseif returnType == "betaOnly"
        return (β=β)
    else
        error("Unknown return type")
    end
end

function dotProductSparBckwd(U::SparseMatrix, 
    x::Vector{<:Union{Int64, Float64, ComplexF64}},
    row::Int64;
    returnType::String="productAndBeta",
    verbose::Bool=false)

    product, β = 0, 0
    
    NVec, MVec, nnzVec = U.NVec, U.MVec, U.nnzVec
    FIR = NVec.FIR
    N = length(FIR)
    
    myprintln(verbose, "Taking the dot product along row $(row) for " *
    "backward substitution.")
    i = row
    if i == N
        myprintln(verbose, "Last row, don't perform any computation. "*
        "returning zero.")
        if returnType == "productAndBeta"
            return (product=product, β=β)
        elseif returnType == "yOnly"
            return (product=product)
        elseif returnType == "betaOnly"
            return (β=β)
        else
            error("Unknown return type")
        end
    end
        
    rowEnded = false
    diagonalBreached = false
    currentElemID = FIR[i]
    if currentElemID == -1 
        rowEnded = true
        myprintln(verbose, "What? Row never starts?")
    else
        nnzElem = nnzVec[currentElemID, :]
        myprintln(verbose, "Starting iterating over row $(i) with elemID = $(currentElemID).")
    end

    if ~rowEnded
        j = nnzElem.NCOL
        myprintln(verbose, "Current element with ID $(currentElemID) is at ($(i), $(j)) "*
        "and has value $(nnzElem.Val)")
        while ~diagonalBreached && ~rowEnded
            if j ≥ i + 1
                diagonalBreached = true
                myprintln(verbose, "Looks like this will be the first " *
                "computation for this row, as diagonal for Q has just been breached.")
            elseif j ≤ i
                myprintln(verbose, "Keep looking to breach the diagonal.")
                currentElemID = nnzElem.NIR
                if currentElemID == -1
                    rowEnded = true
                    myprintln(verbose, "What? Row ended before we even breached the diagonal?")
                else
                    nnzElem = nnzVec[currentElemID, :]
                    j = nnzElem.NCOL
                    myprintln(verbose, "Current element with ID $(currentElemID) is at ($(i), $(j)) "*
                    "and has value $(nnzElem.Val)")
                end
            else
                error("Never supposed to reach here lol.")
            end
        end
    else
        myprintln(verbose, "Looks like there will be no product.")
    end

    while ~rowEnded
        qᵢⱼ = nnzElem.Val
        product += qᵢⱼ*x[j]
        myprintln(verbose, "Product is now $(product) after an increment of " * 
        "$(qᵢⱼ)*$(x[j]) = $(qᵢⱼ*x[j]).")
        β += 1
        myprintln(verbose, "Beta counter increased by one to become $(β).")
        currentElemID = nnzElem.NIR

        if currentElemID == -1
            rowEnded = true
            myprintln(verbose, "Last element in this row reached.")
        else
            nnzElem =  nnzVec[currentElemID, :]
            j = nnzElem.NCOL
            myprintln(verbose, "Current element with ID $(currentElemID) is at ($(i), $(j)) "*
            "and has value $(nnzElem.Val)")
        end
    end

    if returnType == "productAndBeta"
        return (product=product, β=β)
    elseif returnType == "productOnly"
        return (product=product)
    elseif returnType == "betaOnly"
        return (β=β)
    else
        error("Unknown return type")
    end
    
    error("Invalid code zone reached.")

end

function compareSparseAndDense(ASpar::SparseMatrix, ADense::Matrix)
    ASpar2Full = spar2Full(ASpar)
    diff = ADense - ASpar2Full
    non_zero_elements = [(diff[i, j], i, j) for i in 1:size(diff, 1), j in 1:size(diff, 2) if diff[i, j] != 0]
    return non_zero_elements
end

function solveForPowerFlow_Sparse(CDF_DF_List_pu::Vector{DataFrame};
    tolerance::Float64=1e-5,
    itrMax::Int64=30,
    verbose::Bool=false,
    roundDigits::Int64=0)

    powSysData = initializeVectors_pu(CDF_DF_List_pu);
    PSpecified = powSysData.PSpecified;
    QSpecified = powSysData.QSpecified;
    V = powSysData.V;
    δ = powSysData.delta;
    lSlack = powSysData.listOfSlackBuses;
    lPV = powSysData.listOfPVBuses;
    lPQ = powSysData.listOfPQBuses;
    lNonSlack = powSysData.listOfNonSlackBuses;
    nPV = powSysData.nPV;
    nPQ = powSysData.nPQ;
    nSlack = powSysData.nSlack;
    nNonSlack = nPV+nPQ;
    N = nSlack + nNonSlack;
    P = deepcopy(PSpecified)
    Q = deepcopy(QSpecified)
    sparYBus = constructSparseYBus(CDF_DF_List_pu);

    residual = 100
    itr = 1
    while itr ≤ itrMax && residual > tolerance
        myprintln(verbose, "Iteration $(itr) of Power flow.")

        ΔP, ΔQ = computeMismatchesViaSparseYBus(PSpecified, QSpecified, V, δ, sparYBus);
        mismatch = vcat(ΔP[lNonSlack], ΔQ[lPQ])
        myprintln(verbose, "Iteration $(itr): Mismatch = $([round(x, digits=roundDigits) for x in mismatch])")
        P = PSpecified - ΔP;
        Q = QSpecified - ΔQ;
        myprintln(verbose, "Iteration $(itr): P = $([round(x, digits=roundDigits) for x in P])")
        myprintln(verbose, "Iteration $(itr): Q = $([round(x, digits=roundDigits) for x in Q])")

        J = constructSparseJacobian(CDF_DF_List_pu, P, Q, V, δ, sparYBus);

        myprintln(verbose, "Iteration $(itr): Jacobian = $([round(x, digits=roundDigits) for x in spar2Full(J)])")
        correction = solveUsingSparseLU(J, mismatch, verbose=verbose).x
        myprintln(verbose, "Iteration $(itr): Correction = $([round(x, digits=roundDigits) for x in correction])")
        residual = mean(abs.(correction))
        myprintln(verbose, "Iteration $(itr): Residual = $([round(x, digits=6) for x in residual])")
        Δδ = correction[1:nNonSlack];
        ΔVbyV = correction[N:N+nPQ-1];
        δ[lNonSlack] = δ[lNonSlack] + Δδ;
        V[lPQ] .*= (1 .+ (ΔVbyV));
        DataFrame(i=1:N, P=P, Q=Q, V=V, δ=δ)
        itr += 1
    end
    
    if residual ≥ tolerance
        warning("Convergence NOT achieved even after $(itrMax) iterations.\n"*
        "Abandoning run and returning current state variables.")
    elseif residual < tolerance
        myprintln(true, "Convergence achieved after $(itr) iterations.\n") # always true 
    else
        error("No man's land.")
    end

    results = DataFrame(i=1:N, P=P, Q=Q, V=V, δ=δ)
    return results

end;
