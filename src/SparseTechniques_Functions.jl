# SparseTechniquesInPowerSystems.jl

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
compMatrix = DataFrame(val = values, i = rows, j = cols)

# Convert the compressed matrix to a full matrix
matFull = compressed2Full(compMatrix)
"""
function compressed2Full(compMatrix::DataFrame)
	# Use SparseArrays's sparse function to conveniently convert the compressed 
	# matrix (i, j, val) into Compressed Storage Column CSC format
	# I don't care how it does it. 
	# This function is only to be called for testing purposes anyway.
    sparseMatrix = sparse(compMatrix.i, compMatrix.j, compMatrix.val)
	# Convert the sparse matrix into the full matrix.
    fullMatrix = Matrix(sparseMatrix)
end