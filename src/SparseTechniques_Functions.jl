# SparseTechniquesInPowerSystems.jl

using SparseArrays
using DataFrames

include("Helper_Functions.jl")

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
sparmat(compMatrix::DataFrame; verbose::Bool = false)

The `sparmat` function creates two DataFrames, `NVec` and `nnzVec`, which together represent a sparse matrix. It iterates over the `compMatrix` DataFrame, constructs `nnzVec` rows using the `nnzRowConstructor` function, and updates `NVec` and `nnzVec` accordingly.

Arguments:
- `compMatrix::DataFrame`: The DataFrame containing the compressed matrix information.
- `verbose::Bool`: (optional) A Boolean value indicating whether to display verbose output during the construction process. Default is `false`.

Returns:
- A tuple `(NVec, nnzVec)` representing the two DataFrames `NVec` and `nnzVec` that represent the sparse matrix.

Example:
```julia
compMatrix = DataFrame(ID = [1, 2, 3], Val = [0.5, 0.3, 0.8], i = [2, 1, 3], j = [3, 2, 1])
NVec, nnzVec = sparmat(compMatrix, verbose = true)
"""
function sparmat(compMatrix::DataFrame;
	verbose::Bool = false)

	# N = maximum([compMatrix.i compMatrix.j])
    N = maximum(compMatrix.i)
    M = maximum(compMatrix.j)
	# (firs, fics) = (repeat([-1], N), repeat([-1], N))
    (firs, fics) = (repeat([-1], N), repeat([-1], M))


	# NVec = DataFrame(FIR = firs, FIC = fics)
    NVec = DataFrame(FIR = firs)
    MVec = DataFrame(FIC = fics)
	nnzVec = DataFrame(ID = Int64[], Val = ComplexF64[], NROW = Int64[], NCOL = Int64[], NIR = Int64[], NIC = Int64[])

	numElems = size(compMatrix, 1)
	
	for elemNum in 1:numElems
		compElem = compMatrix[elemNum, :]

		nnzElem = nnzRowConstructor(compElem)

		NVec, nnzVec = updateSparse(NVec, MVec, nnzVec, nnzElem, type="replace", verbose=verbose)
	end

	return (NVec=NVec, MVec=MVec, nnzVec=nnzVec)
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
    checkElementIntoRow(NVec::DataFrame, nnzVec::DataFrame, nnzElem::DataFrameRow;
                        type::String = "replace", verbose::Bool = false)

Check and insert an element into a row of a sparse matrix.

This function checks if the given `nnzElem` can be inserted into the row specified by `nnzElem.NROW`
in the `nnzVec` DataFrame. If the row is empty, the element is inserted as the first element. If there
is already an element present, the function determines the correct position for the element based on
the `nnzElem.NCOL` value. The position can be before the incumbent element, after it, or it may result
in a tie, which can be resolved using the `resolveTie` function.

## Arguments
- `NVec::DataFrame`: DataFrame containing additional information about the sparse matrix.
- `nnzVec::DataFrame`: DataFrame representing the non-zero elements of the sparse matrix.
- `nnzElem::DataFrameRow`: DataFrameRow representing the new element to be inserted.
- `type::String`: Optional. The type of tie resolution. Default is "replace". Possible values are "replace" and "add".
- `verbose::Bool`: Optional. If set to `true`, print verbose output. Default is `false`.

## Returns
- `NVec::DataFrame`: Updated DataFrame `NVec` with the modified first-in-row information.
- `nnzVec::DataFrame`: Updated DataFrame `nnzVec` with the newly inserted element or resolved tie.
- `nnzElem::DataFrameRow`: Updated DataFrameRow `nnzElem` with the updated element information.
- `updateFlag::Bool`: A flag indicating if the sparse matrix was updated.

## Example
```julia
nnzVec = DataFrame(ID = 1:5, NROW = [1, 1, 2, 2, 3], NCOL = [1, 2, 1, 2, 1], Val = [1.0, 2.0, 3.0, 4.0, 5.0], NIR = [2, -1, 4, -1, -1], NIC = [-1, -1, -1, -1, -1])
nnzElem = DataFrameRow([0, 3, 2, 3, 1], ["ID", "NROW", "NCOL", "Val", "NIR", "NIC"])
NVec, nnzVec, nnzElem, updateFlag = checkElementIntoRow(NVec, nnzVec, nnzElem, type="replace", verbose=true)
"""
function checkElementIntoRow(sparMat::NamedTuple{(:NVec, :MVec, :nnzVec), Tuple{DataFrame, DataFrame, DataFrame}},
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
    checkElementIntoColumn(NVec::DataFrame, nnzVec::DataFrame, nnzElem::DataFrameRow;
    type::String = "replace", verbose::Bool = false)

Check and insert an element into a column of a sparse matrix.

This function is analogous to the `checkElementIntoRow` function and serves the same purpose,
but it operates on columns instead of rows. It checks if the given `nnzElem` can be inserted
into the column specified by `nnzElem.NCOL` in the `nnzVec` DataFrame. If the column is empty,
the element is inserted as the first element. If there is already an element present, the function
determines the correct position for the element based on the `nnzElem.NROW` value. The position
can be before the incumbent element, after it, or it may result in a tie, which can be resolved
using the `resolveTie` function.

## Arguments
- `NVec::DataFrame`: DataFrame containing additional information about the sparse matrix.
- `nnzVec::DataFrame`: DataFrame representing the non-zero elements of the sparse matrix.
- `nnzElem::DataFrameRow`: DataFrameRow representing the new element to be inserted.
- `type::String`: Optional. The type of tie resolution. Default is "replace". Possible values are "replace" and "add".
- `verbose::Bool`: Optional. If set to `true`, print verbose output. Default is `false`.

## Returns
- `NVec::DataFrame`: Updated DataFrame `NVec` with the modified first-in-column information.
- `nnzVec::DataFrame`: Updated DataFrame `nnzVec` with the newly inserted element or resolved tie.
- `nnzElem::DataFrameRow`: Updated DataFrameRow `nnzElem` with the updated element information.
- `updateFlag::Bool`: A flag indicating if the sparse matrix was updated.

## Example
```julia
nnzVec = DataFrame(ID = 1:5, NROW = [1, 1, 2, 2, 3], NCOL = [1, 2, 1, 2, 1], Val = [1.0, 2.0, 3.0, 4.0, 5.0], NIR = [2, -1, 4, -1, -1], NIC = [-1, -1, -1, -1, -1])
nnzElem = DataFrameRow([0, 3, 2, 3, 1], ["ID", "NROW", "NCOL", "Val", "NIR", "NIC"])
NVec, nnzVec, nnzElem, updateFlag = checkElementIntoColumn(NVec, nnzVec, nnzElem, type="replace", verbose=true)
"""
function checkElementIntoColumn(sparMat::NamedTuple{(:NVec, :MVec, :nnzVec), Tuple{DataFrame, DataFrame, DataFrame}},
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
    updateSparse(NVec::DataFrame, nnzVec::DataFrame, nnzElem::DataFrameRow;
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
- `NVec::DataFrame`: DataFrame containing additional information about the sparse matrix.
- `nnzVec::DataFrame`: DataFrame representing the non-zero elements of the sparse matrix.
- `nnzElem::DataFrameRow`: DataFrameRow representing the new element to be inserted.
- `type::String`: Optional. The type of tie resolution. Default is "replace". Possible values are "replace" and "add".
- `verbose::Bool`: Optional. If set to `true`, print verbose output. Default is `false`.

## Returns
- `NVec::DataFrame`: Updated DataFrame `NVec` with the modified first-in-row and first-in-column information.
- `nnzVec::DataFrame`: Updated DataFrame `nnzVec` with the newly inserted element or resolved tie.

## Example
```julia
NVec = DataFrame(FIR = [-1, -1, -1], FIC = [-1, -1, -1])
nnzVec = DataFrame(ID = 1:5, NROW = [1, 1, 2, 2, 3], NCOL = [1, 2, 1, 2, 1], Val = [1.0, 2.0, 3.0, 4.0, 5.0], NIR = [2, -1, 4, -1, -1], NIC = [-1, -1, -1, -1, -1])
nnzElem = DataFrameRow([0, 3, 2, 3, 1], ["ID", "NROW", "NCOL", "Val", "NIR", "NIC"])
NVec, nnzVec = updateSparse(NVec, nnzVec, nnzElem, type="replace", verbose=true)
"""
function updateSparse(sparMat::NamedTuple{(:NVec, :MVec, :nnzVec), Tuple{DataFrame, DataFrame, DataFrame}},
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

"""
    constructSparseYBus(CDF_DF_List_pu::Vector{DataFrame};
    disableTaps::Bool = false,
    sortBy::String = "busNumbers",
    verbose::Bool = false,
    saveTables::Bool = false,
    saveLocation::String = "processedData/")

This function constructs a sparse YBus representation from the provided 
    `busData` and `branchData` extracted from `CDF_DF_List_pu`.

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

    (firs, fics) = (repeat([-1], N), repeat([-1], N))

	NVec = DataFrame(FIR = firs)
    MVec = DataFrame(FIC = fics)
	nnzVec = DataFrame(ID = Int64[], Val = ComplexF64[], NROW = Int64[], NCOL = Int64[], NIR = Int64[], NIC = Int64[])
    sparMat = (NVec=NVec, MVec=MVec, nnzVec=nnzVec)

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
            
            compElem = compElem[1, :] # I don't know how to directly initialize a DataFrameRow
            # so I awkwardly extract a row from the DataFrame
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
    computeMismatchesViaSparseYBus(PSpecified, QSpecified, V, delta, NYBus, nnzYBus)

The `computeMismatchesViaSparseYBus` function calculates the mismatches between specified active and reactive power values and the computed power values using the sparse YBus matrix.

# Arguments
- `PSpecified::Vector{Float64}`: A vector of specified active power values for each bus.
- `QSpecified::Vector{Float64}`: A vector of specified reactive power values for each bus.
- `V::Vector{Float64}`: A vector of voltage magnitudes for each bus.
- `delta::Vector{Float64}`: A vector of voltage phase angles for each bus.
- `NYBus::DataFrame`: A DataFrame representing the non-zero elements of the YBus matrix, with columns `FIR`, `NIR`, and `Val`.
- `nnzYBus::DataFrame`: A DataFrame containing the non-zero elements of the YBus matrix, with columns `NCOL` and `Val`.

# Returns
- `deltaP::Vector{Float64}`: A vector representing the mismatches between specified and computed active power values.
- `deltaQ::Vector{Float64}`: A vector representing the mismatches between specified and computed reactive power values.

# Details
The function iterates over each bus in the system and calculates the active and reactive power mismatches using the specified power values, voltage magnitudes, voltage phase angles, and the non-zero elements of the YBus matrix.

The active power mismatch `deltaP` is computed as the difference between the specified active power (`PSpecified`) and the computed active power (`P`) for each bus.

The reactive power mismatch `deltaQ` is computed as the difference between the specified reactive power (`QSpecified`) and the computed reactive power (`Q`) for each bus.

The function assumes that the `NYBus` and `nnzYBus` DataFrames contain the necessary information for building the sparse YBus matrix.

# Example
```julia
PSpecified = [100.0, 200.0, 150.0]
QSpecified = [50.0, 100.0, 75.0]
V = [1.0, 1.0, 1.0]
delta = [0.0, 0.0, 0.0]
NYBus = DataFrame(FIR = [1, 2, 3], NIR = [-1, -1, -1], Val = [1.0, 2.0, 3.0])
nnzYBus = DataFrame(NCOL = [1, 2, 3], Val = [0.5, 0.6, 0.7])

deltaP, deltaQ = computeMismatchesViaSparseYBus(PSpecified, QSpecified, V, delta, NYBus, nnzYBus)
"""
function computeMismatchesViaSparseYBus(PSpecified::Vector{Float64}, 
    QSpecified::Vector{Float64}, 
    V::Vector{Float64}, 
    delta::Vector{Float64}, 
    sparMat::NamedTuple{(:NVec, :MVec, :nnzVec), Tuple{DataFrame, DataFrame, DataFrame}})

    NVec = sparMat.NVec
    MVec = sparMat.MVec
    nnzVec = sparMat.nnzVec

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
    nnzYBus::DataFrame,
    NYBus::DataFrame;
    verbose::Bool=false,
    saveTable::Bool=false,
    processedDataFolder::String="processedData/")

    powSysData = initializeVectors_pu(CDF_DF_List_pu)
    nPV = powSysData.nPV
    nPQ = powSysData.nPQ
    lPV = powSysData.listOfPVBuses
    lPQ = powSysData.listOfPQBuses

    J11 = constructSparseJacobianSubMatrix(CDF_DF_List_pu, type="J11")
    J12 = constructSparseJacobianSubMatrix(CDF_DF_List_pu, type="J12")
    J21 = constructSparseJacobianSubMatrix(CDF_DF_List_pu, type="J21")
    J22 = constructSparseJacobianSubMatrix(CDF_DF_List_pu, type="J22")

    J = combineJacobianSubmatrices(J11, J12, J21, J22)
    return J
end

function combineJacobianSubmatrices(J11::Tuple{DataFrame, DataFrame}, 
    J12::Tuple{DataFrame, DataFrame},
    J21::Tuple{DataFrame, DataFrame},
    J22::Tuple{DataFrame, DataFrame};
    verbose=false)

    J1 = hcatSparse(J11, J12)
    J2 = hcatSparse(J21, J22)
    J = vcatSparse(J1, J2)
    return J
end

function hcatSparse(matLeft::Tuple{DataFrame, DataFrame}, 
    matRight::Tuple{DataFrame, DataFrame};
    verbose::Bool = false)

    return matHorz
end

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

function spar2Full(sparMat::NamedTuple{(:NVec, :MVec, :nnzVec), Tuple{DataFrame, DataFrame, DataFrame}};
    readMethod::String = "row-wise",
    verbose::Bool = false)

    NVec = sparMat.NVec
    MVec = sparMat.MVec
    nnzVec = sparMat.nnzVec

    FIR = NVec.FIR
    FIC = MVec.FIC
    N = length(FIR)
    M = length(FIC)
    nnz = length(nnzVec.ID)
    mat = zeros(ComplexF64, N, M)

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