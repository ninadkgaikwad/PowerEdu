using DataFrames
using CSV

"""
    Ybus_Builder(CDF_DF_List_pu; kwargs...)

Generates the Y-bus matrix and related matrices for a power system.

# Arguments
- `CDF_DF_List_pu`: A Vector of DataFrames containing the bus data for the power system.

# Optional keyword arguments
- `disableTaps`: A logical value indicating whether to disable tap ratios (default: `false`).
- `sortBy`: A string specifying the order of the Y-bus matrix ('busNumbers' or 'busTypes', default: 'busNumbers').
- `verbose`: A logical value indicating whether to display verbose output (default: `false`).
- `saveTables`: A logical value indicating whether to save the Y-bus and B-matrix as CSV files (default: `false`).
- `saveLocation`: A string specifying the folder location to save the tables (default: 'processedData/').
- `systemName`: A string specifying the name of the power system (default: 'systemNameNOTSpecified').

# Returns
- `ybus`: The Y-bus matrix representing the power system.
- `BMatrix`: The B-matrix representing the power system.
- `b`: A matrix representing the susceptance values of the branches.
- `A`: A matrix representing the connection between branches and buses.
- `branchNames`: A vector containing the names of the branches.
- `E`: A vector of vectors representing the adjacency list of buses.

"""

# Function code
# function Ybus_Builder(BusDataCard_DF::DataFrame, BranchDataCard_DF::DataFrame;
function Ybus_Builder(CDF_DF_List_pu::Vector{DataFrame};
    disableTaps::Bool = false,
    sortBy::String = "busNumbers",
    verbose::Bool = false,
    saveTables::Bool = false,
    saveLocation::String = "processedData/",
    systemName::String = "systemNameNOTSpecified")
#Later, systemName should be read from CDF_DF_List_pu itself, so there
#shouldn't be a need for specifying systemName

    BusDataCard_DF = CDF_DF_List_pu[2];
    BranchDataCard_DF = CDF_DF_List_pu[3];

    N = size(BusDataCard_DF, 1)
    numBranch = size(BranchDataCard_DF, 1)

    ybus = zeros(ComplexF64, N, N)
    BMatrix = zeros(ComplexF64, N, N)
    E = Array{Vector{Int64}}(undef, N)
    b = zeros(Float64, numBranch, numBranch)

    for bus in 1:N
        E[bus] = Vector{Int64}()
    end

    A = zeros(Float64, numBranch, N)
    branchNames = Vector{String}(undef, numBranch)

    for branch = 1:numBranch
        currentBranch = BranchDataCard_DF[branch, :]
        i = currentBranch.Tap_Bus_Num
        k = currentBranch.Z_Bus_Num
        branchNames[branch] = "$(i) to $(k)"
        A[branch, i] = 1
        A[branch, k] = -1
        b[branch, branch] = currentBranch.B_pu

        if disableTaps
            a = 1
        elseif currentBranch.a != 0
            a = currentBranch.a
        else
            a = 1
        end

        y_ik = 1/(currentBranch.R + im*currentBranch.X)
        ybus[i, i] += y_ik/(a^2) + currentBranch.B_pu / 2
        ybus[k, k] += y_ik + currentBranch.B_pu / 2
        ybus[i, k] = -y_ik/a
        ybus[k, i] = -y_ik/a

        push!(E[i], k)
        push!(E[k], i)
    end

    for bus = 1:N
        ybus[bus, bus] += BusDataCard_DF.G_pu[bus] + im*BusDataCard_DF.B_pu[bus]
    end

    BMatrix = -imag(ybus)
    # Sort Y-bus matrix
    if sortBy == "busNumbers"
        
        rowNamesSingleDigits = ["Gen0$i" for i in 1:1:min(N, 9)]
        rowNamesDoubleDigits = ["Gen$i" for i in 10:1:N]
        rowNames = [vcat(rowNamesSingleDigits, rowNamesDoubleDigits)]
        # @show rowNames = [string(i) for i in 1:N]
        #might wanna change the names to be Gen01, Gen02, ... , Gen14.
        tag = ""
    elseif sortBy == "busTypes"
        ybusByTypes, rowNamesByTypes = sortMatrixByBusTypes(BusDataCard_DF, ybus)
        ybus = ybusByTypes
        rowNames  = rowNamesByTypes
        BMatrixByTypes, rowNamesByTypes = sortMatrixByBusTypes(BusDataCard_DF, BMatrix)
        BMatrix = BMatrixByTypes
        tag = "_sortedByBusTypes"
    end

    @show ybusTable = DataFrame(ybus, Symbol.(rowNames))
    BMatrixTable = DataFrame(BMatrix, Symbol.(rowNames))

    if verbose
        println("Y-bus Matrix:")
        show(stdout, "text/plain", ybus)
        println("\nB-Matrix:")
        show(stdout, "text/plain", BMatrix)
        println("\nBranch Names:")
        show(stdout, "text/plain", branchNames)
        println("\nA-Matrix:")
        show(stdout, "text/plain", A)
        println("\nb-Matrix:")
        show(stdout, "text/plain", b)
        println("\nE (Adjacency list):")
        show(stdout, "text/plain", E)
    end

    if saveTables
        fileType = ".csv"
        filenameYBus = "$saveLocation$systemName/YBus$tag$fileType"
        filenameBMatrix = "$saveLocation$systemName/BMatrix$tag$fileType"
        CSV.write(filenameYBus, ybusTable)
        CSV.write(filenameBMatrix, BMatrixTable)
    end

    return [ybus, BMatrix, b, A, branchNames, E]

end