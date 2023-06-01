# in file Foo.jl
module Foo # Declare the module

export initializeVectors_pu
export sortMatrixByBusTypes
export ybusGenerator # write what will be accessible from outside

# Write your functions...

# function initializeVectors_pu(CDF_DF_List_pu)
function initializeVectors_pu(CDF_DF_List_pu)
    
    busData_pu = CDF_DF_List_pu[2]

    N = size(busData_pu, 1)
    PSpecified = zeros(N)
    QSpecified = zeros(N)
    V = zeros(N)
    delta = zeros(N)
    listOfPQBuses = zeros(Int64, N)
    listOfPVBuses = zeros(Int64, N)
    nPQ = 0
    nPV = 0
    n = 0
    nSlack = 0
    listOfNonSlackBuses = zeros(Int64, N)
    listOfSlackBuses = zeros(Int64, N)

    for i = 1:N
        bus = busData_pu.Bus_Num[i]
        delta[bus] = 0.0000
        if busData_pu.Type[i] == 0
            nPQ += 1
            listOfPQBuses[nPQ] = bus
            n += 1
            listOfNonSlackBuses[n] = bus
            V[bus] = 1.0000
        elseif busData_pu.Type[i] == 2
            nPV += 1
            listOfPVBuses[nPV] = bus
            n += 1
            listOfNonSlackBuses[n] = bus
            V[bus] = busData_pu.Desired_V_pu[i]
        elseif busData_pu.Type[i] == 3
            nSlack += 1
            listOfSlackBuses[nSlack] = bus
            V[bus] = busData_pu.Desired_V_pu[i]
        end
        PSpecified[bus] = busData_pu.Gen_MW[i] - busData_pu.Load_MW[i]
        QSpecified[bus] = busData_pu.Gen_MVAR[i] - busData_pu.Load_MVAR[i]
    end

    listOfSlackBuses = reshape(listOfSlackBuses[1:nSlack], nSlack)
    listOfSlackBuses = reshape(listOfSlackBuses[1:nSlack], nSlack)
    listOfPVBuses = reshape(listOfPVBuses[1:nPV], nPV)
    listOfPQBuses = reshape(listOfPQBuses[1:nPQ], nPQ)

    return [PSpecified, QSpecified, V, delta, listOfSlackBuses, listOfPVBuses, listOfPQBuses, listOfNonSlackBuses, nSlack, nPV, nPQ]
end


function sortMatrixByBusTypes(CDF_DF_List_pu, ybus)
    busData_pu = CDF_DF_List_pu[2]
    outputs  = initializeVectors_pu(CDF_DF_List_pu)
    listOfSlackBuses, listOfPVBuses, listOfPQBuses = outputs[5], outputs[6], outputs[7]
    newOrder = vcat(listOfSlackBuses, listOfPVBuses, listOfPQBuses)
    ybusByTypes = ybus[newOrder, newOrder]

    typeof(newOrder)
    rowNamesByTypes = [string(i) for i in newOrder]

    return ybusByTypes, rowNamesByTypes
end

"""
    extractSystemName(CDF_DF_List::Vector{DataFrame})

Extracts the system name from a vector of DataFrames (assumed to be its CDF_DF_List).

# Arguments
- `CDF_DF_List::Vector{DataFrame}`: A vector of DataFrames.

# Returns
- `bus_name::AbstractString`: The extracted system name, in the format "{prefix}_{bus_number}" if found, or "Bus number not found." otherwise.
"""
function extractSystemName(CDF_DF_List::Vector{DataFrame})
    using 
    header_CDF = CDF_DF_List[1]
    vector_string = header_CDF.Case_ID
    pattern = r"(\D+)\s*(\d+)"

    match_obj = match(pattern, vector_string)
    if match_obj !== nothing
        prefix = strip(match_obj.captures[1])
        bus_number = match_obj.captures[2]
        bus_name = string(prefix, "_", bus_number)
        return bus_name
    else
        return "Bus number not found.\n"
    end
end


"""
    ybusGenerator(busData_pu, branchData_pu; kwargs...)

Generates the Y-bus matrix and related matrices for a power system.

# Arguments
- `busData_pu`: A DataFrame containing the bus data for the power system.
- `branchData_pu`: A DataFrame containing the branch data for the power system.

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

# Example
```julia
using DataFrames, CSV

busData_pu = DataFrame(G=[0.1, 0.05, 0.2], B=[0.2, 0.15, 0.25])
branchData_pu = DataFrame(i=[1, 2, 3], j=[2, 3, 1], R=[0.1, 0.2, 0.15], X=[0.2, 0.3, 0.25], B=[0.05, 0.1, 0.08])
ybus, BMatrix, b, A, branchNames, E = ybusGenerator(busData_pu, branchData_pu, verbose=true, saveTables=true)
"""
function ybusGenerator(CDF_DF_List_pu::Vector{DataFrame};
    disableTaps::Bool = false,
    sortBy::String = "busNumbers",
    verbose::Bool = false,
    saveTables::Bool = false,
    saveLocation::String = "processedData/",
    systemName::String = "systemNameNOTSpecified")

    busData_pu = CDF_DF_List_pu[2]
    branchData_pu = CDF_DF_List_pu[3]
    N = size(busData_pu, 1)
    numBranch = size(branchData_pu, 1)

    ybus = zeros(ComplexF64, N, N)
    BMatrix = zeros(ComplexF64, N, N)
    E = Array{Vector{Int64}}(undef, N)
    b = zeros(Float64, numBranch, numBranch)

    for i in 1:N
        E[i] = Vector{Int64}()
    end
    A = zeros(Float64, numBranch, N)
    branchNames = Vector{String}(undef, numBranch)

    for branch = 1:numBranch
        currentBranch = branchData_pu[branch, :]
        # vscodedisplay(currentBranch)
        i = currentBranch.Tap_Bus_Num
        k = currentBranch.Z_Bus_Num
        branchNames[branch] = "$(i) to $(k)"
        A[branch, i] = 1
        A[branch, k] = -1
        b[branch, branch] = currentBranch.B_pu

        if disableTaps
            a = 1
        elseif currentBranch.Transformer_t != 0
            a = currentBranch.Transformer_t
        else
            a = 1
        end

        y_ik = 1/(currentBranch.R_pu + im*currentBranch.X_pu)
        ybus[i, i] += y_ik/(a^2) + currentBranch.B_pu / 2
        ybus[k, k] += y_ik + currentBranch.B_pu / 2
        ybus[i, k] = -y_ik/a
        ybus[k, i] = -y_ik/a

        push!(E[i], k)
        push!(E[k], i)
    end

    for bus = 1:N
        ybus[bus, bus] += busData_pu.G_pu[bus] + im*busData_pu.B_pu[bus]
    end

    BMatrix = -imag(ybus)
    # Sort Y-bus matrix
    if sortBy == "busNumbers"
        rowNames = [string(i) for i in 1:N]
        #might wanna change the names to be Gen01, Gen02, ... , Gen14.
        tag = ""
    elseif sortBy == "busTypes"
        ybusByTypes, rowNamesByTypes = sortMatrixByBusTypes(CDF_DF_List_pu, ybus)
        ybus = ybusByTypes
        rowNames  = rowNamesByTypes
        BMatrixByTypes, rowNamesByTypes = sortMatrixByBusTypes(CDF_DF_List_pu, BMatrix)
        BMatrix = BMatrixByTypes
        tag = "_sortedByBusTypes"
    end

    ybusTable = DataFrame(ybus, Symbol.(rowNames))
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

# Test case generated by ChatGPT
# busData_pu = DataFrame(G=[0.1, 0.05, 0.2], B=[0.2, 0.15, 0.25])
# branchData_pu = DataFrame(i=[1, 2, 3], j=[2, 3, 1], R=[0.1, 0.2, 0.15], X=[0.2, 0.3, 0.25], B=[0.05, 0.1, 0.08])
# ybus, BMatrix, b, A, branchNames, E = ybusGenerator(busData_pu, branchData_pu, verbose=true, saveTables=true)


end # Foo