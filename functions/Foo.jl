# in file Foo.jl

module Foo # Declare the module

export initializeVectors
export sortMatrixByBusTypes
export ybusGenerator # write what will be accessible from outside

# Write your functions...

# function initializeVectors(busData; MVAb=100)
function initializeVectors(busData; MVAb=100)


    N = size(busData, 1)
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
        bus = busData.bus[i]
        delta[bus] = 0.0000
        if busData.busType[i] == 0
            nPQ += 1
            listOfPQBuses[nPQ] = bus
            n += 1
            listOfNonSlackBuses[n] = bus
            V[bus] = 1.0000
        elseif busData.busType[i] == 2
            nPV += 1
            listOfPVBuses[nPV] = bus
            n += 1
            listOfNonSlackBuses[n] = bus
            V[bus] = busData.Vset[i]
        elseif busData.busType[i] == 3
            nSlack += 1
            listOfSlackBuses[nSlack] = bus
            V[bus] = busData.Vset[i]
        end
        PSpecified[bus] = busData.PG[i] / MVAb - busData.PL[i] / MVAb
        QSpecified[bus] = busData.QG[i] / MVAb - busData.QL[i] / MVAb
    end

    listOfSlackBuses = reshape(listOfSlackBuses[1:nSlack], nSlack)
    listOfSlackBuses = reshape(listOfSlackBuses[1:nSlack], nSlack)
    listOfPVBuses = reshape(listOfPVBuses[1:nPV], nPV)
    listOfPQBuses = reshape(listOfPQBuses[1:nPQ], nPQ)
    # @show typeof(listOfPQBuses)

    return [PSpecified, QSpecified, V, delta, listOfSlackBuses, listOfPVBuses, listOfPQBuses, listOfNonSlackBuses, nSlack, nPV, nPQ]
end


function sortMatrixByBusTypes(busData, ybus)

    @show outputs  = initializeVectors(busData)
    @show listOfSlackBuses, listOfPVBuses, listOfPQBuses = outputs[5], outputs[6], outputs[7]
    # @show typeof(listOfPVBuses)
    @show newOrder = vcat(listOfSlackBuses, listOfPVBuses, listOfPQBuses)
    # @show newOrder = [listOfSlackBuses listOfPVBuses listOfPQBuses]
    ybusByTypes = ybus[newOrder, newOrder]

    @show typeof(newOrder)
    @show rowNamesByTypes = [string(i) for i in newOrder]

    return ybusByTypes, rowNamesByTypes
end

using DataFrames
using CSV

"""
    ybusGenerator(busData, branchData; kwargs...)

Generates the Y-bus matrix and related matrices for a power system.

# Arguments
- `busData`: A DataFrame containing the bus data for the power system.
- `branchData`: A DataFrame containing the branch data for the power system.

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

busData = DataFrame(G=[0.1, 0.05, 0.2], B=[0.2, 0.15, 0.25])
branchData = DataFrame(i=[1, 2, 3], j=[2, 3, 1], R=[0.1, 0.2, 0.15], X=[0.2, 0.3, 0.25], B=[0.05, 0.1, 0.08])
ybus, BMatrix, b, A, branchNames, E = ybusGenerator(busData, branchData, verbose=true, saveTables=true)
"""

# Function code
function ybusGenerator(busData::DataFrame, branchData::DataFrame;
    disableTaps::Bool = false,
    sortBy::String = "busNumbers",
    verbose::Bool = false,
    saveTables::Bool = false,
    saveLocation::String = "processedData/",
    systemName::String = "systemNameNOTSpecified")

    N = size(busData, 1)
    numBranch = size(branchData, 1)

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
        currentBranch = branchData[branch, :]
        i = currentBranch.i
        k = currentBranch.j
        branchNames[branch] = "$(i) to $(k)"
        A[branch, i] = 1
        A[branch, k] = -1
        b[branch, branch] = currentBranch.B

        if disableTaps
            a = 1
        elseif currentBranch.a != 0
            a = currentBranch.a
        else
            a = 1
        end

        y_ik = 1/(currentBranch.R + im*currentBranch.X)
        ybus[i, i] += y_ik/(a^2) + currentBranch.B / 2
        ybus[k, k] += y_ik + currentBranch.B / 2
        ybus[i, k] = -y_ik/a
        ybus[k, i] = -y_ik/a

        push!(E[i], k)
        push!(E[k], i)
    end

    for bus = 1:N
        ybus[bus, bus] += busData.G[bus] + im*busData.B[bus]
    end

    BMatrix = -imag(ybus)
    # Sort Y-bus matrix
    if sortBy == "busNumbers"
        @show rowNames = [string(i) for i in 1:N]
        #might wanna change the names to be Gen01, Gen02, ... , Gen14.
        tag = ""
    elseif sortBy == "busTypes"
        ybusByTypes, rowNamesByTypes = sortMatrixByBusTypes(busData, ybus)
        ybus = ybusByTypes
        rowNames  = rowNamesByTypes
        BMatrixByTypes, rowNamesByTypes = sortMatrixByBusTypes(busData, BMatrix)
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

# busData = DataFrame(G=[0.1, 0.05, 0.2], B=[0.2, 0.15, 0.25])
# branchData = DataFrame(i=[1, 2, 3], j=[2, 3, 1], R=[0.1, 0.2, 0.15], X=[0.2, 0.3, 0.25], B=[0.05, 0.1, 0.08])
# ybus, BMatrix, b, A, branchNames, E = ybusGenerator(busData, branchData, verbose=true, saveTables=true)


end # Foo