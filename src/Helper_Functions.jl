# Helper_Functions.jl

"""
    createFolderIfNotExisting(systemName::AbstractString, folderPath::AbstractString = "processedData/")

Create a folder with the specified system name if it doesn't already exist in the specified folder path.

## Arguments
- `systemName::AbstractString`: The name of the system for which you want to create the folder.
- `folderPath::AbstractString (optional)`: The path where the folder should be created. Default is "processedData/".

## Example
```julia
createFolderIfNotExisting("IEEE_14")
"""
function createFolderIfNotExisting(systemName::AbstractString, folderPath::AbstractString = "processedData/")
    folder_name = folderPath * systemName

    if !isdir(folder_name)
        mkdir(folder_name)
        println("Folder '$folder_name' created successfully.")
    else
        println("Folder '$folder_name' already exists.")
    end
end

"""
myprintln(verbose::Bool, args...)

The `myprintln` function is used to conditionally print a variable number of arguments based on the value of the `verbose` parameter.

Arguments:
- `verbose::Bool`: A Boolean value indicating whether to print the provided arguments or not.
- `args...`: Variable number of arguments to be printed.

Example:
```julia
myprintln(true, "This", "is", "verbose")    # prints: This is verbose
myprintln(false, "This", "is", "quiet")     # does not print anything
"""
function myprintln(verbose::Bool, args...)
    if verbose
        println(args...)
    end
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
    header_CDF = CDF_DF_List[1]
    vector_encapsulating_string = header_CDF.Case_ID #A Vector of strings (size: 1x1)
    fullString = vector_encapsulating_string[1]
    pattern = r"(\D+)\s*(\d+)"

    match_obj = match(pattern, fullString)
    if match_obj !== nothing
        prefix = strip(match_obj.captures[1])
        bus_number = match_obj.captures[2]
        bus_name = string(prefix, "_", bus_number)
        return bus_name
    elsen
        return "Bus number not found.\n"
    end
end

"""
    initializeVectors_pu(CDF_DF_List_pu::Vector{DataFrame};
    busTypes::String = "current") :: NamedTuple{(:PSpecified, :QSpecified, :V, :delta, :listOfSlackBuses, :listOfPVBuses, :listOfPQBuses, :listOfNonSlackBuses, :nSlack, :nPV, :nPQ), Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int64}, Vector{Int64}, Vector{Int64}, Vector{Int64}, Int64, Int64, Int64}}

Initialize vectors for power system analysis in per unit (pu) values.

## Arguments
- `CDF_DF_List_pu::Vector{DataFrame}`: A vector of DataFrames containing bus and branch data in per unit values.
- `busTypes::String`: Optional. The category of bus types to consider. Default is "current". Possible values are "current" and "original".

## Returns
- `PSpecified::Vector{Float64}`: A vector representing the specified active power values for each bus.
- `QSpecified::Vector{Float64}`: A vector representing the specified reactive power values for each bus.
- `V::Vector{Float64}`: A vector representing the voltage magnitudes for each bus.
- `delta::Vector{Float64}`: A vector representing the voltage phase angles for each bus.
- `listOfSlackBuses::Vector{Int64}`: A vector containing the bus numbers of the slack buses.
- `listOfPVBuses::Vector{Int64}`: A vector containing the bus numbers of the PV buses.
- `listOfPQBuses::Vector{Int64}`: A vector containing the bus numbers of the PQ buses.
- `listOfNonSlackBuses::Vector{Int64}`: A vector containing the bus numbers of the non-slack buses.
- `nSlack::Int64`: The number of slack buses.
- `nPV::Int64`: The number of PV buses.
- `nPQ::Int64`: The number of PQ buses.

## Details
The `initializeVectors_pu` function initializes the vectors required for power system analysis in per unit (pu) values. It extracts relevant information from the provided `CDF_DF_List_pu` vector of DataFrames, including bus data.

The function iterates over each bus in the system and populates the vectors `PSpecified`, `QSpecified`, `V`, and `delta` with the corresponding values based on the bus type.

The `busTypes` argument allows specifying the category of bus types to consider. If set to "current" (default), it considers the current bus types. If set to "original", it considers the original bus types.

The resulting vectors are returned as a named tuple containing the initialized vectors and additional information about the bus categories and counts.

```julia
initializeVectors_pu(CDF_DF_List_pu, busTypes="current")
"""
function initializeVectors_pu(CDF_DF_List_pu::Vector{DataFrame};
    busTypes::String="current"):: NamedTuple{(:PSpecified, :QSpecified, :V, :delta, :listOfSlackBuses, :listOfPVBuses, :listOfPQBuses, :listOfNonSlackBuses, :nSlack, :nPV, :nPQ), Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int64}, Vector{Int64}, Vector{Int64}, Vector{Int64}, Int64, Int64, Int64}}
    
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

    if busTypes == "current"
        busType = busData_pu.Type
    elseif busTypes == "original"
        busType = busData_pu.Type_Original
    else
        error("Unknown category of busTypes.")
    end
    
    for i = 1:N
        bus = busData_pu.Bus_Num[i]
        delta[bus] = 0.0000

        if busType[i] == 0
            nPQ += 1
            push!(listOfPQBuses, bus)
            n += 1
            push!(listOfNonSlackBuses, bus)
            V[bus] = 1.0000
        elseif busType[i] == 2
            nPV += 1
            push!(listOfPVBuses, bus)
            n += 1
            push!(listOfNonSlackBuses, bus)
            V[bus] = busData_pu.Desired_V_pu[i]
        elseif busType[i] == 3
            nSlack += 1
            push!(listOfSlackBuses, bus)
            V[bus] = busData_pu.Desired_V_pu[i]
        end

        PSpecified[bus] = busData_pu.Gen_MW[i] - busData_pu.Load_MW[i]
        QSpecified[bus] = busData_pu.Gen_MVAR[i] - busData_pu.Load_MVAR[i]
    end

    listOfSlackBuses = sort!(reshape(listOfSlackBuses[1:nSlack], nSlack))
    listOfNonSlackBuses = sort!(reshape(listOfNonSlackBuses[1:n], n))
    listOfPVBuses = sort!(reshape(listOfPVBuses[1:nPV], nPV))
    listOfPQBuses = sort!(reshape(listOfPQBuses[1:nPQ], nPQ))

    return (PSpecified=PSpecified, QSpecified=QSpecified, V=V, delta=delta, 
    listOfSlackBuses=listOfSlackBuses, listOfPVBuses=listOfPVBuses, 
    listOfPQBuses=listOfPQBuses, listOfNonSlackBuses=listOfNonSlackBuses,
    nSlack=nSlack, nPV=nPV, nPQ=nPQ)
end

"""
    getPVBuses(CDF_Data_List_pu::Vector{DataFrame}; busTypes::String = "current")
    
    Retrieve the list of PV buses from the CDF data.

## Arguments
- `CDF_Data_List_pu::Vector{DataFrame}`: A vector of DataFrames containing the CDF data in per-unit (pu) format.
- `busTypes::String`: Optional. The category of bus types to consider. Default is "current". Possible values are "current" and "original".

## Returns
- `result.listOfPVBuses::Vector{Int64}`: A list of PV buses extracted from the CDF data.
"""
function getPVBuses(CDF_Data_List_pu::Vector{DataFrame}; 
    busTypes::String = "current") :: Vector{Int64}
    result = initializeVectors_pu(CDF_Data_List_pu; busTypes)
    return result.listOfPVBuses
end

"""
    getPQBuses(CDF_Data_List_pu::Vector{DataFrame}; busTypes::String = "current")

Retrieve the list of PQ (Load) buses from the CDF data.

## Arguments
- `CDF_Data_List_pu::Vector{DataFrame}`: A vector of DataFrames containing the CDF data in per-unit (pu) format.
- `busTypes::String`: Optional. The category of bus types to consider. Default is "current". Possible values are "current" and "original".

## Returns
- `result.listOfPQBuses::Vector{Int64}`: A list of PQ (Load) buses extracted from the CDF data.
"""
function getPQBuses(CDF_Data_List_pu::Vector{DataFrame}; busTypes::String = "current") :: Vector{Int64}
    result = initializeVectors_pu(CDF_Data_List_pu; busTypes)
    return result.listOfPQBuses
end

"""
    getNonSlackBuses(CDF_Data_List_pu::Vector{DataFrame}; busTypes::String = "current")

Retrieve the list of non-slack buses (PV and PQ) from the CDF data.

## Arguments
- `CDF_Data_List_pu::Vector{DataFrame}`: A vector of DataFrames containing the CDF data in per-unit (pu) format.
- `busTypes::String`: Optional. The category of bus types to consider. Default is "current". Possible values are "current" and "original".

## Returns
- `result.listOfNonSlackBuses::Vector{Int64}`: A list of non-slack buses (PV and PQ) extracted from the CDF data.
"""
function getNonSlackBuses(CDF_Data_List_pu::Vector{DataFrame}; busTypes::String = "current") :: Vector{Int64}
    result = initializeVectors_pu(CDF_Data_List_pu; busTypes)
    return result.listOfNonSlackBuses
end

"""
    getSlackBuses(CDF_Data_List_pu::Vector{DataFrame}; busTypes::String = "current")

Retrieve the list of slack buses from the CDF data.

## Arguments
- `CDF_Data_List_pu::Vector{DataFrame}`: A vector of DataFrames containing the CDF data in per-unit (pu) format.
- `busTypes::String`: Optional. The category of bus types to consider. Default is "current". Possible values are "current" and "original".

## Returns
- `result.listOfSlackBuses::Vector{Int64}`: A list of slack buses extracted from the CDF data.
"""
function getSlackBuses(CDF_Data_List_pu::Vector{DataFrame}; busTypes::String = "current") :: Vector{Int64}
    result = initializeVectors_pu(CDF_Data_List_pu; busTypes)
    return result.listOfSlackBuses
end

"""
    numSlackBuses(CDF_Data_List_pu::Vector{DataFrame}; busTypes::String = "current")

Retrieve the number of slack buses from the CDF data.

## Arguments
- `CDF_Data_List_pu::Vector{DataFrame}`: A vector of DataFrames containing the CDF data in per-unit (pu) format.
- `busTypes::String`: Optional. The category of bus types to consider. Default is "current". Possible values are "current" and "original".

## Returns
- `result.nSlack::Int64`: The number of slack buses extracted from the CDF data.
"""
function numSlackBuses(CDF_Data_List_pu::Vector{DataFrame}; busTypes::String = "current") :: Int64
    result = initializeVectors_pu(CDF_Data_List_pu; busTypes)
    return result.nSlack
end

"""
    numPVBuses(CDF_Data_List_pu::Vector{DataFrame}; busTypes::String = "current")

Retrieve the number of PV buses from the CDF data.

## Arguments
- `CDF_Data_List_pu::Vector{DataFrame}`: A vector of DataFrames containing the CDF data in per-unit (pu) format.
- `busTypes::String`: Optional. The category of bus types to consider. Default is "current". Possible values are "current" and "original".

## Returns
- `result.nPV::Int64`: The number of PV buses extracted from the CDF data.
"""
function numPVBuses(CDF_Data_List_pu::Vector{DataFrame}; busTypes::String = "current") :: Int64
    result = initializeVectors_pu(CDF_Data_List_pu; busTypes)
    return result.nPV
end

"""
    numPQBuses(CDF_Data_List_pu::Vector{DataFrame}; busTypes::String = "current")

Retrieve the number of PQ (Load) buses from the CDF data.

## Arguments
- `CDF_Data_List_pu::Vector{DataFrame}`: A vector of DataFrames containing the CDF data in per-unit (pu) format.
- `busTypes::String`: Optional. The category of bus types to consider. Default is "current". Possible values are "current" and "original".

## Returns
- `result.nPQ::Int64`: The number of PQ (Load) buses extracted from the CDF data.
"""
function numPQBuses(CDF_Data_List_pu::Vector{DataFrame}; busTypes::String = "current") :: Int64
    result = initializeVectors_pu(CDF_Data_List_pu; busTypes)
    return result.nPQ
end

function get_position(matrix, position::String)
    nrows, ncols = size(matrix)
    positions = Dict(
        "northeast" => (nrows * 0.9, ncols * 0.9),
        "northwest" => (nrows * 0.9, ncols * 0.1),
        "southeast" => (nrows * 0.1, ncols * 0.9),
        "southwest" => (nrows * 0.1, ncols * 0.1)
        # Add more if needed
    )
    return positions[position]
end

