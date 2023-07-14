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
    initializeVectors_pu(CDF_DF_List_pu::Vector{DataFrame})

This function initializes various vectors and variables based on the input `CDF_DF_List_pu`.

# Arguments
- `CDF_DF_List_pu::Vector{DataFrame}`: A vector of DataFrames containing the CDF (Common Data Format) data in per-unit (pu) format.

# Returns
- NamedTuple containing the following fields:
- `PSpecified::Vector{Float64}`: An array of size `N` representing the specified active power injections (in MW) for each bus.
- `QSpecified::Vector{Float64}`: An array of size `N` representing the specified reactive power injections (in MVAR) for each bus.
- `V::Vector{Float64}`: An array of size `N` representing the initial voltage magnitudes (in pu) for each bus.
- `delta::Vector{Float64}`: An array of size `N` representing the initial voltage phase angles (in radians) for each bus.
- `listOfSlackBuses::Vector{Int64}`: An array containing the bus numbers of the slack buses in the system.
- `listOfPVBuses::Vector{Int64}`: An array containing the bus numbers of the PV buses in the system.
- `listOfPQBuses::Vector{Int64}`: An array containing the bus numbers of the PQ buses in the system.
- `listOfNonSlackBuses::Vector{Int64}`: An array containing the bus numbers of the non-slack buses in the system.
- `nSlack::Int64`: The number of slack buses in the system.
- `nPV::Int64`: The number of PV buses in the system.
- `nPQ::Int64`: The number of PQ buses in the system.

# Description
This function takes the input `CDF_DF_List_pu`, which represents the data for power system buses in the per unit (pu) system, and initializes various vectors and variables based on this data. The function iterates over each bus in the system and performs the following steps:

1. Retrieves the bus data for the current bus from `busData_pu`.
2. Sets the initial value of `delta` (voltage phase angle) for the current bus to zero.
3. Checks the type of the current bus:
- If the bus type is 0, it is a PQ bus.
    - Increments the count of PQ buses (`nPQ`).
    - Adds the bus number to `listOfPQBuses`.
    - Adds the bus number to `listOfNonSlackBuses`.
    - Sets the initial value of voltage magnitude `V` for the bus to 1.0000 pu.
- If the bus type is 2, it is a PV bus.
    - Increments the count of PV buses (`nPV`).
    - Adds the bus number to `listOfPVBuses`.
    - Adds the bus number to `listOfNonSlackBuses`.
    - Sets the initial value of voltage magnitude `V` for the bus to the desired voltage magnitude specified in `busData_pu`.
- If the bus type is 3, it is a slack bus.
    - Increments the count of slack buses (`nSlack`).
    - Adds the bus number to `listOfSlackBuses`.
    - Sets the initial value of voltage magnitude `V` for the bus to the desired voltage magnitude specified in `busData_pu`.
4. Calculates the specified active power injection `PSpecified` for the bus by subtracting the load active power from the generator active power specified in `busData_pu`.
5. Calculates the specified reactive power injection `QSpecified` for the bus by subtracting the load reactive power from the generator reactive power specified in `busData_pu`.

After iterating over all buses, the function performs the following additional operations:

- Reshapes `listOfSlackBuses`, `listOfPVBuses`, `listOfPQBuses`, and `listOfNonSlackBuses` to remove any unused elements.
- Returns all the calculated vectors and variables as a NamedTuple.

Please note that the documentation assumes some familiarity with power system analysis terminology and concepts.
"""
function initializeVectors_pu(CDF_DF_List_pu::Vector{DataFrame})
    
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
    getPVBuses(CDF_Data_List_pu::Vector{DataFrame})

This function retrieves the list of PV buses from the CDF_Data_List_pu. It initializes the vectors using the `initializeVectors_pu` function and returns the listOfPVBuses from the result.

Arguments:
- CDF_Data_List_pu::Vector{DataFrame}: A vector of DataFrames containing the CDF (Common Data Format) data in per-unit (pu) format.

Returns:
- result.listOfPVBuses: A list of PV buses extracted from the CDF data.
"""
function getPVBuses(CDF_Data_List_pu::Vector{DataFrame})
    result = initializeVectors_pu(CDF_Data_List_pu)
    return result.listOfPVBuses
end

"""
    getPQBuses(CDF_Data_List_pu::Vector{DataFrame})

This function retrieves the list of PQ (Load) buses from the CDF_Data_List_pu. It initializes the vectors using the `initializeVectors_pu` function and returns the listOfPQBuses from the result.

Arguments:
- CDF_Data_List_pu::Vector{DataFrame}: A vector of DataFrames containing the CDF (Common Data Format) data in per-unit (pu) format.

Returns:
- result.listOfPQBuses: A list of PQ (Load) buses extracted from the CDF data.
"""
function getPQBuses(CDF_Data_List_pu::Vector{DataFrame})
    result = initializeVectors_pu(CDF_Data_List_pu)
    return result.listOfPQBuses
end

"""
    getNonSlackBuses(CDF_Data_List_pu::Vector{DataFrame})

This function retrieves the list of non-slack buses (PV and PQ) from the CDF_Data_List_pu. It initializes the vectors using the `initializeVectors_pu` function and returns the listOfNonSlackBuses from the result.

Arguments:
- CDF_Data_List_pu::Vector{DataFrame}: A vector of DataFrames containing the CDF (Common Data Format) data in per-unit (pu) format.

Returns:
- result.listOfNonSlackBuses: A list of non-slack buses (PV and PQ) extracted from the CDF data.
"""
function getNonSlackBuses(CDF_Data_List_pu::Vector{DataFrame})
    result = initializeVectors_pu(CDF_Data_List_pu)
    return result.listOfNonSlackBuses
end

"""
    getSlackBuses(CDF_Data_List_pu::Vector{DataFrame})

This function retrieves the list of slack buses from the CDF_Data_List_pu. It initializes the vectors using the `initializeVectors_pu` function and returns the listOfSlackBuses from the result.

Arguments:
- CDF_Data_List_pu::Vector{DataFrame}: A vector of DataFrames containing the CDF (Common Data Format) data in per-unit (pu) format.

Returns:
- result.listOfSlackBuses: A list of slack buses extracted from the CDF data.
"""
function getSlackBuses(CDF_Data_List_pu::Vector{DataFrame})
    result = initializeVectors_pu(CDF_Data_List_pu)
    return result.listOfSlackBuses
end

"""
    numSlackBuses(CDF_Data_List_pu::Vector{DataFrame})

This function retrieves the number of slack buses from the CDF_Data_List_pu. It initializes the vectors using the `initializeVectors_pu` function and returns the nSlack from the result.

Arguments:
- CDF_Data_List_pu::Vector{DataFrame}: A vector of DataFrames containing the CDF (Common Data Format) data in per-unit (pu) format.

Returns:
- result.nSlack: The number of slack buses extracted from the CDF data.
"""
function numSlackBuses(CDF_Data_List_pu::Vector{DataFrame})
    result = initializeVectors_pu(CDF_Data_List_pu)
    return result.nSlack
end

"""
    numPVBuses(CDF_Data_List_pu::Vector{DataFrame})

This function retrieves the number of PV buses from the CDF_Data_List_pu. It initializes the vectors using the `initializeVectors_pu` function and returns the nPV from the result.

Arguments:
- CDF_Data_List_pu::Vector{DataFrame}: A vector of DataFrames containing the CDF (Common Data Format) data in per-unit (pu) format.

Returns:
- result.nPV: The number of PV buses extracted from the CDF data.
"""
function numPVBuses(CDF_Data_List_pu::Vector{DataFrame})
    result = initializeVectors_pu(CDF_Data_List_pu)
    return result.nPV
end

"""
    numPQBuses(CDF_Data_List_pu::Vector{DataFrame})

This function retrieves the number of PQ (Load) buses from the CDF_Data_List_pu. It initializes the vectors using the `initializeVectors_pu` function and returns the nPQ from the result.

Arguments:
- CDF_Data_List_pu::Vector{DataFrame}: A vector of DataFrames containing the CDF (Common Data Format) data in per-unit (pu) format.

Returns:
- result.nPQ: The number of PQ (Load) buses extracted from the CDF data.
"""
function numPQBuses(CDF_Data_List_pu::Vector{DataFrame})
    result = initializeVectors_pu(CDF_Data_List_pu)
    return result.nPQ
end
