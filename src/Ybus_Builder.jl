# Ybus_Builder.jl
module Ybus_Builder
#This module requires DataFrames, CSV, DelimitedFiles
using DataFrames 
#For some reason, invoking it from the main file
#throws off an error when trying to include/use this module
#saying that it does not recognize 'any' DataFrame.

export initializeVectors_pu
export sortMatrixByBusTypes
export ybusGenerator 
export Create_Ybus_WithoutTaps
export Create_Ybus_WithTaps

"""
    Create_Ybus_WithoutTaps(CDF_DF_List)

Creates Ybus without taps for a power system network.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
'''
'''
# Output
- 'Ybus_WithoutTaps': A complex array of Ybus elements ordered according to bus
type: Slack->PQ->PV.
'''
"""
function Create_Ybus_WithoutTaps(CDF_DF_List_pu)

    # Getting required data from CDF_DF_List
    BusDataCard_DF = CDF_DF_List_pu[2]
    BranchDataCard_DF = CDF_DF_List_pu[3]

    # Getting Size of Ybus
    Size_Ybus = length(BusDataCard_DF.Bus_Num)

    # Initializing Ybus Complex Array
    Ybus_WithoutTaps = Array{Complex{Float64}}(undef, Size_Ybus,Size_Ybus)

    # Computing Ybus Off-Diagonal elements
    for ii in 1:Size_Ybus # Through Rows

        for jj in 1:1:Size_Ybus # Through Columns

            if (ii == jj) # Diagonal Element

                continue

            elseif (ii < jj) # Off-Diagonal elements upper triangle

                # Getting currentBus Numbers from BusDataCard_DF
                Bus1_Num = BusDataCard_DF.Bus_Num[ii]

                Bus2_Num = BusDataCard_DF.Bus_Num[jj]

                # Finding Row in BranchDataCard_DF based on current Bus Numbers
                BranchDataCard_FilterRow = filter(row -> ((row.Tap_Bus_Num == Bus1_Num) && (row.Z_Bus_Num == Bus2_Num)) || ((row.Tap_Bus_Num == Bus2_Num) && (row.Z_Bus_Num == Bus1_Num)), BranchDataCard_Row)

                BranchDataCard_FilterRow_Num = nrow(BranchDataCard_FilterRow)

                if (BranchDataCard_FilterRow_Num == 0) # There is no connection between buses

                    # Filling up Branch Admittance in Ybus_WithoutTaps
                    Ybus_WithoutTaps[ii,jj] = complex(0,0)

                    # Ybus is Symmetrical
                    Ybus_WithoutTaps[jj,ii] = Ybus_WithoutTaps[ii,jj]

                elseif (BranchDataCard_FilterRow_Num > 0) # There is connection between buses

                    # Creating Line Series Admittance
                    Line_SeriesAdmittance = 1/complex(BranchDataCard_FilterRow.R_pu[1],BranchDataCard_FilterRow.X_pu[1])

                    # Filling up Branch Admittance in Ybus_WithoutTaps
                    Ybus_WithoutTaps[ii,jj] = -Line_SeriesAdmittance

                    # Ybus is Symmetrical
                    Ybus_WithoutTaps[jj,ii] = Ybus_WithoutTaps[ii,jj]

                end


            else # Off-Diagonal elements lower triangle

                continue

            end

        end

    end

    # Computing Ybus Diagonal elements
    for ii in 1:Size_Ybus # Through Diagonal Elements Row-wise

        # Getting effect of Off-Diagonal Terms
        OffDiagonal_Terms = 0

        for jj = 1:Size_Ybus # Thorough Columns

            if (ii == jj) # Diagonal Term

                continue

            else

                OffDiagonal_Terms = OffDiagonal_Terms + (-Ybus_WithoutTaps[ii,jj])

            end

        end

        # Getting Effect of Bus Shunt Admittance
        BusAdmittance_Shunt = complex(BusDataCard_DF.G_pu[ii],-BusDataCard_DF.B_pu[ii])

        # Getting Effect of Line Shunt Admittance connected to the Bus
        BusLineAdmittance_Shunt = 0

        Bus_Num = BusDataCard_DF.Bus_Num[ii]

        BranchDataCard_Filter = filter(row -> (row.Tap_Bus_Num == Bus_Num) || (row.Z_Bus_Num == Bus_Num), BranchDataCard_Row)

        BranchDataCard_Filter_Num = nrow(BranchDataCard_Filter)

        if (BranchDataCard_Filter_Num == 0) # Bus not connected to any other bus through a transmission line

            BusLineAdmittance_Shunt = 0

        elseif (BranchDataCard_Filter_Num > 0) # Bus connected to other buses through a transmission lines

            for kk in 1:length(BranchDataCard_Filter.Tap_Bus_Num)

                BusLineAdmittance_Shunt = BusLineAdmittance_Shunt + complex(0,-(BranchDataCard_Filter.B_pu[kk]/2))

            end

        end

        # Total effect oin Ybus Diagonal term
        Ybus_WithoutTaps[ii,ii] = OffDiagonal_Terms + BusAdmittance_Shunt + BusLineAdmittance_Shunt

    end

    # Rearranging Create_Ybus_WithoutTaps in the order Slack->PQ->PV
    Ybus_WithoutTaps_PQ_PV = Ybus_WithoutTaps[1:end-1,1:end-1]

    Ybus_WithoutTaps_Slack1 = Ybus_WithoutTaps[1:end-1,end]

    Ybus_WithoutTaps_Slack2 = Ybus_WithoutTaps[end,1:end-1]

    Ybus_WithoutTaps_Slack3 = Ybus_WithoutTaps[end,end]

    Ybus_WithoutTaps_1 = vcat(Ybus_WithoutTaps_Slack2,Ybus_WithoutTaps_PQ_PV)

    Ybus_WithoutTaps_2 = vcat(Ybus_WithoutTaps_Slack3,Ybus_WithoutTaps_Slack1)

    Ybus_WithoutTaps = hcat(Ybus_WithoutTaps_1,Ybus_WithoutTaps_2)

    return Ybus_WithoutTaps

end

"""
    Create_Ybus_WithoutTaps(CDF_DF_List)

Creates Ybus with taps for a power system network.

'''
# Arguments
- 'Ybus_WithoutTaps': A complex array of Ybus elements ordered according to bus
type: Slack->PQ->PV.
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [BusDataCard_DF, BranchDataCard_DF,
LossZonesCard_DF, InterchangeDataCard_DF, TieLinesDataCard_DF].
'''
'''
# Output
- 'Ybus_WithTaps': A complex array of Ybus elements ordered according to bus
type: Slack->PV->PQ.
'''
"""
function Create_Ybus_WithTaps(Ybus_WithoutTaps,CDF_DF_List_pu)

    # Getting required data from CDF_DF_List
    BusDataCard_DF = CDF_DF_List_pu[2]

    BranchDataCard_DF = CDF_DF_List_pu[3]

    # Getting Last Row number of BusDataCard_DF to locate Slack Bus Row number
    SlackBus_RowNumber = length(BusDataCard_DF.Bus_Num)

    # Initializing Ybus_WithTaps
    Ybus_WithTaps = Ybus_WithoutTaps

    # Getting Subset of BranchDataCard_DFfor lines with Tap Changing Transformers
    BranchDataCard_Filter = filter(row -> ((row.Transformer_t != 0) || (row.Transformer_ps != 0), BranchDataCard_DF))

    BranchDataCard_Filter_Num = nrow(BranchDataCard_Filter)

    if (BranchDataCard_Filter_Num == 0) # No Tap Changing Transformers present

        Ybus_WithTaps = Ybus_WithoutTaps

    elseif (BranchDataCard_Filter_Num > 0) # Tap Changing Transformers present

        for ii in 1:BranchDataCard_Filter_Num

            # Creating Tap value 'a'
            a = BranchDataCard_Filter.Transformer_t[ii] * cis(deg2rad(BranchDataCard_Filter.Transformer_ps[ii]))

            # Getting Bus Numbers 'i': Z_Bus_Num (Impedance Side) , 'j': Tap_Bus_Num (non-unity tap Side)
            Bus_Num_i = BranchDataCard_Filter.Z_Bus_Num[ii]

            Bus_Num_j = BranchDataCard_Filter.Tap_Bus_Num[ii]

            # Getting associated 'Bus_i_Index' and 'Bus_j_Index' from BusDataCard_DF to access correct location within Ybus_WithoutTaps
            for jj in 1:SlackBus_RowNumber

                if (Bus_Num_i == BusDataCard_DF.Bus_Num[jj])

                    if (jj == SlackBus_RowNumber)

                        Bus_i_Index = 1

                    else

                        Bus_i_Index = jj+1

                    end

                elseif (Bus_Num_j == BusDataCard_DF.Bus_Num[jj])

                    if (jj == SlackBus_RowNumber)

                        Bus_j_Index = 1

                    else

                        Bus_j_Index = jj+1

                    end

                else

                    continue

                end

            end

            # Changing the [Bus_i_Index, Bus_j_Index] in Ybus_WithTaps based on 'a'

            # Changing [Bus_i_Index, Bus_i_Index]
            Ybus_WithTaps[Bus_i_Index, Bus_i_Index] = Ybus_WithoutTaps[Bus_i_Index, Bus_i_Index]

            # Changing [Bus_i_Index, Bus_j_Index]
            Ybus_WithTaps[Bus_i_Index, Bus_j_Index] = Ybus_WithoutTaps[Bus_i_Index, Bus_j_Index]/a

            # Changing [Bus_j_Index, Bus_i_Index]
            Ybus_WithTaps[Bus_j_Index, Bus_i_Index] = Ybus_WithoutTaps[Bus_j_Index, Bus_i_Index]/conj(a)

            # Changing [Bus_j_Index, Bus_j_Index]
            Ybus_WithTaps[Bus_j_Index, Bus_j_Index] = (Ybus_WithoutTaps[Bus_j_Index, Bus_j_Index]) - (-Ybus_WithoutTaps[Bus_j_Index, Bus_i_Index]) + (-Ybus_WithoutTaps[Bus_j_Index, Bus_i_Index]/abs2(a))

        end

    end

    return Ybus_WithTaps

end

"""
initializeVectors_pu

This function initializes various vectors and variables based on the input `CDF_DF_List_pu`.

# Parameters
- `CDF_DF_List_pu`: A list containing the data for power system buses in the per unit (pu) system. It should be a two-dimensional array-like object where each row represents a bus and each column represents a specific attribute of the bus.

# Returns
- A list containing the following vectors and variables:
  - `PSpecified`: An array of size `N` representing the specified active power injections (in MW) for each bus.
  - `QSpecified`: An array of size `N` representing the specified reactive power injections (in MVAR) for each bus.
  - `V`: An array of size `N` representing the initial voltage magnitudes (in pu) for each bus.
  - `delta`: An array of size `N` representing the initial voltage phase angles (in radians) for each bus.
  - `listOfSlackBuses`: An array containing the bus numbers of the slack buses in the system.
  - `listOfPVBuses`: An array containing the bus numbers of the PV buses in the system.
  - `listOfPQBuses`: An array containing the bus numbers of the PQ buses in the system.
  - `listOfNonSlackBuses`: An array containing the bus numbers of the non-slack buses in the system.
  - `nSlack`: The number of slack buses in the system.
  - `nPV`: The number of PV buses in the system.
  - `nPQ`: The number of PQ buses in the system.

# Description
The `initializeVectors_pu` function takes the input `CDF_DF_List_pu`, which represents the data for power system buses in the per unit (pu) system, and initializes various vectors and variables based on this data. The function iterates over each bus in the system and performs the following steps:

1. Retrieves the bus data for the current bus from `busData_pu`.
2. Sets the initial value of `delta` (voltage phase angle) for the current bus to zero.
3. Checks the type of the current bus:
   - If the bus type is 0, it is a PQ bus.
     - Increments the count of PQ buses (`nPQ`).
     - Adds the bus number to the `listOfPQBuses`.
     - Adds the bus number to the `listOfNonSlackBuses`.
     - Sets the initial value of voltage magnitude `V` for the bus to 1.0000 pu.
   - If the bus type is 2, it is a PV bus.
     - Increments the count of PV buses (`nPV`).
     - Adds the bus number to the `listOfPVBuses`.
     - Adds the bus number to the `listOfNonSlackBuses`.
     - Sets the initial value of voltage magnitude `V` for the bus to the desired voltage magnitude specified in `busData_pu`.
   - If the bus type is 3, it is a slack bus.
     - Increments the count of slack buses (`nSlack`).
     - Adds the bus number to the `listOfSlackBuses`.
     - Sets the initial value of voltage magnitude `V` for the bus to the desired voltage magnitude specified in `busData_pu`.
4. Calculates the specified active power injection `PSpecified` for the bus by subtracting the load active power from the generator active power specified in `busData_pu`.
5. Calculates the specified reactive power injection `QSpecified` for the bus by subtracting the load reactive power from the generator reactive power specified in `busData_pu`.

After iterating over all buses, the function performs the following additional operations:

- Reshapes `listOfSlackBuses`, `listOfPVBuses`, `listOfPQBuses`, and `listOfNonSlackBuses` to remove any unused elements.
- Returns all the calculated vectors and variables as a list.

Please note that the documentation assumes some familiarity with power system analysis terminology and concepts.
"""
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

"""
sortMatrixByBusTypes

This function sorts the given `ybus` matrix and row names based on bus types using the `initializeVectors_pu` function.

# Parameters
- `CDF_DF_List_pu`: A list containing the data for power system buses in the per unit (pu) system. It should be a two-dimensional array-like object where each row represents a bus and each column represents a specific attribute of the bus.
- `ybus`: The admittance matrix of the power system.

# Returns
- `ybusByTypes`: The sorted admittance matrix `ybus` based on bus types.
- `rowNamesByTypes`: The sorted row names corresponding to `ybusByTypes`.

# Description
The `sortMatrixByBusTypes` function uses the `initializeVectors_pu` function to obtain the lists of slack buses, PV buses, and PQ buses. It then combines these lists into a new order, which represents the desired sorting order of buses in the `ybus` matrix.

The function creates a new `ybusByTypes` matrix by reordering the rows and columns of `ybus` according to the new order of bus types. The `rowNamesByTypes` is also updated to match the new order of bus types.

Finally, the sorted `ybusByTypes` matrix and `rowNamesByTypes` are returned as the output of the function.

Please note that the `initializeVectors_pu` function is assumed to be defined and implemented separately.
"""
function sortMatrixByBusTypes(CDF_DF_List_pu, ybus)
    # Call initializeVectors_pu to obtain bus type information
    outputs  = initializeVectors_pu(CDF_DF_List_pu)
    listOfSlackBuses, listOfPVBuses, listOfPQBuses = outputs[5], outputs[6], outputs[7]
    # Create a new order based on bus types
    newOrder = vcat(listOfSlackBuses, listOfPVBuses, listOfPQBuses)
    # Reorder the ybus matrix and row names according to the new order
    ybusByTypes = ybus[newOrder, newOrder]

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
ybus, BMatrix, b, A, branchNames, E = ybusGenerator(CDF_DF_List_pu, verbose=true, saveTables=true)
"""
function ybusGenerator(CDF_DF_List_pu::Vector{DataFrame};
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

end
