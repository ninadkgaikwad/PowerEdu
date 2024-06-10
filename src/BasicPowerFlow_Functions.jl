# BasicPowerFlow_Functions.jl

using LaTeXStrings

"""
    Create_Initial_SolutionVector_NR(CDF_DF_List_pu)

Creates Initial Solution Vector for the Newton-Raphson Method.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
'''
'''
# Output
- 'Initial_SolutionVector_NR': Flat start Voltage and Angle at each bus ordered according to bus
type: PQ->PV.
'''
"""
function Create_Initial_SolutionVector_NR(CDF_DF_List_pu)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type == 0) || (row.Type == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

        # Creating Initial_SolutionVector_NR Length
        Initial_SolutionVector_NR_Len = 2 * N_PQ_Bus + N_PV_Bus

        # Initializing Initial_SolutionVector_NR
        Initial_SolutionVector_NR = Array{Float64}(undef,Initial_SolutionVector_NR_Len,1)

        # Creating Initial_SolutionVector_NR Delta part for both PQ and PV Buses
        Initial_SolutionVector_NR[1:(N_PQ_Bus+N_PV_Bus),1] = zeros((N_PQ_Bus+N_PV_Bus),1)

        # Creating Initial_SolutionVector_NR V part for PQ Buses
        Initial_SolutionVector_NR[(N_PQ_Bus+N_PV_Bus+1):end,1] = ones(N_PQ_Bus,1)

        return Initial_SolutionVector_NR

end

"""
    Create_SolutionVector_VDelta_NR(CDF_DF_List_pu, SolutionVector_NR)

Creates separate V and Delta vectors from the solution vector computed by the NR Method and data in the CDF file.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
- 'SolutionVector_NR': Voltage and Angle at each bus ordered according to bus
type: PQ->PV.
'''
'''
# Output
- 'SolutionVector_V': Voltage at each bus ordered according to bus
type: Slack->PQ->PV.
- 'SolutionVector_Delta': Angle at each bus ordered according to bus
type: Slack->PQ->PV.
'''
"""
function Create_SolutionVector_VDelta_NR(CDF_DF_List_pu, SolutionVector_NR)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type == 0) || (row.Type == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

        # Initializing SolutionVector_V, SolutionVector_Delta
        SolutionVector_V = Array{Float64}(undef, N_Bus,1)
        SolutionVector_Delta = Array{Float64}(undef, N_Bus,1)

        # Creating SolutionVector_V, SolutionVector_Delta
        SolutionVector_Delta =  vcat([0.0], SolutionVector_NR[1:(N_Bus-1),1])
        
        SolutionVector_V[1:(N_PQ_Bus+1),1] =  vcat([BusDataCard_DF.Final_V_pu_Original[end]], SolutionVector_NR[(N_Bus-1)+1:end,1])
        
        for ii in N_PQ_Bus+1+1:N_Bus # For each current PV Bus

                if (BusDataCard_DF.Type_Original[ii-1] == 1) # PQ Converted to PV

                        SolutionVector_V[ii,1] = BusDataCard_DF.MVAR_V_Limit[ii-1]

                elseif (BusDataCard_DF.Type_Original[ii-1] == 2) # Originally PV

                        SolutionVector_V[ii,1] = BusDataCard_DF.Final_V_pu_Original[ii-1]

                end

        end

        return SolutionVector_V, SolutionVector_Delta

end

"""
    Create_SolutionVector_NR(CDF_DF_List_pu, SolutionVector_V, SolutionVector_Delta)

Creates Initial Solution Vector for the Newton-Raphson Method.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
- 'SolutionVector_V': Voltage at each bus ordered according to bus
type: Slack->PQ->PV.
- 'SolutionVector_Delta': Angle at each bus ordered according to bus
type: Slack->PQ->PV.
'''
'''
# Output
- 'SolutionVector_NR': Voltage and Angle at each bus ordered according to bus
type: PQ->PV.
'''
"""
function Create_SolutionVector_NR(CDF_DF_List_pu, SolutionVector_V, SolutionVector_Delta)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type == 0) || (row.Type == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

        # Creating SolutionVector_NR Length
        SolutionVector_NR_Len = 2*N_PQ_Bus+N_PV_Bus

        # Initializing SolutionVector_NR
        SolutionVector_NR = Array{Float64}(undef, SolutionVector_NR_Len ,1)

        # Creating SolutionVector_NR Delta part for both PQ and PV Buses
        SolutionVector_NR[1:(N_PQ_Bus+N_PV_Bus),1] = SolutionVector_Delta[2:end,1]

        # Creating Initial_SolutionVector_NR V part for PQ Buses
        SolutionVector_NR = vcat(SolutionVector_NR[1:(N_PQ_Bus+N_PV_Bus),1], SolutionVector_V[2:N_PQ_Bus+1,1])

        return SolutionVector_NR

end

"""
    Compute_PQ_BusArray(Ybus, SolutionVector_V, SolutionVector_Delta)

Computes P and Q at each bus of a power system network.

'''
# Arguments
- 'Ybus': A complex array of Ybus elements ordered according to bus
type: Slack->PQ->PV.
- 'SolutionVector_V': Voltage at each bus ordered according to bus
type: Slack->PQ->PV.
- 'SolutionVector_Delta': Voltage Angle at each bus ordered according to bus
type: Slack->PQ->PV.
'''
'''
# Output
- 'PQ_BusArray': An array (N*2) for P and Q vectors ordered according to bus
type: Slack->PQ->PV.
'''
"""
function Compute_PQ_BusArray(Ybus, SolutionVector_V, SolutionVector_Delta)
#     @show Ybus
    #@show SolutionVector_Delta
    #@show typeof(SolutionVector_Delta)
    #@show size(SolutionVector_Delta)
    #@show SolutionVector_V
    #@show typeof(SolutionVector_V)
    #@show size(SolutionVector_Delta)
    # Initializing PQ_BusArray
    PQ_BusArray = Array{Float64}(undef, length(SolutionVector_V), 2)

    # Computing P-Q Values
    for ii in 1:length(SolutionVector_V)

            # Initializing P_Bus and Q_Bus
            P_Bus = 0
            Q_Bus = 0

            # Computing P_Bus  and Q_Bus
            for jj in 1:length(SolutionVector_V)

                    # Computing P_Bus
                    P_Bus = P_Bus + (SolutionVector_V[ii]*SolutionVector_V[jj]*abs(Ybus[ii,jj])*cos(angle(Ybus[ii,jj])+deg2rad(SolutionVector_Delta[jj])-deg2rad(SolutionVector_Delta[ii])))

                    # Computing Q_Bus
                    Q_Bus = Q_Bus - (SolutionVector_V[ii]*SolutionVector_V[jj]*abs(Ybus[ii,jj])*sin(angle(Ybus[ii,jj])+deg2rad(SolutionVector_Delta[jj])-deg2rad(SolutionVector_Delta[ii])))

            end

            # Filling up PQ_BusArray
            PQ_BusArray[ii, 1:2] = [P_Bus, Q_Bus]

    end

    # Addressing Machine Precision Problem
    for ii in 1:size(PQ_BusArray)[1]

        if (abs(PQ_BusArray[ii,1]) < 1e-12)
                PQ_BusArray[ii,1] = 0
        end

        if (abs(PQ_BusArray[ii,2]) < 1e-12)
                PQ_BusArray[ii,2] = 0
        end

    end

    return PQ_BusArray

end

"""
    Compute_PQ_MismatchVector(CDF_DF_List_pu, PQ_BusArray, SolutionVector_V,
    NR_Type)

Computes P and Q mistmatch at each bus of a power system network for
Newton-Raphson Method.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
- 'PQ_BusArray': An array (N*2) for P and Q vectors ordered according to bus
type: Slack->PQ->PV.
- 'SolutionVector_V': Voltage at each bus ordered according to bus
type: Slack->PQ->PV.
- 'NR_Type': 1 -> Full Newton-Raphson, 2-> Decoupled Newton-Raphson,
3 -> Fast Decoupled Newton-Raphson
'''
'''
# Output
- 'PQ_MismatchVector': An complex array of P-Q Mistmatch elements ordered
according to bus type: PQ->PV.
'''
"""
function Compute_PQ_MismatchVector(CDF_DF_List_pu, PQ_BusArray, SolutionVector_V, NR_Type)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type == 0) || (row.Type == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

        # Length of PQ_MismatchVector
        Len_PQ_MismatchVector = 2*N_PQ_Bus+N_PV_Bus

        # Initializing PQ_MismatchVector
        PQ_MismatchVector = Array{Float64}(undef, Len_PQ_MismatchVector, 1)

        # Initializing PQ_PV_Counter
        PQ_PV_Counter = 0

        # Computing P-Q Mismatch
        for ii in 1:Len_PQ_MismatchVector

                if (ii <= (N_Bus-1)) # Excluding Slack Bus


                        if ((NR_Type == 1) || (NR_Type == 2)) # Full/Decoupled NR

                                # Computing P Mismatch
                                PQ_MismatchVector[ii] = (BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii]) - PQ_BusArray[ii+1,1]

                        elseif (NR_Type == 3) # Fast-Decoupled NR

                                # Computing P Mismatch
                                PQ_MismatchVector[ii] = ((BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii]) - PQ_BusArray[ii+1,1])/(SolutionVector_V[ii+1])

                        end

                else

                        # Incrementing PQ_PV_Counter
                        PQ_PV_Counter = PQ_PV_Counter +1

                        if ((NR_Type == 1) || (NR_Type == 2)) # Full/Decoupled NR

                                if ((BusDataCard_DF.Type[PQ_PV_Counter] != BusDataCard_DF.Type_Original[PQ_PV_Counter]) && (BusDataCard_DF.MVAR_V_Limit[PQ_PV_Counter]  != -9999)) # Q Limit Violation

                                        # Computing Q Mismatch
                                        # PQ_MismatchVector[ii] = (BusDataCard_DF.Gen_MVAR[PQ_PV_Counter] - BusDataCard_DF.Load_MVAR[PQ_PV_Counter] )- BusDataCard_DF.MVAR_V_Limit[PQ_PV_Counter]

                                        PQ_MismatchVector[ii] = BusDataCard_DF.MVAR_V_Limit[PQ_PV_Counter] - PQ_BusArray[PQ_PV_Counter,2]

                                else # No Q Limit Violation

                                        # Computing Q Mismatch
                                        PQ_MismatchVector[ii] = (BusDataCard_DF.Gen_MVAR[PQ_PV_Counter] - BusDataCard_DF.Load_MVAR[PQ_PV_Counter] )- PQ_BusArray[PQ_PV_Counter+1,2]

                                end

                        elseif (NR_Type == 3) # Fast-Decoupled NR

                                if ((BusDataCard_DF.Type[PQ_PV_Counter] != BusDataCard_DF.Type_Original[PQ_PV_Counter]) && (BusDataCard_DF.MVAR_V_Limit[PQ_PV_Counter]  != -9999)) # Q Limit Violation

                                        # Computing Q Mismatch
                                        # PQ_MismatchVector[ii] = ((BusDataCard_DF.Gen_MVAR[PQ_PV_Counter] - BusDataCard_DF.Load_MVAR[PQ_PV_Counter]) - BusDataCard_DF.MVAR_V_Limit[PQ_PV_Counter]) / (SolutionVector_V[PQ_PV_Counter+1])

                                        PQ_MismatchVector[ii] = (BusDataCard_DF.MVAR_V_Limit[PQ_PV_Counter] - PQ_BusArray[PQ_PV_Counter+1,2]) / (SolutionVector_V[PQ_PV_Counter+1])

                                else # No Q Limit Violation

                                        # Computing Q Mismatch
                                        PQ_MismatchVector[ii] = ((BusDataCard_DF.Gen_MVAR[PQ_PV_Counter] - BusDataCard_DF.Load_MVAR[PQ_PV_Counter] )- PQ_BusArray[PQ_PV_Counter+1,2] )/(SolutionVector_V[PQ_PV_Counter+1])

                                end

                        end
                end

        end

        return PQ_MismatchVector

end

function Compute_PQ_MismatchVector1(CDF_DF_List_pu, PQ_BusArray, SolutionVector_V, NR_Type)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type == 0) || (row.Type == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

        # Length of PQ_MismatchVector
        Len_PQ_MismatchVector = 2*N_PQ_Bus+N_PV_Bus

        # Initializing PQ_MismatchVector
        PQ_MismatchVector = Array{Float64}(undef, Len_PQ_MismatchVector, 1)

        # Initializing PQ_PV_Counter
        PQ_PV_Counter = 0

        # Computing P-Q Mismatch
        for ii in 1:Len_PQ_MismatchVector

                if (ii <= (N_Bus-1)) # Excluding Slack Bus


                        if ((NR_Type == 1) || (NR_Type == 2)) # Full/Decoupled NR

                                # Computing P Mismatch
                                PQ_MismatchVector[ii] = (BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii]) - PQ_BusArray[ii+1,1]

                        elseif (NR_Type == 3) # Fast-Decoupled NR

                                # Computing P Mismatch
                                PQ_MismatchVector[ii] = ((BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii]) - PQ_BusArray[ii+1,1])/(SolutionVector_V[ii+1])

                        end

                else

                        # Incrementing PQ_PV_Counter
                        PQ_PV_Counter = PQ_PV_Counter +1

                        if ((NR_Type == 1) || (NR_Type == 2)) # Full/Decoupled NR

                                

                                        # Computing Q Mismatch
                                        PQ_MismatchVector[ii] = (BusDataCard_DF.Gen_MVAR[PQ_PV_Counter] - BusDataCard_DF.Load_MVAR[PQ_PV_Counter] )- PQ_BusArray[PQ_PV_Counter+1,2]

                                

                        elseif (NR_Type == 3) # Fast-Decoupled NR

                                

                                        # Computing Q Mismatch
                                        PQ_MismatchVector[ii] = ((BusDataCard_DF.Gen_MVAR[PQ_PV_Counter] - BusDataCard_DF.Load_MVAR[PQ_PV_Counter] )- PQ_BusArray[PQ_PV_Counter+1,2] )/(SolutionVector_V[PQ_PV_Counter+1])

                                

                        end
                end

        end

        return PQ_MismatchVector

end

"""
    PQ_PV_Bus_Check_Modify(CDF_DF_List_pu, NR_Type, SolutionVector_NR)

Checks for and modifies PQ-PV buses according to their V/MVAR limits.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
- 'Ybus': A complex array of Ybus elements ordered according to bus
type: Slack->PQ->PV.
- 'Ybus_Type': 1 -> Without Taps , 2 -> With Taps
- 'NR_Type': 1 -> Full Newton-Raphson, 2-> Decoupled Newton-Raphson,
3 -> Fast Decoupled Newton-Raphson
- 'SolutionVector_NR': Voltage and Angle at each bus ordered according to bus
type: PQ->PV.
- 'BusSwitching': true -> Bus Switching Employed, false -> Bus Switching not Employed
'''
'''
# Output
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF]. Changed with changes in Bus Types
- 'Ybus': A complex array of Ybus elements ordered according to bus
type: Slack->PQ->PV. Changed with changes in Bus Types
- 'PQ_BusArray': An array (N*2) for P and Q vectors ordered according to bus
type: Slack->PQ->PV. Changed with changes in Bus Type
- 'PQ_MismatchVector': An complex array of P-Q Mistmatch elements ordered
according to bus type: PQ->PV.
- 'SolutionVector_V': Voltage at each bus ordered according to bus
type: Slack->PQ->PV.
- 'SolutionVector_Delta': Angle at each bus ordered according to bus
type: Slack->PQ->PV.
'''
"""
function PQ_PV_Bus_Check_Modify(CDF_DF_List_pu, Ybus, Ybus_Type, NR_Type, SolutionVector_NR, BusSwitching)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type == 0) || (row.Type == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

        # Creating SolutionVector_V and SolutionVector_Delta
        SolutionVector_V, SolutionVector_Delta = Create_SolutionVector_VDelta_NR(CDF_DF_List_pu, SolutionVector_NR)

        # Updating Final_V_pu and Final_A_deg
        for ii in 1:N_Bus

                if (ii == N_Bus)

                        BusDataCard_DF.Final_V_pu[ii] = SolutionVector_V[1,1]

                        BusDataCard_DF.Final_A_deg[ii] = SolutionVector_Delta[1,1]

                else

                        BusDataCard_DF.Final_V_pu[ii] = SolutionVector_V[ii+1,1]

                        BusDataCard_DF.Final_A_deg[ii] = SolutionVector_Delta[ii+1,1]

                end

        end
        

        if (BusSwitching == true) # Employ Bus Switching

                #  Compute P-Q at Buses
                PQ_BusArray = Compute_PQ_BusArray(Ybus, SolutionVector_V, SolutionVector_Delta)

                # Checking and Modifying Bus Types

                for ii in 1:(N_Bus-1)

                        if (BusDataCard_DF.Type_Original[ii] == 1) # Hold MVAR generation within voltage limits, (PQ)

                                # Getting V_Bus
                                V_Bus = SolutionVector_V[ii+1,1]

                                # Getting Q_set/Q_i
                                Q_set = BusDataCard_DF.Gen_MVAR[ii] - BusDataCard_DF.Load_MVAR[ii]
                                Q_i = PQ_BusArray[ii+1,2]

                                # Getting VAR Limits
                                Max_V = BusDataCard_DF.Max_MVAR_V[ii]
                                Min_V = BusDataCard_DF.Min_MVAR_V[ii]

                                # Checking V Limit Violation
                                if (BusDataCard_DF.Type[ii] == 1) # Currently PQ Bus is PQ

                                        if (V_Bus <= Min_V)

                                                # Switiching PQ to PV
                                                BusDataCard_DF.Type[ii] = 2

                                                # Change MVAR_V_Limit to V_Min
                                                BusDataCard_DF.MVAR_V_Limit[ii] = Min_V

                                                BusDataCard_DF.Final_V_pu[ii] = Min_V


                                        elseif (V_Bus >= Max_V)

                                                # Switiching PQ to PV
                                                BusDataCard_DF.Type[ii] = 2

                                                # Change MVAR_V_Limit to V_Max
                                                BusDataCard_DF.MVAR_V_Limit[ii] = Max_V

                                                BusDataCard_DF.Final_V_pu[ii] = Max_V

                                        else

                                                # Keep the current Type to Original Type
                                                BusDataCard_DF.Type[ii] = 1

                                                # Reinitialize MVAR_V_Limit
                                                BusDataCard_DF.MVAR_V_Limit[ii] = -9999

                                                BusDataCard_DF.Final_V_pu[ii] = BusDataCard_DF.Final_V_pu[ii]

                                        end

                                elseif (BusDataCard_DF.Type[ii] == 2) # Currently PQ Bus is PV

                                        if (BusDataCard_DF.MVAR_V_Limit[ii] == Min_V)

                                                if (Q_i >= Q_set)

                                                        # Switiching PQ to PV
                                                        BusDataCard_DF.Type[ii] = 2

                                                        # Change MVAR_V_Limit to V_Min
                                                        BusDataCard_DF.MVAR_V_Limit[ii] = Min_V

                                                        BusDataCard_DF.Final_V_pu[ii] = Min_V

                                                elseif (Q_i < Q_set)

                                                        if (V_Bus <= Min_V)

                                                                # Switiching PQ to PV
                                                                BusDataCard_DF.Type[ii] = 2

                                                                # Change MVAR_V_Limit to V_Min
                                                                BusDataCard_DF.MVAR_V_Limit[ii] = Min_V

                                                                BusDataCard_DF.Final_V_pu[ii] = Min_V

                                                        elseif (V_Bus >= Max_V)

                                                                # Switiching PV to PQ
                                                                BusDataCard_DF.Type[ii] = 2

                                                                # Change MVAR_V_Limit to V_Max
                                                                BusDataCard_DF.MVAR_V_Limit[ii] = Max_V

                                                                BusDataCard_DF.Final_V_pu[ii] = Max_V

                                                        else

                                                                # Keep the current Type to Original Type
                                                                BusDataCard_DF.Type[ii] = 1

                                                                # Reinitialize MVAR_V_Limit
                                                                BusDataCard_DF.MVAR_V_Limit[ii] = -9999

                                                                BusDataCard_DF.Final_V_pu[ii] = BusDataCard_DF.Final_V_pu[ii]

                                                        end

                                                end

                                        elseif (BusDataCard_DF.MVAR_V_Limit[ii] == Max_V)

                                                if (Q_i <= Q_set)

                                                        # Switiching PQ to PV
                                                        BusDataCard_DF.Type[ii] = 2

                                                        # Change MVAR_V_Limit to V_Max
                                                        BusDataCard_DF.MVAR_V_Limit[ii] = Max_V

                                                        BusDataCard_DF.Final_V_pu[ii] = Max_V

                                                elseif (Q_i > Q_set)

                                                        if (V_Bus <= Min_V)

                                                                # Switiching PQ to PV
                                                                BusDataCard_DF.Type[ii] = 2

                                                                # Change MVAR_V_Limit to V_Min
                                                                BusDataCard_DF.MVAR_V_Limit[ii] = Min_V

                                                                BusDataCard_DF.Final_V_pu[ii] = Min_V

                                                        elseif (V_Bus >= Max_V)

                                                                # Switiching PQ to PV
                                                                BusDataCard_DF.Type[ii] = 2

                                                                # Change MVAR_V_Limit to V_Max
                                                                BusDataCard_DF.MVAR_V_Limit[ii] = Max_V

                                                                BusDataCard_DF.Final_V_pu[ii] = Max_V

                                                        else

                                                                # Keep the current Type to Original Type
                                                                BusDataCard_DF.Type[ii] = 1

                                                                # Reinitialize MVAR_V_Limit
                                                                BusDataCard_DF.MVAR_V_Limit[ii] = -9999

                                                                BusDataCard_DF.Final_V_pu[ii] = BusDataCard_DF.Final_V_pu[ii]

                                                        end

                                                end

                                        end

                                end  


                        elseif (BusDataCard_DF.Type_Original[ii] == 0) # Unregulated (load, PQ)

                                # No changes required

                        elseif (BusDataCard_DF.Type_Original[ii] == 2) # Hold voltage within VAR limits (gen, PV)

                                # Getting Q_Bus
                                Q_Bus = PQ_BusArray[ii+1,2]

                                # Getting V_set/V_i
                                V_set = BusDataCard_DF.Final_V_pu_Original[ii]
                                V_i = BusDataCard_DF.Final_V_pu[ii]
                                
                                # Getting VAR Limits
                                Max_VAR = BusDataCard_DF.Max_MVAR_V[ii]
                                Min_VAR = BusDataCard_DF.Min_MVAR_V[ii]

                                # Checking VAR Limit Violation
                                if (BusDataCard_DF.Type[ii] == 2) # Currently PV Bus is PV

                                        if (Q_Bus <= Min_VAR)

                                                # Switiching PV to PQ
                                                BusDataCard_DF.Type[ii] = 0

                                                # Change MVAR_V_Limit to MVAR_Min
                                                BusDataCard_DF.MVAR_V_Limit[ii] = Min_VAR


                                        elseif (Q_Bus >= Max_VAR)

                                                # Switiching PV to PQ
                                                BusDataCard_DF.Type[ii] = 0

                                                # Change MVAR_V_Limit to MVAR_Max
                                                BusDataCard_DF.MVAR_V_Limit[ii] = Max_VAR

                                        else

                                                # Keep the current Type to Original Type
                                                BusDataCard_DF.Type[ii] = 2

                                                # Reinitialize MVAR_V_Limit
                                                BusDataCard_DF.MVAR_V_Limit[ii] = -9999

                                                BusDataCard_DF.Final_V_pu[ii] = V_set

                                        end

                                elseif (BusDataCard_DF.Type[ii] == 0) # Currently PV Bus is PQ

                                        if (BusDataCard_DF.MVAR_V_Limit[ii] == Min_VAR)

                                                if (V_i >= V_set)

                                                        # Switiching PV to PQ
                                                        BusDataCard_DF.Type[ii] = 0

                                                        # Change MVAR_V_Limit to MVAR_Min
                                                        BusDataCard_DF.MVAR_V_Limit[ii] = Min_VAR

                                                elseif (V_i < V_set)

                                                        if (Q_Bus <= Min_VAR)

                                                                # Switiching PV to PQ
                                                                BusDataCard_DF.Type[ii] = 0

                                                                # Change MVAR_V_Limit to MVAR_Min
                                                                BusDataCard_DF.MVAR_V_Limit[ii] = Min_VAR

                                                        elseif (Q_Bus >= Max_VAR)

                                                                # Switiching PV to PQ
                                                                BusDataCard_DF.Type[ii] = 0

                                                                # Change MVAR_V_Limit to MVAR_Max
                                                                BusDataCard_DF.MVAR_V_Limit[ii] = Max_VAR

                                                        else

                                                                # Keep the current Type to Original Type
                                                                BusDataCard_DF.Type[ii] = 2

                                                                # Reinitialize MVAR_V_Limit
                                                                BusDataCard_DF.MVAR_V_Limit[ii] = -9999

                                                                BusDataCard_DF.Final_V_pu[ii] = V_set

                                                        end

                                                end

                                        elseif (BusDataCard_DF.MVAR_V_Limit[ii] == Max_VAR)

                                                if (V_i <= V_set)

                                                        # Switiching PV to PQ
                                                        BusDataCard_DF.Type[ii] = 0

                                                        # Change MVAR_V_Limit to MVAR_Min
                                                        BusDataCard_DF.MVAR_V_Limit[ii] = Max_VAR

                                                elseif (V_i > V_set)

                                                        if (Q_Bus <= Min_VAR)

                                                                # Switiching PV to PQ
                                                                BusDataCard_DF.Type[ii] = 0

                                                                # Change MVAR_V_Limit to MVAR_Min
                                                                BusDataCard_DF.MVAR_V_Limit[ii] = Min_VAR

                                                        elseif (Q_Bus >= Max_VAR)

                                                                # Switiching PV to PQ
                                                                BusDataCard_DF.Type[ii] = 0

                                                                # Change MVAR_V_Limit to MVAR_Max
                                                                BusDataCard_DF.MVAR_V_Limit[ii] = Max_VAR

                                                        else

                                                                # Keep the current Type to Original Type
                                                                BusDataCard_DF.Type[ii] = 2

                                                                # Reinitialize MVAR_V_Limit
                                                                BusDataCard_DF.MVAR_V_Limit[ii] = -9999

                                                                BusDataCard_DF.Final_V_pu[ii] = V_set

                                                        end                                                       

                                                end

                                        end

                                end                        
                                
                        end

                end 


                # Making Changes based on Bus_Modification_Counter        

                # ReOrdering BusDataCard_DF: PQ->PV->Slack
                sort!(BusDataCard_DF, [order(:Type)])

                CDF_DF_List_pu[2] = BusDataCard_DF

                # Create New Ybus
                if (Ybus_Type == 1) # Without Taps

                        Ybus = Create_Ybus_WithoutTaps(CDF_DF_List_pu)

                elseif (Ybus_Type == 2) # With Taps

                        Ybus_WithoutTaps = Create_Ybus_WithoutTaps(CDF_DF_List_pu)

                        Ybus = Create_Ybus_WithTaps(Ybus_WithoutTaps,CDF_DF_List_pu)

                end

                # Get SolutionVector_V
                SolutionVector_V = BusDataCard_DF.Final_V_pu
                SolutionVector_V = vcat(SolutionVector_V[end], SolutionVector_V[1:(end-1)])

                # Get SolutionVector_Delta
                SolutionVector_Delta = BusDataCard_DF.Final_A_deg
                SolutionVector_Delta = vcat(SolutionVector_Delta[end], SolutionVector_Delta[1:(end-1)])

                #  Compute P-Q at Buses
                PQ_BusArray = Compute_PQ_BusArray(Ybus, SolutionVector_V, SolutionVector_Delta)

                #  Compute PQ Mistmatch
                PQ_MismatchVector = Compute_PQ_MismatchVector(CDF_DF_List_pu, PQ_BusArray, SolutionVector_V, NR_Type)

        else

                #  Compute P-Q at Buses
                PQ_BusArray = Compute_PQ_BusArray(Ybus, SolutionVector_V, SolutionVector_Delta)

                #  Compute PQ Mistmatch
                PQ_MismatchVector = Compute_PQ_MismatchVector(CDF_DF_List_pu, PQ_BusArray, SolutionVector_V, NR_Type)

        end

        return CDF_DF_List_pu, Ybus, PQ_BusArray, PQ_MismatchVector, SolutionVector_V, SolutionVector_Delta

end

"""
    Compute_LineFlows(CDF_DF_List_pu, Ybus, SolutionVector_V, SolutionVector_Delta)

Computes P and Q at each bus of a power system network.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
- 'Ybus': A complex array of Ybus elements ordered according to bus
type: Slack->PQ->PV.
- 'SolutionVector_NR': Voltage and Angle at each bus ordered according to bus
type: PQ->PV.
'''
'''
# Output
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF]. Changed with changes in Line flows
- 'LineFlow_Array': An array (N_Lines*4) for P and Q line flows ordered
according to Branch Data Card.
'''
"""
function Compute_LineFlows(CDF_DF_List_pu, Ybus, SolutionVector_NR)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        BranchDataCard_DF = CDF_DF_List_pu[3]

        # Creating SolutionVector_V and SolutionVector_Delta
        SolutionVector_V, SolutionVector_Delta = Create_SolutionVector_VDelta_NR(CDF_DF_List_pu, SolutionVector_NR)

        # Number of Buses and Lines
        N_Bus = nrow(BusDataCard_DF)

        N_Lines = nrow(BranchDataCard_DF)

        # Updating Final_V_pu and Final_A_deg
        for ii in 1:N_Bus

                if (ii == N_Bus)

                        BusDataCard_DF.Final_V_pu[ii] = SolutionVector_V[1,1]

                        BusDataCard_DF.Final_A_deg[ii] = SolutionVector_Delta[1,1]

                else

                        BusDataCard_DF.Final_V_pu[ii] = SolutionVector_V[ii+1,1]

                        BusDataCard_DF.Final_A_deg[ii] = SolutionVector_Delta[ii+1,1]

                end

        end        

        # Initializing LineFlow_Array
        LineFlow_Array = Array{Float64}(undef, N_Lines,4)

        # Initializing Bus_i_Index/Bus_j_Index - for getting reference outside For Loop
        Bus_i_Index = 0
        Bus_j_Index = 0

        # Computing Elements of LineFlow_Array
        for ii in 1:N_Lines # Through Rows of BranchDataCard_DF/LineFlow_Array

                # Getting Bus Numbers for the current Branch
                Bus_i_Num = BranchDataCard_DF.Tap_Bus_Num[ii]

                Bus_j_Num = BranchDataCard_DF.Z_Bus_Num[ii]

                # Getting Bus_i_Index
                for jj in 1:N_Bus # Through Rows of BusDataCard_DF

                        if (Bus_i_Num == BusDataCard_DF.Bus_Num[jj])

                                if (jj == N_Bus) # Accounting for Slack Bus at the end of BusDataCard_DF

                                        Bus_i_Index = 1

                                else

                                        Bus_i_Index = jj+1

                                end

                        end

                end

                # Getting Bus_j_Index
                for jj in 1:N_Bus # Through Rows of BusDataCard_DF

                        if (Bus_j_Num == BusDataCard_DF.Bus_Num[jj])

                                if (jj == N_Bus) # Accounting for Slack Bus at the end of BusDataCard_DF

                                        Bus_j_Index = 1

                                else

                                        Bus_j_Index = jj+1

                                end

                        end

                end

                # Computing P_ij
                P_ij = -(SolutionVector_V[Bus_i_Index]^(2)*real(Ybus[Bus_i_Index, Bus_j_Index])) +(SolutionVector_V[Bus_i_Index]*SolutionVector_V[Bus_j_Index]*abs(Ybus[Bus_i_Index, Bus_j_Index]))*cos(angle(Ybus[Bus_i_Index, Bus_j_Index])+deg2rad(SolutionVector_Delta[Bus_j_Index])-deg2rad(SolutionVector_Delta[Bus_i_Index]))

                # Computing Q_ij
                Q_ij = -((SolutionVector_V[Bus_i_Index]^(2)*((BranchDataCard_DF.B_pu[ii]/2) - imag(Ybus[Bus_i_Index, Bus_j_Index])))) - (SolutionVector_V[Bus_i_Index]*SolutionVector_V[Bus_j_Index]*abs(Ybus[Bus_i_Index, Bus_j_Index]))*sin(angle(Ybus[Bus_i_Index, Bus_j_Index])+deg2rad(SolutionVector_Delta[Bus_j_Index])-deg2rad(SolutionVector_Delta[Bus_i_Index]))

                # Computing P_ji
                P_ji = -(SolutionVector_V[Bus_j_Index]^(2)*real(Ybus[Bus_j_Index, Bus_i_Index])) +(SolutionVector_V[Bus_j_Index]*SolutionVector_V[Bus_i_Index]*abs(Ybus[Bus_j_Index, Bus_i_Index]))*cos(angle(Ybus[Bus_j_Index, Bus_i_Index])+deg2rad(SolutionVector_Delta[Bus_i_Index])-deg2rad(SolutionVector_Delta[Bus_j_Index]))

                # Computing Q_ij
                Q_ji = -((SolutionVector_V[Bus_j_Index]^(2)*((BranchDataCard_DF.B_pu[ii]/2) - imag(Ybus[Bus_j_Index, Bus_i_Index])))) - (SolutionVector_V[Bus_j_Index]*SolutionVector_V[Bus_i_Index]*abs(Ybus[Bus_j_Index, Bus_i_Index]))*sin(angle(Ybus[Bus_j_Index, Bus_i_Index])+deg2rad(SolutionVector_Delta[Bus_i_Index])-deg2rad(SolutionVector_Delta[Bus_j_Index]))

                # Addressing Machine Precision Problem and Filling up LineFlow_Array/Updating BranchDataCard_DF
                if (abs(P_ij) < 1e-12)
                        LineFlow_Array[ii,1] = 0
                        BranchDataCard_DF.Line_Flow_P_ij[ii] = 0
                else
                        LineFlow_Array[ii,1] = P_ij
                        BranchDataCard_DF.Line_Flow_P_ij[ii] = P_ij
                end

                if (abs(P_ji) < 1e-12)
                        LineFlow_Array[ii,2] = 0
                        BranchDataCard_DF.Line_Flow_P_ji[ii] = 0
                else
                        LineFlow_Array[ii,2] = P_ji
                        BranchDataCard_DF.Line_Flow_P_ji[ii] = P_ji
                end

                if (abs(Q_ij) < 1e-12)
                        LineFlow_Array[ii,3] = 0
                        BranchDataCard_DF.Line_Flow_Q_ij[ii] = 0
                else
                        LineFlow_Array[ii,3] = Q_ij
                        BranchDataCard_DF.Line_Flow_Q_ij[ii] = Q_ij
                end

                if (abs(Q_ji) < 1e-12)
                        LineFlow_Array[ii,4] = 0
                        BranchDataCard_DF.Line_Flow_Q_ji[ii] = 0
                else
                        LineFlow_Array[ii,4] = Q_ji
                        BranchDataCard_DF.Line_Flow_Q_ji[ii] = Q_ji
                end

        end


        return CDF_DF_List_pu, LineFlow_Array

end

"""
    Compute_ToleranceSatisfaction(Tolerance, Corrections_Vector)

Computes if the correction vector satifies the tolerance.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
- 'Ybus': A complex array of Ybus elements ordered according to bus
type: Slack->PQ->PV.
- 'NR_Type': 1 -> Full Newton-Raphson, 2-> Decoupled Newton-Raphson,
3 -> Fast Decoupled Newton-Raphson
- 'Tolerance': Tolerance level for stopping criterion of Newton-Raphson Method.
- 'SolutionVector_NR': Voltage and Angle at each bus ordered according to bus
type: PQ->PV.
'''
'''
# Output
- 'Tolerance_Satisfaction': A boolean value indicating tolerance satifaction
(True) and disatisfaction (False).
'''
"""
function Compute_ToleranceSatisfaction(CDF_DF_List_pu, Ybus, NR_Type, Tolerance, SolutionVector_NR)
        
    # Creating SolutionVector_V and SolutionVector_Delta
    SolutionVector_V, SolutionVector_Delta = Create_SolutionVector_VDelta_NR(CDF_DF_List_pu, SolutionVector_NR)

    # @show SolutionVector_V
    # @show SolutionVector_Delta

    #  Compute P-Q at Buses
    PQ_BusArray = Compute_PQ_BusArray(Ybus, SolutionVector_V, SolutionVector_Delta)

    #  Compute PQ Mistmatch
    PQ_MismatchVector = Compute_PQ_MismatchVector(CDF_DF_List_pu, PQ_BusArray, SolutionVector_V, NR_Type)
    
    # Initializing Tolerance Counter
    ToleranceCounter = 0

    # Computing Tolerance Counter
    for ii in 1: length(PQ_MismatchVector)

            if (abs(PQ_MismatchVector[ii]) > Tolerance)

                    # Incrementing the Tolerance Counter
                    ToleranceCounter = ToleranceCounter + 1

            else

                    # Do nothing

            end

    end

    # Computing Tolerance Satisfaction
    if (ToleranceCounter > 0)

          Tolerance_Satisfaction = false

    else

            Tolerance_Satisfaction = true

    end

    return Tolerance_Satisfaction

end

"""
    Compute_Corrected_CorrectionVector(CDF_DF_List_pu, Correction_Vector_NR,
    SolutionVector_NR, NR_Type)

Computes if the correction vector satifies the tolerance.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
- 'Correction_Vector_NR': The correction vector computed during Newton-Raphson
Method.
- 'SolutionVector_NR': Voltage and Angle at each bus ordered according to bus
type: PQ->PV.
- 'NR_Type': 1 -> Full Newton-Raphson, 2-> Decoupled Newton-Raphson,
3 -> Fast Decoupled Newton-Raphson
'''
'''
# Output
- 'Correction_Vector_NR_Corrected': The correction vector computed during Newton-Raphson
Method corrected for NR_Type.
'''
"""
function Compute_Corrected_CorrectionVector(CDF_DF_List_pu, Correction_Vector_NR, SolutionVector_NR, NR_Type)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type == 0) || (row.Type == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

        # Initializing Correction_Vector_NR_Corrected
        Correction_Vector_NR_Corrected = Correction_Vector_NR

        # If Else Loop: For NR_Type
        if ((NR_Type == 1) || (NR_Type == 2)) # Full NR Decoupled NR

                for ii in 1 : (N_PQ_Bus + N_PV_Bus)

                        Correction_Vector_NR_Corrected[ii] = rad2deg(Correction_Vector_NR_Corrected[ii])

                end
 
 
                for ii in (N_PQ_Bus + N_PV_Bus + 1) : length(SolutionVector_NR)

                        Correction_Vector_NR_Corrected[ii] = Correction_Vector_NR_Corrected[ii]*SolutionVector_NR[ii]

                end

        elseif (NR_Type == 3) # Fast Decoupled NR

                Correction_Vector_NR_Corrected = Correction_Vector_NR

                for ii in 1 : (N_PQ_Bus + N_PV_Bus)

                        Correction_Vector_NR_Corrected[ii] = rad2deg(Correction_Vector_NR_Corrected[ii])

                end

        end

        return Correction_Vector_NR_Corrected

end


"""
    Compute_Corrected_Matrix(Matrix)

Computes corrected matrix for elements close to zero for machine precision.

'''
# Arguments
- 'Matrix': 
'''
'''
# Output
- 'Matrix': Elements close to zero by machine precision are converted to zero.
'''
"""
function Compute_Corrected_Matrix(Matrix)

        # Define Machine precision
        Precision = 1e-12

        # Getting number of Rows and Columns of Matrix
        Row_Num = size(Matrix)[1]
        Col_Num = size(Matrix)[2]

        # For Loop: For correcting each element in Matrix
        for ii in 1:Row_Num # Through Rows

                for jj in 1:Col_Num # Through Columns 

                        # If Loop: For correcting current Matrix element
                        if (abs(Matrix[ii,jj]) < Precision)

                                Matrix[ii,jj] = 0.0 

                        end

                end

        end

        return Matrix

end



function plotBuswiseDifferences(CDF_DF_List_pu::Vector{DataFrame},
        results::DataFrame;
        savePlots::Bool=true,
        processedDataFolder::String="processedData/",
        fileExtension::String=".png")

        busData = CDF_DF_List_pu[2]
        systemName = extractSystemName(CDF_DF_List_pu)
        N = size(results, 1)
        # Compute the relative differences in % for V
        rel_diff_V = 100 * (results.V - busData.Final_V_pu_Original) ./ busData.Final_V_pu_Original
        # Compute the absolute differences for δ
        abs_diff_δ = 180/π * results.δ - busData.Final_A_deg_Original
        
        # Bar plots
        p1 = bar(1:N, rel_diff_V, color=:green,
        alpha=0.7,
        linewidth=1.5,
        ylabel="% Difference",
        xlabel="Bus Number",
        title="Relative Difference in Official and Computed Voltages\n"*
        "for the $(systemName) Bus System.",
        titlefontsize=12,
        xticks=1:N,
        label=L"$ \frac{ΔV}{V} * 100 \% $ where $ΔV = V_{computed} - V_{CDF}$",
        legendfontsize=8);

        p2 = bar(1:N, abs_diff_δ, color=:orange,
        alpha=0.7,
        linewidth=1.5,
        ylabel="Absolute Difference in Degrees",
        xlabel="Bus Number",
        title="Absolute Difference in Official and Computed Angles\n"*
        "for the $(systemName) Bus System.",
        titlefontsize=12,
        xticks=1:N,
        label=L"$\Delta \delta \, [\degree]$ where $\Delta \delta = \delta_{computed} - \delta_{CDF}$",
        legendfontsize=8);
        # Try specifying a larger size
        p = plot(p1, p2, layout=(2,1), size=(800, 800))
        # bar(p1)

        if savePlots
                filename = processedDataFolder*"results_"*systemName*"_sparse"*fileExtension;
                savefig(p, filename)   
        end
end

