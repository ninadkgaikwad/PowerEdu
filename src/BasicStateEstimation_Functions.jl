# BasicStateEstimation_Functions.jl

"""
    CDF_AddMeasurements_SE(CDF_DF_List, Bus_Measurement_Array,
    Branch_Measurement_Array)

Adds Bus and Branch Measurements to CDF_DF_List for State Estimation.

'''
# Arguments
- 'CDF_DF_List': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
- 'Bus_Measurement_Array': Bus measurement information in Array format -
[['Bus Number', 'V Measurement', 'P Measurement', 'Q Measurement',
Error Variance V, Error Variance P, Error Variance Q],...] if
measurement not present -9999 should be added.
- 'Branch_Measurement_Array': Branch measurement information in Array format -
[['Bus i Number', 'Bus j Number', 'P + Measurement', 'P - Measurement',
'Q + Measurement', 'Q - Measurement', Error Variance P Pos,
Error Variance P Neg, Error Variance Q Pos, Error Variance Q Neg],...]
if measurement not present -9999 should be added, '+' means i->j and '-' means
j->i.
'''
'''
# Output
- 'CDF_DF_List': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF] with State Estimation Measurements added.
'''
"""
function CDF_AddMeasurements_SE(CDF_DF_List, Bus_Measurement_Array, Branch_Measurement_Array)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List[2]
        BranchDataCard_DF = CDF_DF_List[3]

        # Getting Row Numbers of the DataFrames
        N_Bus = nrow(BusDataCard_DF)
        N_Branch = nrow(BranchDataCard_DF)

        # Getting length of Bus/Branch Measurement Arrays        
        N_Bus_Measurement = size(Bus_Measurement_Array)[1]
        N_Branch_Measurement = size(Branch_Measurement_Array)[1]

        # For Loop: For adding Bus Measurements
        for ii in 1:N_Bus_Measurement

                # Getting Current Bus Measurements
                Current_Bus_Measurement = Bus_Measurement_Array[ii,:]

                # Getting Components of Current_Bus_Measurement
                Current_BusNumber = Current_Bus_Measurement[1]
                Current_V_Measured = Current_Bus_Measurement[2]
                Current_P_Measured = Current_Bus_Measurement[3]
                Current_Q_Measured = Current_Bus_Measurement[4]
                Current_V_EV = Current_Bus_Measurement[5]
                Current_P_EV = Current_Bus_Measurement[6]
                Current_Q_EV = Current_Bus_Measurement[7]

                # For Loop: For going through rows of BusDataCard_DF
                for jj in 1:N_Bus

                        if (BusDataCard_DF.Bus_Num[jj] == Current_BusNumber)

                                if (Current_V_Measured != -9999)

                                        BusDataCard_DF.V_Measured[jj] = Current_V_Measured

                                end

                                if (Current_P_Measured != -9999)

                                        BusDataCard_DF.P_Measured[jj] = Current_P_Measured

                                end

                                if (Current_Q_Measured != -9999)

                                        BusDataCard_DF.Q_Measured[jj] = Current_Q_Measured

                                end

                                if (Current_V_EV != -9999)

                                        BusDataCard_DF.ErrorVariance_V[jj] = Current_V_EV

                                end

                                if (Current_P_EV != -9999)

                                        BusDataCard_DF.ErrorVariance_P[jj] = Current_P_EV

                                end

                                if (Current_Q_EV != -9999)

                                        BusDataCard_DF.ErrorVariance_Q[jj] = Current_Q_EV

                                end

                        end

                end

        end

        # For Loop: For adding Branch Measurements
        for ii in 1:N_Branch_Measurement

                # Getting Current Bus Measurements
                Current_Branch_Measurement = Branch_Measurement_Array[ii,:]

                # Getting Components of Current_Bus_Measurement
                Current_BusNumber_i = Current_Branch_Measurement[1]
                Current_BusNumber_j = Current_Branch_Measurement[2]
                Current_P_Pos_Measured = Current_Branch_Measurement[3]
                Current_P_Neg_Measured = Current_Branch_Measurement[4]
                Current_Q_Pos_Measured = Current_Branch_Measurement[5]
                Current_Q_Neg_Measured = Current_Branch_Measurement[6]
                Current_P_Pos_EV = Current_Branch_Measurement[7]
                Current_P_Neg_EV = Current_Branch_Measurement[8]
                Current_Q_Pos_EV = Current_Branch_Measurement[9]
                Current_Q_Neg_EV = Current_Branch_Measurement[10]

                # For Loop: For going through rows of BranchDataCard_DF
                for jj in 1:N_Branch

                        if ((BranchDataCard_DF.Tap_Bus_Num[jj] == Current_BusNumber_i) && (BranchDataCard_DF.Z_Bus_Num[jj] == Current_BusNumber_j))

                                if (Current_P_Pos_Measured != -9999)

                                        BranchDataCard_DF.Line_Flow_Pos_P[jj] = Current_P_Pos_Measured

                                end

                                if (Current_P_Neg_Measured != -9999)

                                        BranchDataCard_DF.Line_Flow_Neg_P[jj] = Current_P_Neg_Measured

                                end

                                if (Current_Q_Pos_Measured != -9999)

                                        BranchDataCard_DF.Line_Flow_Pos_Q[jj] = Current_Q_Pos_Measured

                                end

                                if (Current_Q_Neg_Measured != -9999)

                                        BranchDataCard_DF.Line_Flow_Neg_Q[jj] = Current_Q_Neg_Measured

                                end

                                if (Current_P_Pos_EV != -9999)

                                        BranchDataCard_DF.ErrorVariance_Pos_P[jj] = Current_P_Pos_EV

                                end

                                if (Current_P_Neg_EV != -9999)

                                        BranchDataCard_DF.ErrorVariance_Neg_P[jj] = Current_P_Neg_EV

                                end

                                if (Current_Q_Pos_EV != -9999)

                                        BranchDataCard_DF.ErrorVariance_Pos_Q[jj] = Current_Q_Pos_EV

                                end

                                if (Current_Q_Neg_EV != -9999)

                                        BranchDataCard_DF.ErrorVariance_Neg_Q[jj] = Current_Q_Neg_EV

                                end

                        elseif ((BranchDataCard_DF.Tap_Bus_Num[jj] == Current_BusNumber_j) && (BranchDataCard_DF.Z_Bus_Num[jj] == Current_BusNumber_i))

                                if (Current_P_Pos_Measured != -9999)

                                        BranchDataCard_DF.Line_Flow_Neg_P[jj] = Current_P_Pos_Measured

                                end

                                if (Current_P_Neg_Measured != -9999)

                                        BranchDataCard_DF.Line_Flow_Pos_P[jj] = Current_P_Neg_Measured

                                end

                                if (Current_Q_Pos_Measured != -9999)

                                        BranchDataCard_DF.Line_Flow_Neg_Q[jj] = Current_Q_Pos_Measured

                                end

                                if (Current_Q_Pos_Measured != -9999)

                                        BranchDataCard_DF.Line_Flow_Neg_Q[jj] = Current_Q_Neg_Measured

                                end

                                if (Current_P_Pos_EV != -9999)

                                        BranchDataCard_DF.ErrorVariance_Neg_P[jj] = Current_P_Pos_EV

                                end

                                if (Current_P_Neg_EV != -9999)

                                        BranchDataCard_DF.ErrorVariance_Pos_P[jj] = Current_P_Neg_EV

                                end

                                if (Current_Q_Pos_EV != -9999)

                                        BranchDataCard_DF.ErrorVariance_Neg_Q[jj] = Current_Q_Pos_EV

                                end

                                if (Current_Q_Neg_EV != -9999)

                                        BranchDataCard_DF.ErrorVariance_Pos_Q[jj] = Current_Q_Neg_EV

                                end

                        end

                end

        end

        # Putting New BusDataCard_DF and BranchDataCard_DF in CDF_DF_List
        CDF_DF_List[2] = BusDataCard_DF
        CDF_DF_List[3] = BranchDataCard_DF

        return CDF_DF_List

end

"""
    CDF_pu_Converter_SE(CDF_DF_List)

Converts CDF_DF_List objects to per unit (pu) for State Estimation.

'''
# Arguments
- 'CDF_DF_List': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
'''
'''
# Output
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file coverted to pu: [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF] with SE Measurements.
'''
"""
function CDF_pu_Converter_SE(CDF_DF_List)

        # Getting required data from CDF_DF_List
        TitleCard_DF = CDF_DF_List[1]
        BusDataCard_DF = CDF_DF_List[2]
        BranchDataCard_DF = CDF_DF_List[3]
        LossZonesCard_DF = CDF_DF_List[4]
        InterchangeDataCard_DF = CDF_DF_List[5]
        TieLinesDataCard_DF = CDF_DF_List[6]

        # Getting Base MVA
        Base_MVA = TitleCard_DF.MVA_Base[1]

        # Converting Actual Value Columns in BusDataCard_DF to pu
        for ii in 1:length(BusDataCard_DF.Bus_Num)

                # Converting Load MW/MVAR and Gen MW/MVAR to pu
                BusDataCard_DF.Load_MW[ii] = BusDataCard_DF.Load_MW[ii]/Base_MVA

                BusDataCard_DF.Load_MVAR[ii] = BusDataCard_DF.Load_MVAR[ii]/Base_MVA

                BusDataCard_DF.Gen_MW[ii] = BusDataCard_DF.Gen_MW[ii]/Base_MVA

                BusDataCard_DF.Gen_MVAR[ii] = BusDataCard_DF.Gen_MVAR[ii]/Base_MVA

                # Converting Max_MVAR_V and Min_MVAR_V to pu
                if (BusDataCard_DF.Type[ii] == 0) # Unregulated PQ

                        # Do Nothing

                elseif (BusDataCard_DF.Type[ii] == 1) # Voltage Limit PQ

                        BusDataCard_DF.Max_MVAR_V[ii] = BusDataCard_DF.Max_MVAR_V[ii]/BusDataCard_DF.Base_KV[ii]

                        BusDataCard_DF.Min_MVAR_V[ii] = BusDataCard_DF.Min_MVAR_V[ii]/BusDataCard_DF.Base_KV[ii]

                elseif (BusDataCard_DF.Type[ii] == 2) # VAR Limit PV

                        BusDataCard_DF.Max_MVAR_V[ii] = BusDataCard_DF.Max_MVAR_V[ii]/Base_MVA

                        BusDataCard_DF.Min_MVAR_V[ii] = BusDataCard_DF.Min_MVAR_V[ii]/Base_MVA

                elseif (BusDataCard_DF.Type[ii] == 3) # Slack

                        # Do Nothing

                end

                # Converting SE Measurements to Pu
                if (BusDataCard_DF.V_Measured[ii] != -9999)

                        BusDataCard_DF.V_Measured[ii] = BusDataCard_DF.V_Measured[ii]/BusDataCard_DF.Base_KV[ii]

                end

                if (BusDataCard_DF.P_Measured[ii] != -9999)

                        BusDataCard_DF.P_Measured[ii] = BusDataCard_DF.P_Measured[ii]/Base_MVA

                end

                if (BusDataCard_DF.Q_Measured[ii] != -9999)

                        BusDataCard_DF.Q_Measured[ii] = BusDataCard_DF.Q_Measured[ii]/Base_MVA

                end


        end

        # Converting Actual Value Columns in BranchDataCard_DF to pu
        for ii in 1:length(BranchDataCard_DF.Tap_Bus_Num)

                # Converting Max_MVAR_MW_V and Min_MVAR_MW_V to pu
                if (BranchDataCard_DF.Type[ii] == 0) # Transmission Line

                        # Do Nothing

                elseif (BranchDataCard_DF.Type[ii] == 1) # Fixed Tap

                        # Do Nothing

                elseif (BranchDataCard_DF.Type[ii] == 2) # Variable Tap V Control

                        # Getting Base KV from BusDataCard_DF
                        Base_KV = filter(row -> row.Bus_Num == BranchDataCard_DF.Tap_Bus_Num[ii], BusDataCard_DF).Base_KV

                        # Converting to pu
                        BranchDataCard_DF.Max_MVAR_MW_V[ii] = BranchDataCard_DF.Max_MVAR_MW_V[ii]/Base_KV

                        BranchDataCard_DF.Min_MVAR_MW_V[ii] = BranchDataCard_DF.Min_MVAR_MW_V[ii]/Base_KV

                elseif (BranchDataCard_DF.Type[ii] == 3) # Variable Tap MVAR Control

                        BranchDataCard_DF.Max_MVAR_MW_V[ii] = BranchDataCard_DF.Max_MVAR_MW_V[ii]/Base_MVA

                        BranchDataCard_DF.Min_MVAR_MW_V[ii] = BranchDataCard_DF.Min_MVAR_MW_V[ii]/Base_MVA

                elseif (BranchDataCard_DF.Type[ii] == 4) # Variable Phase Angle MW Control

                        BranchDataCard_DF.Max_MVAR_MW_V[ii] = BranchDataCard_DF.Max_MVAR_MW_V[ii]/Base_MVA

                        BranchDataCard_DF.Min_MVAR_MW_V[ii] = BranchDataCard_DF.Min_MVAR_MW_V[ii]/Base_MVA

                end

                # Converting SE Measurements to Pu
                if (BranchDataCard_DF.Line_Flow_Pos_P[ii] != -9999)

                        BranchDataCard_DF.Line_Flow_Pos_P[ii] = BranchDataCard_DF.Line_Flow_Pos_P[ii]/Base_MVA

                end

                if (BranchDataCard_DF.Line_Flow_Neg_P[ii] != -9999)

                        BranchDataCard_DF.Line_Flow_Neg_P[ii] = BranchDataCard_DF.Line_Flow_Neg_P[ii]/Base_MVA

                end

                if (BranchDataCard_DF.Line_Flow_Pos_Q[ii] != -9999)

                        BranchDataCard_DF.Line_Flow_Pos_Q[ii] = BranchDataCard_DF.Line_Flow_Pos_Q[ii]/Base_MVA

                end

                if (BranchDataCard_DF.Line_Flow_Neg_Q[ii] != -9999)

                        BranchDataCard_DF.Line_Flow_Neg_Q[ii] = BranchDataCard_DF.Line_Flow_Neg_Q[ii]/Base_MVA

                end

        end

        # Converting Actual Value Columns in InterchangeDataCard_DF to pu
        for ii in 1:length(InterchangeDataCard_DF.Area_Num)

                # Converting Area_Inter_export and Area_Inter_tol to pu
                InterchangeDataCard_DF.Area_Inter_export[ii] = InterchangeDataCard_DF.Area_Inter_export[ii]/Base_MVA

                InterchangeDataCard_DF.Area_Inter_tol[ii] = InterchangeDataCard_DF.Area_Inter_tol[ii]/Base_MVA

        end

        return CDF_DF_List_pu = [TitleCard_DF, BusDataCard_DF, BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF, TieLinesDataCard_DF]

end


"""
    Create_Initial_SolutionVector_SE(CDF_DF_List_pu)

Creates Initial Solution Vector for the Power System State Estimation.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
'''
'''
# Output
- 'Initial_SolutionVector_SE': Power Flow Solution Voltage and Angle at each
bus ordered according to bus type: Slack->PQ->PV.
'''
"""
function Create_Initial_SolutionVector_SE(CDF_DF_List_pu)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type_Original == 0) || (row.Type_Original == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type_Original == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type_Original == 3), BusDataCard_DF))

        # Creating Initial_SolutionVector_SE Length
        Initial_SolutionVector_SE_Len = (2*N_PQ_Bus)+(2*N_PV_Bus)+1

        # Initializing Initial_SolutionVector_SE
        Initial_SolutionVector_SE = Array{Float64}(undef, Initial_SolutionVector_SE_Len ,1)

        # Creating Initial_SolutionVector_SE Delta part for both PQ and PV Buses
        Initial_SolutionVector_SE[1:(N_PQ_Bus+N_PV_Bus),1] = BusDataCard_DF.Final_A_deg[1:end-1]

        # Creating Initial_SolutionVector_SE V part for Slack-PQ-PV Buses
        Initial_SolutionVector_SE[(N_PQ_Bus+N_PV_Bus+1+1):end,1] = BusDataCard_DF.Final_V_pu[1:end-1]      
        
        Initial_SolutionVector_SE[(N_PQ_Bus+N_PV_Bus+1),1] = BusDataCard_DF.Final_V_pu[end]   

        return Initial_SolutionVector_SE

end

"""
Create_SolutionVector_VDelta_SE(SolutionVector_NR)

Creates separate V and Delta vectors from the solution vector computed by the State Estimation.

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
function Create_SolutionVector_VDelta_SE(CDF_DF_List_pu, SolutionVector_SE)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type_Original == 0) || (row.Type_Original == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type_Original == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type_Original == 3), BusDataCard_DF))

        # Getting Length of SolutionVector_SES
        Len_VDelta_Vector = N_Bus        

        # Initializing SolutionVector_V, SolutionVector_Delta
        SolutionVector_V = Array{Float64}(undef, Len_VDelta_Vector,1)
        SolutionVector_Delta = Array{Float64}(undef, Len_VDelta_Vector,1)

        # Creating SolutionVector_V, SolutionVector_Delta
        SolutionVector_Delta[2:end,1] =  SolutionVector_SE[1:Len_VDelta_Vector-1,1]

        # Adding Slack bus angle [Reference Angle for State Estimation]
        SolutionVector_Delta[1,1] =  0
        
        SolutionVector_V[1:end,1] =  SolutionVector_SE[Len_VDelta_Vector:end,1]        
        
        return SolutionVector_V, SolutionVector_Delta

end


""" 
Create_BusBranch_Measurement_Arrays_SE(CDF_DF_List_pu, LineFlow_Array, Measurement_Error_Variance)

Creates Initial Solution Vector for the Power System State Estimation.

'''     
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
- 'LineFlow_Array': An array (N_Lines*2) for P and Q line flows ordered
according to Branch Data Card.
- 'PQ_BusArray': An array (N*4) for P and Q vectors ordered according to bus
type: Slack->PQ->PV.
- 'Measurement_Error_Variance': Vector of measurement error variance in order V -> P -> Q.
'''
'''
# Output
- 'Bus_Measurement_Array': Bus measurement information in Array format -
[['Bus Number', 'V Measurement', 'P Measurement', 'Q Measurement',
Error Variance V, Error Variance P, Error Variance Q],...] if
measurement not present -9999 should be added.
- 'Branch_Measurement_Array': Branch measurement information in Array format -
[['Bus i Number', 'Bus j Number', 'P + Measurement', 'P - Measurement',
'Q + Measurement', 'Q - Measurement', Error Variance P Pos,
Error Variance P Neg, Error Variance Q Pos, Error Variance Q Neg],...]
if measurement not present -9999 should be added, '+' means i->j and '-' means
j->i.
'''
"""
function Create_BusBranch_Measurement_Arrays_SE(CDF_DF_List_pu, LineFlow_Array, Bus_PQ_Array, Measurement_Error_Variance)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]
        BranchDataCard_DF = CDF_DF_List_pu[3]

        # Number of Buses and Lines
        N_Bus = nrow(BusDataCard_DF)
        N_Lines = nrow(BranchDataCard_DF)

        # Getting Error Variances
        V_Normal_Dis = Normal(0, Measurement_Error_Variance[1])
        P_Normal_Dis = Normal(0, Measurement_Error_Variance[2])
        Q_Normal_Dis = Normal(0, Measurement_Error_Variance[3])

        # Initializing LineFlow_Array
        Bus_Measurement_Array = Array{Float64}(undef, N_Bus,7)

        Branch_Measurement_Array = Array{Float64}(undef, N_Lines,10)

        # For Loop: For each bus in BusDataCard_DF
        for ii in 1:N_Bus

                if (ii == N_Bus)

                        # Inserting Measurments in Bus_Measurement_Array
                        Bus_Measurement_Array[ii,1] = BusDataCard_DF.Bus_Num[ii]
                        Bus_Measurement_Array[ii,2] = BusDataCard_DF.Final_V_pu[ii] + rand(V_Normal_Dis, 1)[1]
                        Bus_Measurement_Array[ii,3] = Bus_PQ_Array[1,1] + rand(P_Normal_Dis, 1)[1]
                        Bus_Measurement_Array[ii,4] = Bus_PQ_Array[1,2] + rand(Q_Normal_Dis, 1)[1]
                        Bus_Measurement_Array[ii,5] = Measurement_Error_Variance[1]
                        Bus_Measurement_Array[ii,6] = Measurement_Error_Variance[2]
                        Bus_Measurement_Array[ii,7] = Measurement_Error_Variance[3]

                else

                        # Inserting Measurments in Bus_Measurement_Array
                        Bus_Measurement_Array[ii,1] = BusDataCard_DF.Bus_Num[ii]
                        Bus_Measurement_Array[ii,2] = BusDataCard_DF.Final_V_pu[ii] + rand(V_Normal_Dis, 1)[1]
                        Bus_Measurement_Array[ii,3] = Bus_PQ_Array[ii+1,1] + rand(P_Normal_Dis, 1)[1]
                        Bus_Measurement_Array[ii,4] = Bus_PQ_Array[ii+1,2] + rand(Q_Normal_Dis, 1)[1]
                        Bus_Measurement_Array[ii,5] = Measurement_Error_Variance[1]
                        Bus_Measurement_Array[ii,6] = Measurement_Error_Variance[2]
                        Bus_Measurement_Array[ii,7] = Measurement_Error_Variance[3]

                end
                

        end

        # For Loop: For each line in BranchDataCard_DF
        for ii in 1:N_Lines

                # Inserting Measurments in Branch_Measurement_Array
                Branch_Measurement_Array[ii,1] = BranchDataCard_DF.Tap_Bus_Num[ii]
                Branch_Measurement_Array[ii,2] = BranchDataCard_DF.Z_Bus_Num[ii]
                Branch_Measurement_Array[ii,3] = LineFlow_Array[ii,1] + rand(P_Normal_Dis, 1)[1]
                Branch_Measurement_Array[ii,4] = LineFlow_Array[ii,2] + rand(P_Normal_Dis, 1)[1]
                Branch_Measurement_Array[ii,5] = LineFlow_Array[ii,3] + rand(Q_Normal_Dis, 1)[1]
                Branch_Measurement_Array[ii,6] = LineFlow_Array[ii,4] + rand(P_Normal_Dis, 1)[1]
                Branch_Measurement_Array[ii,7] = Measurement_Error_Variance[2]
                Branch_Measurement_Array[ii,8] = Measurement_Error_Variance[2]
                Branch_Measurement_Array[ii,9] = Measurement_Error_Variance[3]     
                Branch_Measurement_Array[ii,10] = Measurement_Error_Variance[3]                  

        end

        return Bus_Measurement_Array, Branch_Measurement_Array 

end

""" 
Create_Bad_BusBranch_Measurement_Arrays_SE(Bus_Measurement_Array, Branch_Measurement_Array, Bad_Bus_Measurement_Input, Bad_Branch_Measurement_Input)

Creates Initial Solution Vector for the Power System State Estimation.

'''     
# Arguments
- 'Bus_Measurement_Array': Bus measurement information in Array format -
[['Bus Number', 'V Measurement', 'P Measurement', 'Q Measurement',
Error Variance V, Error Variance P, Error Variance Q],...] if
measurement not present -9999 should be added.
- 'Branch_Measurement_Array': Branch measurement information in Array format -
[['Bus i Number', 'Bus j Number', 'P + Measurement', 'P - Measurement',
'Q + Measurement', 'Q - Measurement', Error Variance P Pos,
Error Variance P Neg, Error Variance Q Pos, Error Variance Q Neg],...]
if measurement not present -9999 should be added, '+' means i->j and '-' means
j->i.
- 'Bad_Bus_Measurement_Input': Array of format -  [[Bus_Num, V_scale, P_scale, Q_scale], []] when bad data present, 
otherwise 'nothing'; the measurements are scaled through scale variables to make bad data.
- 'Bad_Branch_Measurement_Input': Array of format - [[Bus_Num_i, Bus_Num_j, P_+_scale, P_-_scale, Q_+_scale, Q_-_scale], []] when bad data present, 
otherwise 'nothing'; the measurements are scaled through scale variables to make bad data.
'''
'''
# Output
- 'Bus_Measurement_Array': Bus measurement information in Array format with bad data injected -
[['Bus Number', 'V Measurement', 'P Measurement', 'Q Measurement',
Error Variance V, Error Variance P, Error Variance Q],...] if
measurement not present -9999 should be added.
- 'Branch_Measurement_Array': Branch measurement information in Array format with bad data injected -
[['Bus i Number', 'Bus j Number', 'P + Measurement', 'P - Measurement',
'Q + Measurement', 'Q - Measurement', Error Variance P Pos,
Error Variance P Neg, Error Variance Q Pos, Error Variance Q Neg],...]
if measurement not present -9999 should be added, '+' means i->j and '-' means
j->i.
'''
"""
function Create_Bad_BusBranch_Measurement_Arrays_SE(Bus_Measurement_Array, Branch_Measurement_Array, Bad_Bus_Measurement_Input, Bad_Branch_Measurement_Input)

        # If ElseIf Else Loop: For varying bad data types
        if ((Bad_Bus_Measurement_Input == nothing) && (Bad_Branch_Measurement_Input == nothing))  # No Bad Data

                Bad_Bus_Measurement_Array = copy(Bus_Measurement_Array)
                Bad_Branch_Measurement_Array = copy(Branch_Measurement_Array)

        elseif (Bad_Branch_Measurement_Input == nothing) # Only Bad Bus Data

                # Number of Buses and Lines
                N_Bus = size(Bus_Measurement_Array)[1]
                N_Lines = size(Branch_Measurement_Array)[1]  

                # Number of Bad Buses and Lines
                N_Bad_Bus = size(Bad_Bus_Measurement_Input)[1]
                N_Bad_Lines = size(Bad_Branch_Measurement_Input)[1]  

                # For Loop: For each bus in Bad_Bus_Measurement_Input
                for ii in 1:N_Bad_Bus

                        # Getting current bad bus number
                        Current_Bad_Bus_Num = Bad_Bus_Measurement_Input[ii,1]

                        # For Loop: For each bus in Bus_Measurement_Array
                        for jj in 1:N_Bus 

                                # If Loop: For checking Bus numbers
                                if (Current_Bad_Bus_Num == Bus_Measurement_Array[jj,1])

                                        # For Loop: For each measurement indicator in Bad_Bus_Measurement_Input
                                        for kk in 1:3

                                                if (Bad_Bus_Measurement_Input[ii,kk+1] != 1)

                                                        Bus_Measurement_Array[jj,kk+1] = Bad_Bus_Measurement_Input[ii,kk+1] * Bus_Measurement_Array[jj,kk+1]

                                                end

                                        end

                                end

                        end

                end

                Bad_Bus_Measurement_Array = copy(Bus_Measurement_Array)
                Bad_Branch_Measurement_Array = copy(Branch_Measurement_Array)

        elseif (Bad_Bus_Measurement_Input == nothing) # Only Bad Branch Data

                # Number of Buses and Lines
                N_Bus = size(Bus_Measurement_Array)[1]
                N_Lines = size(Branch_Measurement_Array)[1]  

                # Number of Bad Buses and Lines
                N_Bad_Bus = size(Bad_Bus_Measurement_Input)[1]
                N_Bad_Lines = size(Bad_Branch_Measurement_Input)[1] 

                # For Loop: For each line in in Bad_Branch_Measurement_Input
                for ii in 1:N_Bad_Lines
                                
                        # Getting current bad bus number
                        Current_Bad_Bus_Num_i = Bad_Branch_Measurement_Input[ii,1]
                        Current_Bad_Bus_Num_j = Bad_Branch_Measurement_Input[ii,2]
        
                        # For Loop: For each branch in Branch_Measurement_Array
                        for jj in 1:N_Lines 
        
                                # If Loop: For checking Branch numbers
                                if ((Current_Bad_Bus_Num_i == Branch_Measurement_Array[jj,1]) && (Current_Bad_Bus_Num_j == Branch_Measurement_Array[jj,2])) 
        
                                        # For Loop: For each measurement indicator in Bad_Branch_Measurement_Input
                                        for kk in 1:4
        
                                                if (Bad_Branch_Measurement_Input[ii,kk+2] != 1)
        
                                                        Branch_Measurement_Array[jj,kk+2] = Bad_Branch_Measurement_Input[ii,kk+2] * Branch_Measurement_Array[jj,kk+2]
        
                                                end
        
                                        end
        
                                end
        
                        end                               
        
                end
        
                Bad_Bus_Measurement_Array = copy(Bus_Measurement_Array)
                Bad_Branch_Measurement_Array = copy(Branch_Measurement_Array)

        else  # Bad Bus/Branch Data

                # Number of Buses and Lines
                N_Bus = size(Bus_Measurement_Array)[1]
                N_Lines = size(Branch_Measurement_Array)[1]  

                # Number of Bad Buses and Lines
                N_Bad_Bus = size(Bad_Bus_Measurement_Input)[1]
                N_Bad_Lines = size(Bad_Branch_Measurement_Input)[1]   

                # For Loop: For each bus in Bad_Bus_Measurement_Input
                for ii in 1:N_Bad_Bus

                        # Getting current bad bus number
                        Current_Bad_Bus_Num = Bad_Bus_Measurement_Input[ii,1]

                        # For Loop: For each bus in Bus_Measurement_Array
                        for jj in 1:N_Bus 

                                # If Loop: For checking Bus numbers
                                if (Current_Bad_Bus_Num == Bus_Measurement_Array[jj,1])

                                        # For Loop: For each measurement indicator in Bad_Bus_Measurement_Input
                                        for kk in 1:3

                                                if (Bad_Bus_Measurement_Input[ii,kk+1] != 1)

                                                        Bus_Measurement_Array[jj,kk+1] = Bad_Bus_Measurement_Input[ii,kk+1] * Bus_Measurement_Array[jj,kk+1]

                                                end

                                        end

                                end

                        end

                end

                # For Loop: For each line in in Bad_Branch_Measurement_Input
                for ii in 1:N_Bad_Lines
                                
                        # Getting current bad bus number
                        Current_Bad_Bus_Num_i = Bad_Branch_Measurement_Input[ii,1]
                        Current_Bad_Bus_Num_j = Bad_Branch_Measurement_Input[ii,2]

                        # For Loop: For each branch in Branch_Measurement_Array
                        for jj in 1:N_Lines 

                                # If Loop: For checking Branch numbers
                                if ((Current_Bad_Bus_Num_i == Branch_Measurement_Array[jj,1]) && (Current_Bad_Bus_Num_j == Branch_Measurement_Array[jj,2])) 

                                        # For Loop: For each measurement indicator in Bad_Branch_Measurement_Input
                                        for kk in 1:4

                                                if (Bad_Branch_Measurement_Input[ii,kk+2] != 1)

                                                        Branch_Measurement_Array[jj,kk+2] = Bad_Branch_Measurement_Input[ii,kk+2] * Branch_Measurement_Array[jj,kk+2]

                                                end

                                        end

                                end

                        end                               

                end

                Bad_Bus_Measurement_Array = copy(Bus_Measurement_Array)
                Bad_Branch_Measurement_Array = copy(Branch_Measurement_Array)

        end

        return Bad_Bus_Measurement_Array, Bad_Branch_Measurement_Array 

end

""" 
Compute_AMat_BadDataVec_SE(CDF_DF_List_pu)

Creates Initial Solution Vector for the Power System State Estimation.

'''     
# Arguments 

'''
'''
# Output

'''
"""
function Compute_AMat_BadDataVec_SE(CDF_DF_List_pu)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]
        BranchDataCard_DF = CDF_DF_List_pu[3]

        # Number of Buses and Lines
        N_Bus = nrow(BusDataCard_DF)
        N_Lines = nrow(BranchDataCard_DF)

        ## Creating Incidence Matrix

        # Initializing LineFlow_Array
        IncidenceMatrix_A = zeros(N_Lines, N_Bus)

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

                # Computing A_ij
                IncidenceMatrix_A[ii, Bus_i_Index] = 1
                IncidenceMatrix_A[ii, Bus_j_Index] = -1
                

        end


        ## Creating Initial Bad Data Vector
        V_Vec = ones(N_Bus, 1)
        P_Vec = ones(N_Bus, 1)
        Q_Vec = ones(N_Bus, 1)

        Pij_Vec = ones(N_Lines, 1)
        Pji_Vec = ones(N_Lines, 1)
        Qij_Vec = ones(N_Lines, 1)
        Qji_Vec = ones(N_Lines, 1)

        Detected_BadData_Vector_Ini = [V_Vec, P_Vec, Q_Vec, Pij_Vec, Pji_Vec, Qij_Vec, Qji_Vec,]

        return IncidenceMatrix_A, Detected_BadData_Vector_Ini 

end

""" 
Compute_R_Inv_Matrix_SE(CDF_DF_List_pu, Detected_BadData_Vector)

Creates Initial Solution Vector for the Power System State Estimation.

'''     
# Arguments 

'''
'''
# Output

'''
"""
function Compute_R_Inv_Matrix_SE(CDF_DF_List_pu, Detected_BadData_Vector)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]
        BranchDataCard_DF = CDF_DF_List_pu[3]

        # Number of Buses and Lines
        N_Bus = nrow(BusDataCard_DF)
        N_Lines = nrow(BranchDataCard_DF)

        # Computing Size of RInv_Matrix
        Size_R_Inv_Matrix = 0  # Initialization
        for ii in 1:length(Detected_BadData_Vector)  # For each element array in Detected_BadData_Vector

                for jj in 1:length(Detected_BadData_Vector[ii])  # For each element in Detected_BadData_Vector[ii]

                        # If Loop: To check for bad data
                        if (Detected_BadData_Vector[ii][jj] == 0)  # Bad Data present

                                # Do not Increment Size_R_Inv_Matrix

                        else  # Bad Data not present

                                # Increment Size_R_Inv_Matrix
                                Size_R_Inv_Matrix = Size_R_Inv_Matrix + 1

                        end

                end

        end

        # Initializing R_Inv_Matrix
        R_Inv_Matrix = zeros(Size_R_Inv_Matrix,Size_R_Inv_Matrix)

        # Filling up R_Inv_Matrix
        R_Inv_Matrix_Size_Counter = 0
        for ii in 1:length(Detected_BadData_Vector)  # For each element array in Detected_BadData_Vector

                for jj in 1:length(Detected_BadData_Vector[ii])  # For each element in Detected_BadData_Vector[ii]

                        # If Eleseif Loop: For different types of measurements
                        if (ii == 1)  # V Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment R_Inv_Matrix_Size_Counter
                                        R_Inv_Matrix_Size_Counter = R_Inv_Matrix_Size_Counter + 1

                                        # Filling the Diagonal element of R_Inv_Matrix_Size_Counter
                                        R_Inv_Matrix[R_Inv_Matrix_Size_Counter,R_Inv_Matrix_Size_Counter] = (1/(BusDataCard_DF.ErrorVariance_V[jj])^(2))

                                end

                        elseif (ii == 2)  # P Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment R_Inv_Matrix_Size_Counter
                                        R_Inv_Matrix_Size_Counter = R_Inv_Matrix_Size_Counter + 1

                                        # Filling the Diagonal element of R_Inv_Matrix_Size_Counter
                                        R_Inv_Matrix[R_Inv_Matrix_Size_Counter,R_Inv_Matrix_Size_Counter] = (1/(BusDataCard_DF.ErrorVariance_P[jj])^(2))

                                end

                        elseif (ii == 3)  # Q Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment R_Inv_Matrix_Size_Counter
                                        R_Inv_Matrix_Size_Counter = R_Inv_Matrix_Size_Counter + 1

                                        # Filling the Diagonal element of R_Inv_Matrix_Size_Counter
                                        R_Inv_Matrix[R_Inv_Matrix_Size_Counter,R_Inv_Matrix_Size_Counter] = (1/(BusDataCard_DF.ErrorVariance_Q[jj])^(2))

                                end

                        elseif (ii == 4)  # Pij Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment R_Inv_Matrix_Size_Counter
                                        R_Inv_Matrix_Size_Counter = R_Inv_Matrix_Size_Counter + 1

                                        # Filling the Diagonal element of R_Inv_Matrix_Size_Counter
                                        R_Inv_Matrix[R_Inv_Matrix_Size_Counter,R_Inv_Matrix_Size_Counter] = (1/(BranchDataCard_DF.ErrorVariance_Pos_P[jj])^(2))

                                end

                        elseif (ii == 5)  # Pji Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment R_Inv_Matrix_Size_Counter
                                        R_Inv_Matrix_Size_Counter = R_Inv_Matrix_Size_Counter + 1

                                        # Filling the Diagonal element of R_Inv_Matrix_Size_Counter
                                        R_Inv_Matrix[R_Inv_Matrix_Size_Counter,R_Inv_Matrix_Size_Counter] = (1/(BranchDataCard_DF.ErrorVariance_Neg_P[jj])^(2))

                                end

                        elseif (ii == 6)  # Qij Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment R_Inv_Matrix_Size_Counter
                                        R_Inv_Matrix_Size_Counter = R_Inv_Matrix_Size_Counter + 1

                                        # Filling the Diagonal element of R_Inv_Matrix_Size_Counter
                                        R_Inv_Matrix[R_Inv_Matrix_Size_Counter,R_Inv_Matrix_Size_Counter] = (1/(BranchDataCard_DF.ErrorVariance_Pos_Q[jj])^(2))

                                end

                        elseif (ii == 7)  # Qji Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment R_Inv_Matrix_Size_Counter
                                        R_Inv_Matrix_Size_Counter = R_Inv_Matrix_Size_Counter + 1

                                        # Filling the Diagonal element of R_Inv_Matrix_Size_Counter
                                        R_Inv_Matrix[R_Inv_Matrix_Size_Counter,R_Inv_Matrix_Size_Counter] = (1/(BranchDataCard_DF.ErrorVariance_Neg_Q[jj])^(2))

                                end

                        end

                end

        end


        return R_Inv_Matrix

end

""" 
Compute_Z_Measured_Vector_SE(CDF_DF_List_pu, Detected_BadData_Vector)

Creates Initial Solution Vector for the Power System State Estimation.

'''     
# Arguments 

'''
'''
# Output

'''
"""
function Compute_Z_Measured_Vector_SE(CDF_DF_List_pu, Detected_BadData_Vector)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]
        BranchDataCard_DF = CDF_DF_List_pu[3]

        # Number of Buses and Lines
        N_Bus = nrow(BusDataCard_DF)
        N_Lines = nrow(BranchDataCard_DF)

        # Computing Size of Z_Measured_Vector
        Size_Z_Measured_Vector = 0  # Initialization
        for ii in 1:length(Detected_BadData_Vector)  # For each element array in Detected_BadData_Vector

                for jj in 1:length(Detected_BadData_Vector[ii])  # For each element in Detected_BadData_Vector[ii]

                        # If Loop: To check for bad data
                        if (Detected_BadData_Vector[ii][jj] == 0)  # Bad Data present

                                # Do not Increment Size_R_Inv_Matrix

                        else  # Bad Data not present

                                # Increment Size_R_Inv_Matrix
                                Size_Z_Measured_Vector = Size_Z_Measured_Vector + 1

                        end

                end

        end

        # Initializing Z_Measured_Vector
        Z_Measured_Vector = zeros(Size_Z_Measured_Vector,1)

        # Filling up Z_Measured_Vector
        Z_Measured_Vector_Size_Counter = 0
        for ii in 1:length(Detected_BadData_Vector)  # For each element array in Detected_BadData_Vector

                for jj in 1:length(Detected_BadData_Vector[ii])  # For each element in Detected_BadData_Vector[ii]

                        # If Eleseif Loop: For different types of measurements
                        if (ii == 1)  # V Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment Z_Measured_Vector_Size_Counter
                                        Z_Measured_Vector_Size_Counter = Z_Measured_Vector_Size_Counter + 1

                                        # Filling the Diagonal element of Z_Measured_Vector
                                        Z_Measured_Vector[Z_Measured_Vector_Size_Counter,1] = BusDataCard_DF.V_Measured[jj]

                                end

                        elseif (ii == 2)  # P Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment Z_Measured_Vector_Size_Counter
                                        Z_Measured_Vector_Size_Counter = Z_Measured_Vector_Size_Counter + 1

                                        # Filling the Diagonal element of Z_Measured_Vector
                                        Z_Measured_Vector[Z_Measured_Vector_Size_Counter,1] = BusDataCard_DF.P_Measured[jj]

                                end

                        elseif (ii == 3)  # Q Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment Z_Measured_Vector_Size_Counter
                                        Z_Measured_Vector_Size_Counter = Z_Measured_Vector_Size_Counter + 1

                                        # Filling the Diagonal element of Z_Measured_Vector
                                        Z_Measured_Vector[Z_Measured_Vector_Size_Counter,1] = BusDataCard_DF.Q_Measured[jj]

                                end

                        elseif (ii == 4)  # Pij Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment Z_Measured_Vector_Size_Counter
                                        Z_Measured_Vector_Size_Counter = Z_Measured_Vector_Size_Counter + 1

                                        # Filling the Diagonal element of Z_Measured_Vector
                                        Z_Measured_Vector[Z_Measured_Vector_Size_Counter,1] = BranchDataCard_DF.Line_Flow_Pos_P[jj]

                                end

                        elseif (ii == 5)  # Pji Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment Z_Measured_Vector_Size_Counter
                                        Z_Measured_Vector_Size_Counter = Z_Measured_Vector_Size_Counter + 1

                                        # Filling the Diagonal element of Z_Measured_Vector
                                        Z_Measured_Vector[Z_Measured_Vector_Size_Counter,1] = BranchDataCard_DF.Line_Flow_Neg_P[jj]

                                end

                        elseif (ii == 6)  # Qij Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment Z_Measured_Vector_Size_Counter
                                        Z_Measured_Vector_Size_Counter = Z_Measured_Vector_Size_Counter + 1

                                        # Filling the Diagonal element of Z_Measured_Vector
                                        Z_Measured_Vector[Z_Measured_Vector_Size_Counter,1] = BranchDataCard_DF.Line_Flow_Pos_Q[jj]

                                end

                        elseif (ii == 7)  # Qji Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment Z_Measured_Vector_Size_Counter
                                        Z_Measured_Vector_Size_Counter = Z_Measured_Vector_Size_Counter + 1

                                        # Filling the Diagonal element of Z_Measured_Vector
                                        Z_Measured_Vector[Z_Measured_Vector_Size_Counter,1] = BranchDataCard_DF.Line_Flow_Neg_Q[jj]

                                end

                        end

                end

        end
        
        return Z_Measured_Vector 

end

""" 
Compute_Z_Calculated_Vector_SE(SolutionVector_V, SolutionVector_Delta, Ybus, Detected_BadData_Vector)

Creates Initial Solution Vector for the Power System State Estimation.

'''     
# Arguments 

'''
'''
# Output

'''
"""
function Compute_Z_Calculated_Vector_SE(CDF_DF_List_pu, SolutionVector_V, SolutionVector_Delta, Ybus, IncidenceMatrix_A, Detected_BadData_Vector)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]
        BranchDataCard_DF = CDF_DF_List_pu[3]

        # Number of Buses and Lines
        N_Bus = nrow(BusDataCard_DF)
        N_Lines = nrow(BranchDataCard_DF)

        # Computing Size of Z_Calculated_Vector
        Size_Z_Calculated_Vector = 0  # Initialization
        for ii in 1:length(Detected_BadData_Vector)  # For each element array in Detected_BadData_Vector

                for jj in 1:length(Detected_BadData_Vector[ii])  # For each element in Detected_BadData_Vector[ii]

                        # If Loop: To check for bad data
                        if (Detected_BadData_Vector[ii][jj] == 0)  # Bad Data present

                                # Do not Increment Size_R_Inv_Matrix

                        else  # Bad Data not present

                                # Increment Size_R_Inv_Matrix
                                Size_Z_Calculated_Vector = Size_Z_Calculated_Vector + 1

                        end

                end

        end

        # Initializing Z_MCalculated_Vector
        Z_Calculated_Vector = zeros(Size_Z_Calculated_Vector,1)

        # Filling up Z_Calculated_Vector
        Z_Calculated_Vector_Size_Counter = 0
        for ii in 1:length(Detected_BadData_Vector)  # For each element array in Detected_BadData_Vector

                for jj in 1:length(Detected_BadData_Vector[ii])  # For each element in Detected_BadData_Vector[ii]

                        # If Eleseif Loop: For different types of measurements
                        if (ii == 1)  # V Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment Z_Calculated_Vector_Size_Counter
                                        Z_Calculated_Vector_Size_Counter = Z_Calculated_Vector_Size_Counter + 1

                                        # Filling the Diagonal element of Z_Calculated_Vector
                                        Z_Calculated_Vector[Z_Calculated_Vector_Size_Counter,1] = SolutionVector_V[jj,1]

                                end

                        elseif (ii == 2)  # P Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment Z_Calculated_Vector_Size_Counter
                                        Z_Calculated_Vector_Size_Counter = Z_Calculated_Vector_Size_Counter + 1

                                        # Filling the Diagonal element of Z_Calculated_Vector
                                        for kk in 1:length(SolutionVector_V)  # For each Node in the system

                                                Z_Calculated_Vector[Z_Calculated_Vector_Size_Counter,1] = Z_Calculated_Vector[Z_Calculated_Vector_Size_Counter,1] + (SolutionVector_V[jj]*SolutionVector_V[kk]*abs(Ybus[jj,kk])*cos(angle(Ybus[jj,kk])+deg2rad(SolutionVector_Delta[kk])-deg2rad(SolutionVector_Delta[jj])))

                                        end

                                end

                        elseif (ii == 3)  # Q Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment Z_Calculated_Vector_Size_Counter
                                        Z_Calculated_Vector_Size_Counter = Z_Calculated_Vector_Size_Counter + 1

                                        # Filling the Diagonal element of Z_Calculated_Vector
                                        for kk in 1:length(SolutionVector_V)  # For each Node in the system

                                                Z_Calculated_Vector[Z_Calculated_Vector_Size_Counter,1] = Z_Calculated_Vector[Z_Calculated_Vector_Size_Counter,1] - (SolutionVector_V[jj]*SolutionVector_V[kk]*abs(Ybus[jj,kk])*sin(angle(Ybus[jj,kk])+deg2rad(SolutionVector_Delta[kk])-deg2rad(SolutionVector_Delta[jj])))

                                        end

                                end

                        elseif (ii == 4)  # Pij Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment Z_Calculated_Vector_Size_Counter
                                        Z_Calculated_Vector_Size_Counter = Z_Calculated_Vector_Size_Counter + 1

                                        # Getting current measurement node indices 
                                        Current_A_Row = IncidenceMatrix_A[jj,:]  
                                        
                                        Bus_i_Index = 0
                                        Bus_j_Index = 0
                                        
                                        for kk in 1:length(Current_A_Row)  # For Bus_i_Index

                                                if (Current_A_Row[kk] == 1)

                                                        Bus_i_Index = kk

                                                        break

                                                end

                                        end

                                        for kk in 1:length(Current_A_Row)  # For Bus_j_Index

                                                if (Current_A_Row[kk] == -1)

                                                        Bus_j_Index = kk

                                                        break

                                                end

                                        end

                                        # Filling the Diagonal element of Z_Calculated_Vector
                                        Z_Calculated_Vector[Z_Calculated_Vector_Size_Counter,1] = -(SolutionVector_V[Bus_i_Index]^(2)*real(Ybus[Bus_i_Index, Bus_j_Index])) +(SolutionVector_V[Bus_i_Index]*SolutionVector_V[Bus_j_Index]*abs(Ybus[Bus_i_Index, Bus_j_Index]))*cos(angle(Ybus[Bus_i_Index, Bus_j_Index])+deg2rad(SolutionVector_Delta[Bus_j_Index])-deg2rad(SolutionVector_Delta[Bus_i_Index]))

                                end

                        elseif (ii == 5)  # Pji Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment Z_Calculated_Vector_Size_Counter
                                        Z_Calculated_Vector_Size_Counter = Z_Calculated_Vector_Size_Counter + 1

                                        # Getting current measurement node indices 
                                        Current_A_Row = IncidenceMatrix_A[jj,:]  
                                        
                                        Bus_i_Index = 0
                                        Bus_j_Index = 0

                                        for kk in 1:length(Current_A_Row)  # For Bus_i_Index

                                                if (Current_A_Row[kk] == 1)

                                                        Bus_i_Index = kk

                                                        break

                                                end

                                        end

                                        for kk in 1:length(Current_A_Row)  # For Bus_j_Index

                                                if (Current_A_Row[kk] == -1)

                                                        Bus_j_Index = kk

                                                        break

                                                end

                                        end

                                        # Filling the Diagonal element of Z_Calculated_Vector
                                        Z_Calculated_Vector[Z_Calculated_Vector_Size_Counter,1] = -(SolutionVector_V[Bus_j_Index]^(2)*real(Ybus[Bus_j_Index, Bus_i_Index])) +(SolutionVector_V[Bus_j_Index]*SolutionVector_V[Bus_i_Index]*abs(Ybus[Bus_j_Index, Bus_i_Index]))*cos(angle(Ybus[Bus_j_Index, Bus_i_Index])+deg2rad(SolutionVector_Delta[Bus_i_Index])-deg2rad(SolutionVector_Delta[Bus_j_Index]))

                                end

                        elseif (ii == 6)  # Qij Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment Z_Calculated_Vector_Size_Counter
                                        Z_Calculated_Vector_Size_Counter = Z_Calculated_Vector_Size_Counter + 1

                                        # Getting current measurement node indices 
                                        Current_A_Row = IncidenceMatrix_A[jj,:] 
                                        
                                        Bus_i_Index = 0
                                        Bus_j_Index = 0

                                        for kk in 1:length(Current_A_Row) # For Bus_i_Index

                                                if (Current_A_Row[kk] == 1)

                                                        Bus_i_Index = kk

                                                        break

                                                end

                                        end

                                        for kk in 1:length(Current_A_Row)  # For Bus_j_Index

                                                if (Current_A_Row[kk] == -1)

                                                        Bus_j_Index = kk

                                                        break

                                                end

                                        end

                                        # Filling the Diagonal element of Z_Calculated_Vector
                                        Z_Calculated_Vector[Z_Calculated_Vector_Size_Counter,1] = -((SolutionVector_V[Bus_i_Index]^(2)*((BranchDataCard_DF.B_pu[jj]/2) - imag(Ybus[Bus_i_Index, Bus_j_Index])))) - (SolutionVector_V[Bus_i_Index]*SolutionVector_V[Bus_j_Index]*abs(Ybus[Bus_i_Index, Bus_j_Index]))*sin(angle(Ybus[Bus_i_Index, Bus_j_Index])+deg2rad(SolutionVector_Delta[Bus_j_Index])-deg2rad(SolutionVector_Delta[Bus_i_Index]))

                                end

                        elseif (ii == 7)  # Qji Vector

                                if (Detected_BadData_Vector[ii][jj] != 0)  # Bad data not present

                                        # Increment Z_Calculated_Vector_Size_Counter
                                        Z_Calculated_Vector_Size_Counter = Z_Calculated_Vector_Size_Counter + 1

                                        # Getting current measurement node indices 
                                        Current_A_Row = IncidenceMatrix_A[jj,:]  
                                        
                                        Bus_i_Index = 0
                                        Bus_j_Index = 0

                                        for kk in 1:length(Current_A_Row)  # For Bus_i_Index

                                                if (Current_A_Row[kk] == 1)

                                                        Bus_i_Index = kk

                                                        break

                                                end

                                        end

                                        for kk in 1:length(Current_A_Row)  # For Bus_j_Index

                                                if (Current_A_Row[kk] == -1)

                                                        Bus_j_Index = kk

                                                        break

                                                end

                                        end

                                        # Filling the Diagonal element of Z_Calculated_Vector
                                        Z_Calculated_Vector[Z_Calculated_Vector_Size_Counter,1] = -((SolutionVector_V[Bus_j_Index]^(2)*((BranchDataCard_DF.B_pu[jj]/2) - imag(Ybus[Bus_j_Index, Bus_i_Index])))) - (SolutionVector_V[Bus_j_Index]*SolutionVector_V[Bus_i_Index]*abs(Ybus[Bus_j_Index, Bus_i_Index]))*sin(angle(Ybus[Bus_j_Index, Bus_i_Index])+deg2rad(SolutionVector_Delta[Bus_i_Index])-deg2rad(SolutionVector_Delta[Bus_j_Index]))

                                end

                        end

                end

        end

        # Addressing Machine Precision Problem
        for ii in 1:length(Z_Calculated_Vector)

                if (abs(Z_Calculated_Vector[ii,1]) < 1e-12)
                        
                        Z_Calculated_Vector[ii,1] = 0

                end

        end

        return Z_Calculated_Vector  

end

""" 
Compute_H_Matrix_SE(CDF_DF_List_pu, SolutionVector_V, SolutionVector_Delta, Ybus, IncidenceMatrix_A, Detected_BadData_Vector)

Creates Initial Solution Vector for the Power System State Estimation.

'''     
# Arguments 

'''
'''
# Output

'''
"""
function Compute_H_Matrix_SE(CDF_DF_List_pu, SolutionVector_V, SolutionVector_Delta, Ybus, IncidenceMatrix_A, Detected_BadData_Vector)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]
        BranchDataCard_DF = CDF_DF_List_pu[3]

        # Number of Buses and Lines
        N_Bus = nrow(BusDataCard_DF)
        N_Lines = nrow(BranchDataCard_DF)       
        
        return H_Matrix 

end

""" 
Compute_Corrected_MisMatchVector_SE(CDF_DF_List_pu, StateEstimation_MisMatch_Vector)

Creates Initial Solution Vector for the Power System State Estimation.

'''     
# Arguments 

'''
'''
# Output

'''
"""
function Compute_Corrected_MisMatchVector_SE(CDF_DF_List_pu, StateEstimation_MisMatch_Vector)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses and Lines
        N_Bus = nrow(BusDataCard_DF)

        # Initializing
        StateEstimation_MisMatch_Vector_Corrected = copy(StateEstimation_MisMatch_Vector)

        # For Loop: For each delta in StateEstimation_MisMatch_Vector
        for ii in 1:(N_Bus - 1)

                # Correcting delta from Radian to Degree
                StateEstimation_MisMatch_Vector_Corrected[ii] = rad2deg(StateEstimation_MisMatch_Vector_Corrected[ii])

        end

        return StateEstimation_MisMatch_Vector_Corrected

end

""" 
Compute_ToleranceSatisfaction_SE(StateEstimation_MisMatch_Vector, Tolerance_SE)

Creates Initial Solution Vector for the Power System State Estimation.

'''     
# Arguments 

'''
'''
# Output

'''
"""
function Compute_ToleranceSatisfaction_SE(StateEstimation_MisMatch_Vector, Tolerance_SE)

        # Initializing Tolerance Counter
        ToleranceCounter = 0

        # Computing Tolerance Counter
        for ii in 1: length(StateEstimation_MisMatch_Vector)

                if (abs(StateEstimation_MisMatch_Vector[ii]) > Tolerance)

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

        return Tolerance_Satisfaction_SE 

end

""" 
Compute_StateEstimation(CDF_DF_List_pu, SolutionVector_SE, Detected_Bad_Data)

Creates Initial Solution Vector for the Power System State Estimation.

'''     
# Arguments 

'''
'''
# Output

'''
"""
function Compute_StateEstimation(CDF_DF_List_pu, SolutionVector_SE_Ini, Ybus, IncidenceMatrix_A, Detected_BadData_Vector, Tolerance_SE)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]
        BranchDataCard_DF = CDF_DF_List_pu[3]

        # Number of Buses and Lines
        N_Bus = nrow(BusDataCard_DF)
        N_Lines = nrow(BranchDataCard_DF)

        # Number of Buses
        N_PQ_Bus = nrow(filter(row -> ((row.Type_Original == 0) || (row.Type_Original == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type_Original == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type_Original == 3), BusDataCard_DF))        

        # Create R Inverse Matrix - Based on Detected Bad Data
        R_Inv_Matrix = Compute_R_Inv_Matrix_SE(CDF_DF_List_pu, Detected_BadData_Vector)

        # Create Z Vector Measured - Based on Detected Bad Data
        Z_Measured_Vector = Compute_Z_Measured_Vector_SE(CDF_DF_List_pu, Detected_BadData_Vector)

        # Initializing Tolerance_Satisfaction
        Tolerance_Satisfaction = false

        # InItializing 
        WhileLoop_Counter = 0

        SolutionVector_SE = zeros(N_Bus-1, 1)

        # While Loop: For performing State Estimation
        while (!Tolerance_Satisfaction)

                # Incrementing WhileLoop_Counter
                WhileLoop_Counter = WhileLoop_Counter + 1

                 # If Else Loop: For checking First Iteration
                if (WhileLoop_Counter == 1)

                        # Getting Solution Vectors divided in Voltages and Deltas
                        SolutionVector_V, SolutionVector_Delta = Create_SolutionVector_VDelta_SE(CDF_DF_List_pu, SolutionVector_SE_Ini)
                
                else

                        # Getting Solution Vectors divided in Voltages and Deltas
                        SolutionVector_V, SolutionVector_Delta = Create_SolutionVector_VDelta_SE(CDF_DF_List_pu, SolutionVector_SE)
                
                end

                # Compute Z Vector for calculated measurements - Based on Detected Bad Data
                Z_Calculated_Vector  = Compute_Z_Calculated_Vector_SE(CDF_DF_List_pu, SolutionVector_V, SolutionVector_Delta, Ybus, IncidenceMatrix_A, Detected_BadData_Vector)

                # Compute the H Matrix - Based on Detected Bad Data
                H_Matrix = Compute_H_Matrix_SE(CDF_DF_List_pu, SolutionVector_V, SolutionVector_Delta, Ybus, IncidenceMatrix_A, Detected_BadData_Vector)

                # Compute MisMatch Vector
                MisMatch_Vector = Z_Measured_Vector - Z_Calculated_Vector 

                # Compute G1 Matrix
                G1_Matrix = H_Matrix' * RInv_Matrix * H_Matrix

                # Compute G2 Matrix
                G2_Matrix = H_Matrix' * RInv_Matrix * MisMatch_Vector

                # Computing State Estimation Mismatch Vector
                StateEstimation_MisMatch_Vector = G1_Matrix \ G2_Matrix

                # Correcting StateEstimation_MisMatch_Vector for Radians to Degree
                StateEstimation_MisMatch_Vector_Corrected = Compute_Corrected_MisMatchVector_SE(CDF_DF_List_pu, StateEstimation_MisMatch_Vector)

                # Updating State Estimation Vector
                SolutionVector_SE = SolutionVector_SE + StateEstimation_MisMatch_Vector_Corrected

                # Compute Tolerance Satisfaction
                Tolerance_Satisfaction = Compute_ToleranceSatisfaction_SE(StateEstimation_MisMatch_Vector_Corrected, Tolerance_SE)

        end

        return SolutionVector_SE, CDF_DF_List_pu, RInv_Matrix 

end