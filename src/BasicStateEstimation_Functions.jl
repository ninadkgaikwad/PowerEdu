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
[Error Variance V, Error Variance P, Error Variance Q]],...] if
measurement not present -9999 should be added.
- 'Branch_Measurement_Array': Branch measurement information in Array format -
[['Bus i Number', 'Bus j Number', 'P + Measurement', 'P - Measurement',
'Q + Measurement', 'Q - Measurement', [Error Variance P Pos,
Error Variance P Neg, Error Variance Q Pos, Error Variance Q Neg]],...]
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

        # For Loop: For adding Bus Measurements
        for ii in 1:length(Bus_Measurement_Array)

                # Getting Current Bus Measurements
                Current_Bus_Measurement = Bus_Measurement_Array[ii]

                # Getting Components of Current_Bus_Measurement
                Current_BusNumber = Current_Bus_Measurement[1]
                Current_V_Measured = Current_Bus_Measurement[2]
                Current_P_Measured = Current_Bus_Measurement[3]
                Current_Q_Measured = Current_Bus_Measurement[4]
                Current_V_EV = Current_Bus_Measurement[5][1]
                Current_P_EV = Current_Bus_Measurement[5][2]
                Current_Q_EV = Current_Bus_Measurement[5][3]

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

                                        BusDataCard_DF.ErrorVariance_V [jj] = Current_V_EV

                                end

                                if (Current_P_EV != -9999)

                                        BusDataCard_DF.ErrorVariance_P [jj] = Current_P_EV

                                end

                                if (Current_Q_EV != -9999)

                                        BusDataCard_DF.ErrorVariance_Q [jj] = Current_Q_EV

                                end

                        end

                end

        end

        # For Loop: For adding Branch Measurements
        for ii in 1:length(Branch_Measurement_Array)

                # Getting Current Bus Measurements
                Current_Branch_Measurement = Branch_Measurement_Array[ii]

                # Getting Components of Current_Bus_Measurement
                Current_BusNumber_i = Current_Bus_Measurement[1]
                Current_BusNumber_j = Current_Bus_Measurement[2]
                Current_P_Pos_Measured = Current_Bus_Measurement[3]
                Current_P_Neg_Measured = Current_Bus_Measurement[4]
                Current_Q_Pos_Measured = Current_Bus_Measurement[5]
                Current_Q_Neg_Measured = Current_Bus_Measurement[6]
                Current_P_Pos_EV = Current_Bus_Measurement[7][1]
                Current_P_Neg_EV = Current_Bus_Measurement[7][2]
                Current_Q_Pos_EV = Current_Bus_Measurement[7][3]
                Current_Q_Neg_EV = Current_Bus_Measurement[7][4]

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

                                        BranchDataCard_DF.ErrorVariance_Neg_Pjj] = Current_P_Neg_EV

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

                                        BranchDataCard_DF.ErrorVariance_Pos_Pjj] = Current_P_Neg_EV

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

                        BusDataCard_DF.V_Measured[ii] = BusDataCard_DF.V_Measured[ii]/BusDataCard_DF.Base_KV [ii]

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
