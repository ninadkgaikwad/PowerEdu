# BasicContinuationPowerFlow_Functions.jl


"""
    Create_Initial_SolutionVector_CPF(CDF_DF_List_pu)

Creates Initial Solution Vector for the Continuation Power Flow Method.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
'''
'''
# Output
- 'Initial_SolutionVector_CPF': Power Flow Solution Voltage and Angle at each
bus ordered according to bus type: PQ->PV->Lambda.
'''
"""
function Create_Initial_SolutionVector_CPF(CDF_DF_List_pu)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type_Original == 0) || (row.Type_Original == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type_Original == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type_Original == 3), BusDataCard_DF))

        # Creating Initial_SolutionVector_CPF Length
        Initial_SolutionVector_CPF_Len = 2*N_PQ_Bus+N_PV_Bus+1

        # Initializing Initial_SolutionVector_CPF
        Initial_SolutionVector_CPF = Array{Float64}(undef, Initial_SolutionVector_CPF_Len ,1)

        # Creating Initial_SolutionVector_CPF Delta part for both PQ and PV Buses
        Initial_SolutionVector_CPF[1:(N_PQ_Bus+N_PV_Bus),1] = BusDataCard_DF.Final_A_deg[1:end-1]

        # Creating Initial_SolutionVector_CPF V part for PQ Buses
        Initial_SolutionVector_CPF[(N_PQ_Bus+N_PV_Bus+1):end-1,1] = BusDataCard_DF.Final_V_pu[1:N_PQ_Bus]

        # Creating Initial_SolutionVector_CPF Lambda part
        Initial_SolutionVector_CPF[end,1] = 0

        return Initial_SolutionVector_CPF

end

"""
    Create_SolutionVector_VDelta_CPF(CDF_DF_List_pu, SolutionVector_NR)

Creates Initial Solution Vector for the Continuation Power Flow.

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
function Create_SolutionVector_VDelta_CPF(CDF_DF_List_pu, SolutionVector_CPF)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type_Original == 0) || (row.Type_Original == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type_Original == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type_Original == 3), BusDataCard_DF))

        # Getting SolutionVector_NR
        SolutionVector_NR = SolutionVector_CPF[1:end-1,1]

        # Initializing SolutionVector_V, SolutionVector_Delta
        SolutionVector_V = Array{Float64}(undef, N_Bus,1)
        SolutionVector_Delta = Array{Float64}(undef, N_Bus,1)

        # Creating SolutionVector_V, SolutionVector_Delta
        SolutionVector_Delta =  vcat([0.0], SolutionVector_NR[1:(N_Bus-1),1])
        SolutionVector_V[1:(N_PQ_Bus+1),1] =  vcat([BusDataCard_DF.Final_V_pu_Original[end]], SolutionVector_NR[(N_Bus-1)+1:end,1])

        for ii in N_PQ_Bus+1+1:N_Bus # For each current PV Bus

                SolutionVector_V[ii,1] = BusDataCard_DF.Final_V_pu_Original[ii-1]

        end

        return SolutionVector_V, SolutionVector_Delta

end

"""
    Create_SolutionVector_CPF(CDF_DF_List_pu, SolutionVector_V, SolutionVector_Delta)

Creates Initial Solution Vector for the Continuation Power Flow.

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
function Create_SolutionVector_CPF(CDF_DF_List_pu, SolutionVector_V, SolutionVector_Delta)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type_Original == 0) || (row.Type_Original == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type_Original == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type_Original == 3), BusDataCard_DF))

        # Creating SolutionVector_NR Length
        SolutionVector_NR_Len = 2*N_PQ_Bus+N_PV_Bus

        # Initializing SolutionVector_NR
        SolutionVector_NR = Array{Float64}(undef, Initial_SolutionVector_NR_Len ,1)

        # Creating SolutionVector_NR Delta part for both PQ and PV Buses
        SolutionVector_NR[1:(N_PQ_Bus+N_PV_Bus),1] = SolutionVector_Delta[2:end,1]

        # Creating Initial_SolutionVector_NR V part for PQ Buses
        SolutionVector_NR = vcat(SolutionVector_NR, SolutionVector_V[2:N_PQ_Bus+1,1])

        return SolutionVector_NR

end

"""
    Create_KVector_CPF(CDF_DF_List_pu, PQ_V_Curve_Tuple)

Creates Initial Solution Vector for the Continuation Power Flow Method.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
-'PQ_V_Curve_Tuple': A Tuple containing (BusNumber,P/Q to be parametrized)
'''
'''
# Output
- 'K_Vector': Consists of a Vector all zero except for the bus used for
parameterization filled with actual either P_i or Q_i.
'''
"""
function Create_KVector_CPF(CDF_DF_List_pu, PQ_V_Curve_Tuple)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type_Original == 0) || (row.Type_Original == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type_Original == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type_Original == 3), BusDataCard_DF))

        # Initializing K_Vector
        K_Vector_Len = 2*N_PQ_Bus+N_PV_Bus
        K_Vector = zeros(K_Vector_Len,1)

        # Initialization
        Bus_Plot_Index = 0

        # For Loop: For filling up K_Vector
        for ii in 1:N_Bus-1

                if (BusDataCard_DF.Bus_Num[ii] == PQ_V_Curve_Tuple[1])

                        if (PQ_V_Curve_Tuple[2] == 'P')

                                K_Vector[ii] = (BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii])

                        elseif (PQ_V_Curve_Tuple[2] == 'Q')

                                K_Vector[N_PQ_Bus+N_PV_Bus+ii] = (BusDataCard_DF.Gen_MVAR[ii] - BusDataCard_DF.Load_MVAR[ii])

                        end

                end

        end

        for ii in 1:N_Bus-1

                if (BusDataCard_DF.Bus_Num[ii] == PQ_V_Curve_Tuple[1])

                        Bus_Plot_Index = N_PQ_Bus+N_PV_Bus + ii
                        
                end

        end

        return K_Vector, Bus_Plot_Index 

end

"""
    Create_KVector_1_CPF(CDF_DF_List_pu, PQ_V_Curve_Tuple)

Creates Initial Solution Vector for the Continuation Power Flow Method.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
-'PQ_V_Curve_Tuple': A Tuple containing (BusNumber,P/Q to be parametrized)
'''
'''
# Output
- 'K_Vector': Consists of a Vector all zero except for the bus used for
parameterization filled with actual either P_i or Q_i.
'''
"""
function Create_KVector_1_CPF(CDF_DF_List_pu, PQ_V_Curve_Tuple)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type_Original == 0) || (row.Type_Original == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type_Original == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type_Original == 3), BusDataCard_DF))

        # Initializing K_Vector
        K_Vector_Len = 2*N_PQ_Bus+N_PV_Bus
        K_Vector = zeros(K_Vector_Len,1)

        # Initialization
        Bus_Plot_Index = 0

        # For Loop: For filling up K_Vector
        for ii in 1:N_Bus-1

                if (BusDataCard_DF.Bus_Num[ii] == PQ_V_Curve_Tuple[1])

                        K_Vector[ii] = (BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii])

                        K_Vector[N_PQ_Bus+N_PV_Bus+ii] = (BusDataCard_DF.Gen_MVAR[ii] - BusDataCard_DF.Load_MVAR[ii])

                end

        end

        for ii in 1:N_Bus-1

                if (BusDataCard_DF.Bus_Num[ii] == PQ_V_Curve_Tuple[1])

                        Bus_Plot_Index = N_PQ_Bus+N_PV_Bus + ii

                end

        end

        return K_Vector, Bus_Plot_Index

end

"""
    Create_KVector_2_CPF(CDF_DF_List_pu, PQ_V_Curve_Tuple)

Creates Initial Solution Vector for the Continuation Power Flow Method.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
-'PQ_V_Curve_Tuple': A Tuple containing (BusNumber,P/Q to be parametrized)
'''
'''
# Output
- 'K_Vector': Consists of a Vector all zero except for the bus used for
parameterization filled with actual either P_i or Q_i.
'''
"""
function Create_KVector_2_CPF(CDF_DF_List_pu, PQ_V_Curve_Tuple)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type_Original == 0) || (row.Type_Original == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type_Original == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type_Original == 3), BusDataCard_DF))

        # Initializing K_Vector
        K_Vector_Len = 2*N_PQ_Bus+N_PV_Bus
        K_Vector = ones(K_Vector_Len,1)

        # Initialization
        Bus_Plot_Index = 0

        # For Loop: For filling up K_Vector
        for ii in 1:K_Vector_Len
               
                if (ii <= (N_PQ_Bus+N_PV_Bus))

                        K_Vector[ii] = (BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii])

                else

                        K_Vector[ii] = (BusDataCard_DF.Gen_MVAR[ii-(N_PQ_Bus+N_PV_Bus)] - BusDataCard_DF.Load_MVAR[ii-(N_PQ_Bus+N_PV_Bus)])

                end
                
        end

        for ii in 1:N_Bus-1

                if (BusDataCard_DF.Bus_Num[ii] == PQ_V_Curve_Tuple[1])

                        Bus_Plot_Index = N_PQ_Bus+N_PV_Bus + ii
                        
                end

        end

        return K_Vector, Bus_Plot_Index

end

"""
    Compute_Tangent_Vector_CPF(Jacobian_CPF_Predict, Tol_Num,
    Predict_b_Value_CPF, Index_CPF)

Computes Tangent Vector for the Predict Step of Continuation Power Flow.

'''
# Arguments
- 'Jacobian_CPF_Predict': Jacobian Matrix for Predict Step of Continuation Power
Flow: PQ->PV->Lambda.
- 'Tol_Num': Tolerance for being near zero
- 'Predict_b_Value_CPF': Scalar +1 or -1 based on choice of Continuation
Parameter
- 'Index_CPF': Continuation Parameter Index
'''
'''
# Output
- 'Tangent_Vector': Tangent Vector for Predict Step of Continuation Power
Flow: PQ->PV->Lambda.
'''
"""
function Compute_Tangent_Vector_CPF(Jacobian_CPF_Predict, Tol_Num, Predict_b_Value_CPF, Index_CPF)

        # Creating b Vector for Predictor Step
        Predict_b_Vector = zeros(size(Jacobian_CPF_Predict)[1],1)

        Predict_b_Vector[Index_CPF,1] = Predict_b_Value_CPF

        # Computing Tangent Vector using LU Factorization
        Tangent_Vector = PLU_Solve(Jacobian_CPF_Predict, Predict_b_Vector, Tol_Num)

    return Tangent_Vector

end

"""
    Check_CriticalPoint_CPF(Tangent_Vector, PostCriticalPoint_Counter)

Checks passing of Critical Point for Continuation Power Flow.

'''
# Arguments
- 'Tangent_Vector': Tangent Vector for Predict Step of Continuation Power
Flow: PQ->PV->Lambda.
- 'PostCriticalPoint_Counter': Integer count of iterations after passing of
Critical Point
'''
'''
# Output
- 'PostCriticalPoint_Counter_New': Integer count of iterations after passing of
Critical Point for next iteration.
'''
"""
function Check_CriticalPoint_CPF(Tangent_Vector, PostCriticalPoint_Counter)

        Differential_Lambda = Tangent_Vector[end,1]

        if (Differential_Lambda < 0) # Lambda started to decrease

                # Incrementing PostCriticalPoint_Counter
                PostCriticalPoint_Counter = PostCriticalPoint_Counter + 1

        else# Lambda is increasing or stable

                # Reinitializing PostCriticalPoint_Counter
                PostCriticalPoint_Counter = 0

        end

    return PostCriticalPoint_Counter

end

"""
    Choose_ContinuationParameter_CPF(Tangent_Vector, PostCriticalPoint_Counter)

Selects Continuation Parameter for Continuation Power Flow.

'''
# Arguments
- 'Tangent_Vector': Tangent Vector for Predict Step of Continuation Power
Flow: PQ->PV->Lambda.
'''
'''
# Output
- 'Index_CPF': Continuation Parameter Index
- 'Predict_b_Value_CPF': Scalar +1 or -1 based on choice of Continuation
Parameter
'''
"""
function Choose_ContinuationParameter_CPF(Tangent_Vector)

        Abs_Tangent_Vector = broadcast(abs, Tangent_Vector)

        Abs_Tangent_Vector_Max, Index_CPF = findmax(Abs_Tangent_Vector)

        Index_CPF = Index_CPF[1]

        if (Tangent_Vector[Index_CPF,1] >= 0 ) # Parameter increasing or stable

                Predict_b_Value_CPF = 1

        else # Parameter decreasing

                Predict_b_Value_CPF = -1

        end

    return Index_CPF, Predict_b_Value_CPF

end

"""
    Choose_ContinuationParameter_CPF(Tangent_Vector, PostCriticalPoint_Counter)

Selects Continuation Parameter for Continuation Power Flow.

'''
# Arguments
- 'SolutionVector_CPF': Power Flow Solution Voltage and Angle at each
bus ordered according to bus type: PQ->PV->Lambda.
- 'Tangent_Vector': Tangent Vector for Predict Step of Continuation Power
Flow: PQ->PV->Lambda.
- 'StepSize_CPF': Step size for Continuation Power Flow Predictor Step
'''
'''
# Output
- 'CPF_Predictor_Vector': Continuation Power Flow Solution Vector computed from
the Predict Step
'''
"""
function Compute_PredictVector_CPF(SolutionVector_CPF, Tangent_Vector, StepSize_CPF)

        CPF_Predictor_Vector = SolutionVector_CPF + (StepSize_CPF * Tangent_Vector)

    return CPF_Predictor_Vector

end

"""
    Compute_Corrected_TangentVector_CPF(CDF_DF_List_pu, TangentVector)

Creates corrected tangent vector by converting Deltas from radians to degrees.

'''
# Arguments
- 'CDF_FilePath': File path to the IEEE CDF text file.
- 'Tangent_Vector': Tangent Vector for Predict Step of Continuation Power
Flow: PQ->PV->Lambda.
'''
'''
# Output
- 'Tangent_Vector_Corrected': Corrected tangent vector for deltas from radians to degrees
'''
"""
function Compute_Corrected_TangentVector_CPF(CDF_DF_List_pu, Tangent_Vector)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type_Original == 0) || (row.Type == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type_Original == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type_Original == 3), BusDataCard_DF))

        N_Slack_Bus = nrow(filter(row -> (row.Type_Original == 3), BusDataCard_DF))

        # Initializing Tangent_Vector_Corrected
        Tangent_Vector_Corrected = Tangent_Vector

        # If Else Loop: For NR_Type       

        for ii in 1 : (N_PQ_Bus + N_PV_Bus)

                Tangent_Vector_Corrected[ii] = rad2deg(Tangent_Vector_Corrected[ii])

        end           
        
        return Tangent_Vector_Corrected

end

"""
    Compute_PQ_MismatchVector_CPF(CDF_DF_List_pu, PQ_BusArray, SolutionVector_V,
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
- 'K_Vector': K Vector of CPF.
- 'Lambda': A scalar Lambda value for CPF
'''
'''
# Output
- 'PQ_MismatchVector': An complex array of P-Q Mistmatch elements along with CPF
Corrector additional constraint ordered according to bus type: PQ->PV.
'''
"""
function Compute_PQ_MismatchVector_CPF(CDF_DF_List_pu, PQ_BusArray, SolutionVector_V, NR_Type, K_Vector, Lambda)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type_Original == 0) || (row.Type_Original == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type_Original == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type_Original == 3), BusDataCard_DF))

        # Length of PQ_MismatchVector
        Len_PQ_MismatchVector = 2*N_PQ_Bus+N_PV_Bus

        # Initializing PQ_MismatchVector
        PQ_MismatchVector = Array{Float64}(undef, Len_PQ_MismatchVector, 1)

        # Initializing PQ_PV_Counter
        PQ_PV_Counter = 0

        # Computing P-Q Mismatch
        for ii in 1:Len_PQ_MismatchVector

                if (ii <= (N_Bus-1)) # Excluding Slack Bus

                        if (K_Vector[ii,1] != 0) # Effect of Lambda

                                if ((NR_Type == 1) || (NR_Type == 2)) # Full/Decoupled NR

                                        # Computing P Mismatch
                                        PQ_MismatchVector[ii] = Lambda*(BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii] )- PQ_BusArray[ii+1,1]

                                elseif (NR_Type == 3) # Fast-Decoupled NR

                                        # Computing P Mismatch
                                        PQ_MismatchVector[ii] = (Lambda*(BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii] ))- (PQ_BusArray[ii+1,1])/(SolutionVector_V[ii+1])
                                end

                        else

                                if ((NR_Type == 1) || (NR_Type == 2)) # Full/Decoupled NR

                                        # Computing P Mismatch
                                        PQ_MismatchVector[ii] = (BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii] )- PQ_BusArray[ii+1,1]

                                elseif (NR_Type == 3) # Fast-Decoupled NR

                                        # Computing P Mismatch
                                        PQ_MismatchVector[ii] = ((BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii] )- PQ_BusArray[ii+1,1])/(SolutionVector_V[ii+1])

                                end


                        end


                else

                        # Incrementing PQ_PV_Counter
                        PQ_PV_Counter = PQ_PV_Counter +1

                        if (K_Vector[ii,1] != 0) # Effect of Lambda

                                if ((NR_Type == 1) || (NR_Type == 2)) # Full/Decoupled NR

                                        # Computing Q Mismatch
                                        PQ_MismatchVector[ii] = Lambda*(BusDataCard_DF.Gen_MVAR[PQ_PV_Counter] - BusDataCard_DF.Load_MVAR[PQ_PV_Counter] )- PQ_BusArray[PQ_PV_Counter+1,2]

                                elseif (NR_Type == 3) # Fast-Decoupled NR

                                        # Computing Q Mismatch
                                        PQ_MismatchVector[ii] = (Lambda*(BusDataCard_DF.Gen_MVAR[PQ_PV_Counter] - BusDataCard_DF.Load_MVAR[PQ_PV_Counter] ))- (PQ_BusArray[PQ_PV_Counter+1,2] )/(SolutionVector_V[PQ_PV_Counter+1])

                                end

                        else

                                if ((NR_Type == 1) || (NR_Type == 2)) # Full/Decoupled NR

                                        # Computing Q Mismatch
                                        PQ_MismatchVector[ii] = (BusDataCard_DF.Gen_MVAR[PQ_PV_Counter] - BusDataCard_DF.Load_MVAR[PQ_PV_Counter] )- PQ_BusArray[PQ_PV_Counter+1,2]

                                elseif (NR_Type == 3) # Fast-Decoupled NR

                                        # Computing Q Mismatch
                                        PQ_MismatchVector[ii] = ((BusDataCard_DF.Gen_MVAR[PQ_PV_Counter] - BusDataCard_DF.Load_MVAR[PQ_PV_Counter] )- PQ_BusArray[PQ_PV_Counter+1,2] )/(SolutionVector_V[PQ_PV_Counter+1])

                                end

                        end

                end

        end

        # Adding CPF Corrector Constraint
        PQ_MismatchVector = vcat(PQ_MismatchVector,0)

        return PQ_MismatchVector

end

"""
    Compute_PQ_MismatchVector_Phase_CPF(CDF_DF_List_pu, PQ_BusArray, SolutionVector_V,
    NR_Type, Phase_CPF)

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
- 'K_Vector': K Vector of CPF.
- 'Lambda': A scalar Lambda value for CPF
'''
'''
# Output
- 'PQ_MismatchVector': An complex array of P-Q Mistmatch elements along with CPF
Corrector additional constraint ordered according to bus type: PQ->PV.
'''
"""
function Compute_PQ_MismatchVector_Phase_CPF(CDF_DF_List_pu, PQ_BusArray, SolutionVector_V, NR_Type, K_Vector, Lambda, Phase_CPF)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type_Original == 0) || (row.Type_Original == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type_Original == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type_Original == 3), BusDataCard_DF))

        # Length of PQ_MismatchVector
        Len_PQ_MismatchVector = 2*N_PQ_Bus+N_PV_Bus

        # Initializing PQ_MismatchVector
        PQ_MismatchVector = Array{Float64}(undef, Len_PQ_MismatchVector, 1)

        # Initializing PQ_PV_Counter
        PQ_PV_Counter = 0

        # Computing P-Q Mismatch
        for ii in 1:Len_PQ_MismatchVector

                if (ii <= (N_Bus-1)) # Excluding Slack Bus

                        if (K_Vector[ii,1] != 0) # Effect of Lambda

                                if ((NR_Type == 1) || (NR_Type == 2)) # Full/Decoupled NR

                                        # Computing P Mismatch
                                        PQ_MismatchVector[ii] = Lambda*(BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii] )- PQ_BusArray[ii+1,1]

                                elseif (NR_Type == 3) # Fast-Decoupled NR

                                        # Computing P Mismatch
                                        PQ_MismatchVector[ii] = (Lambda*(BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii] ))- (PQ_BusArray[ii+1,1])/(SolutionVector_V[ii+1])
                                end

                        else

                                if ((NR_Type == 1) || (NR_Type == 2)) # Full/Decoupled NR

                                        # Computing P Mismatch
                                        PQ_MismatchVector[ii] = (BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii] )- PQ_BusArray[ii+1,1]

                                elseif (NR_Type == 3) # Fast-Decoupled NR

                                        # Computing P Mismatch
                                        PQ_MismatchVector[ii] = ((BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii] )- PQ_BusArray[ii+1,1])/(SolutionVector_V[ii+1])

                                end


                        end


                else

                        # Incrementing PQ_PV_Counter
                        PQ_PV_Counter = PQ_PV_Counter +1

                        if (K_Vector[ii,1] != 0) # Effect of Lambda

                                if ((NR_Type == 1) || (NR_Type == 2)) # Full/Decoupled NR

                                        # Computing Q Mismatch
                                        PQ_MismatchVector[ii] = Lambda*(BusDataCard_DF.Gen_MVAR[PQ_PV_Counter] - BusDataCard_DF.Load_MVAR[PQ_PV_Counter] )- PQ_BusArray[PQ_PV_Counter+1,2]

                                elseif (NR_Type == 3) # Fast-Decoupled NR

                                        # Computing Q Mismatch
                                        PQ_MismatchVector[ii] = (Lambda*(BusDataCard_DF.Gen_MVAR[PQ_PV_Counter] - BusDataCard_DF.Load_MVAR[PQ_PV_Counter] ))- (PQ_BusArray[PQ_PV_Counter+1,2] )/(SolutionVector_V[PQ_PV_Counter+1])

                                end

                        else

                                if ((NR_Type == 1) || (NR_Type == 2)) # Full/Decoupled NR

                                        # Computing Q Mismatch
                                        PQ_MismatchVector[ii] = (BusDataCard_DF.Gen_MVAR[PQ_PV_Counter] - BusDataCard_DF.Load_MVAR[PQ_PV_Counter] )- PQ_BusArray[PQ_PV_Counter+1,2]

                                elseif (NR_Type == 3) # Fast-Decoupled NR

                                        # Computing Q Mismatch
                                        PQ_MismatchVector[ii] = ((BusDataCard_DF.Gen_MVAR[PQ_PV_Counter] - BusDataCard_DF.Load_MVAR[PQ_PV_Counter] )- PQ_BusArray[PQ_PV_Counter+1,2] )/(SolutionVector_V[PQ_PV_Counter+1])

                                end

                        end

                end

        end

        # Adding CPF Corrector Constraint while checkin Phase of CPF        
        if ((Phase_CPF) == 1 || (Phase_CPF == 3))

                # Do Nothing

            elseif (Phase_CPF == 2)

                # Add zero to PQ_MismatchVector
                PQ_MismatchVector = vcat(PQ_MismatchVector,0)

            end

        return PQ_MismatchVector

end

"""
    Compute_PQ_MismatchVector_Ini_CPF(CDF_DF_List_pu, PQ_BusArray, SolutionVector_V,
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
- 'K_Vector': K Vector of CPF.
- 'Lambda': A scalar Lambda value for CPF
'''
'''
# Output
- 'PQ_MismatchVector': An complex array of P-Q Mistmatch elements along with CPF
Corrector additional constraint ordered according to bus type: PQ->PV.
'''
"""
function Compute_PQ_MismatchVector_Ini_CPF(CDF_DF_List_pu, PQ_BusArray, SolutionVector_V, NR_Type, K_Vector, Lambda)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type_Original == 0) || (row.Type_Original == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type_Original == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type_Original == 3), BusDataCard_DF))

        # Length of PQ_MismatchVector
        Len_PQ_MismatchVector = 2*N_PQ_Bus+N_PV_Bus

        # Initializing PQ_MismatchVector
        PQ_MismatchVector = Array{Float64}(undef, Len_PQ_MismatchVector, 1)

        # Initializing PQ_PV_Counter
        PQ_PV_Counter = 0

        # Computing P-Q Mismatch
        for ii in 1:Len_PQ_MismatchVector

                if (ii <= (N_Bus-1)) # Excluding Slack Bus

                        if (K_Vector[ii,1] != 0) # Effect of Lambda

                                if ((NR_Type == 1) || (NR_Type == 2)) # Full/Decoupled NR

                                        # Computing P Mismatch
                                        PQ_MismatchVector[ii] = Lambda*(BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii] )- PQ_BusArray[ii+1,1]

                                elseif (NR_Type == 3) # Fast-Decoupled NR

                                        # Computing P Mismatch
                                        PQ_MismatchVector[ii] = (Lambda*(BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii] ))- (PQ_BusArray[ii+1,1])/(SolutionVector_V[ii+1])
                                end

                        else

                                if ((NR_Type == 1) || (NR_Type == 2)) # Full/Decoupled NR

                                        # Computing P Mismatch
                                        PQ_MismatchVector[ii] = (BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii] )- PQ_BusArray[ii+1,1]

                                elseif (NR_Type == 3) # Fast-Decoupled NR

                                        # Computing P Mismatch
                                        PQ_MismatchVector[ii] = ((BusDataCard_DF.Gen_MW[ii] - BusDataCard_DF.Load_MW[ii] )- PQ_BusArray[ii+1,1])/(SolutionVector_V[ii+1])

                                end


                        end


                else

                        # Incrementing PQ_PV_Counter
                        PQ_PV_Counter = PQ_PV_Counter +1

                        if (K_Vector[ii,1] != 0) # Effect of Lambda

                                if ((NR_Type == 1) || (NR_Type == 2)) # Full/Decoupled NR

                                        # Computing Q Mismatch
                                        PQ_MismatchVector[ii] = Lambda*(BusDataCard_DF.Gen_MVAR[PQ_PV_Counter] - BusDataCard_DF.Load_MVAR[PQ_PV_Counter] )- PQ_BusArray[PQ_PV_Counter+1,2]

                                elseif (NR_Type == 3) # Fast-Decoupled NR

                                        # Computing Q Mismatch
                                        PQ_MismatchVector[ii] = (Lambda*(BusDataCard_DF.Gen_MVAR[PQ_PV_Counter] - BusDataCard_DF.Load_MVAR[PQ_PV_Counter] ))- (PQ_BusArray[PQ_PV_Counter+1,2] )/(SolutionVector_V[PQ_PV_Counter+1])

                                end

                        else

                                if ((NR_Type == 1) || (NR_Type == 2)) # Full/Decoupled NR

                                        # Computing Q Mismatch
                                        PQ_MismatchVector[ii] = (BusDataCard_DF.Gen_MVAR[PQ_PV_Counter] - BusDataCard_DF.Load_MVAR[PQ_PV_Counter] )- PQ_BusArray[PQ_PV_Counter+1,2]

                                elseif (NR_Type == 3) # Fast-Decoupled NR

                                        # Computing Q Mismatch
                                        PQ_MismatchVector[ii] = ((BusDataCard_DF.Gen_MVAR[PQ_PV_Counter] - BusDataCard_DF.Load_MVAR[PQ_PV_Counter] )- PQ_BusArray[PQ_PV_Counter+1,2] )/(SolutionVector_V[PQ_PV_Counter+1])

                                end

                        end

                end

        end

        # Adding CPF Corrector Constraint
        #PQ_MismatchVector = vcat(PQ_MismatchVector,0)

        return PQ_MismatchVector

end

"""
    Compute_Corrected_CorrectionVector_CPF(CDF_DF_List_pu, Correction_Vector_NR, SolutionVector_NR, NR_Type)

Creates corrected correction vector for NR Power-Flow based on NR_Type.

'''
# Arguments
- 'CDF_FilePath': File path to the IEEE CDF text file.
- 'Correction_Vector_NR': Correction Vector computed in the NR Power-Flow iteration voltage, angle and 
Lambda at each bus ordered according to bus type: PQ->PV .
- 'SolutionVector_NR': Voltage, Angle and Lambda at each bus ordered according to bus
type: PQ->PV.
'NR_Type': 1 -> Full Newton-Raphson, 2-> Decoupled Newton-Raphson,
3 -> Fast Decoupled Newton-Raphson
'''
'''
# Output
- 'Correction_Vector_NR_Corrected': Corrected correction vector NR Power-Flow based on NR_Type
'''
"""
function Compute_Corrected_CorrectionVector_CPF(CDF_DF_List_pu, Correction_Vector_NR, SolutionVector_NR, NR_Type)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type_Original == 0) || (row.Type == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type_Original == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type_Original == 3), BusDataCard_DF))

        N_Slack_Bus = nrow(filter(row -> (row.Type_Original == 3), BusDataCard_DF))

        # Initializing Correction_Vector_NR_Corrected
        Correction_Vector_NR_Corrected = Correction_Vector_NR

        # If Else Loop: For NR_Type
        if ((NR_Type == 1) || (NR_Type == 2)) # Full NR Decoupled NR

                for ii in 1 : (N_PQ_Bus + N_PV_Bus)

                        Correction_Vector_NR_Corrected[ii] = rad2deg(Correction_Vector_NR_Corrected[ii])

                end

                #= for ii in (N_PQ_Bus + N_PV_Bus + 1) : (length(SolutionVector_NR)-1)

                        Correction_Vector_NR_Corrected[ii] = Correction_Vector_NR_Corrected[ii]*SolutionVector_NR[ii]

                end =#

        elseif (NR_Type == 3) # Fast Decoupled NR

                Correction_Vector_NR_Corrected = Correction_Vector_NR

                for ii in 1 : (N_PQ_Bus + N_PV_Bus)

                        Correction_Vector_NR_Corrected[ii] = rad2deg(Correction_Vector_NR_Corrected[ii])

                end

        end

        return Correction_Vector_NR_Corrected

end

"""
    Compute_ToleranceSatisfaction_CPF(Tolerance, Corrections_Vector)

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
- 'K_Vector': K Vector of CPF.
- 'Lambda': A scalar Lambda value for CPF
'''
'''
# Output
- 'Tolerance_Satisfaction': A boolean value indicating tolerance satifaction
(True) and disatisfaction (False).
'''
"""
function Compute_ToleranceSatisfaction_CPF(CDF_DF_List_pu, Ybus, NR_Type, Tolerance, SolutionVector_NR, K_Vector, Current_Lambda)
        
    # Creating SolutionVector_V and SolutionVector_Delta
    SolutionVector_V, SolutionVector_Delta = Create_SolutionVector_VDelta_CPF(CDF_DF_List_pu, SolutionVector_NR)

    # @show SolutionVector_V
    # @show SolutionVector_Delta

    #  Compute P-Q at Buses
    PQ_BusArray = Compute_PQ_BusArray(Ybus, SolutionVector_V, SolutionVector_Delta)

    #  Compute PQ Mistmatch (DEBUG add CPF version of Compute PQ mismatch vector)
    PQ_MismatchVector = Compute_PQ_MismatchVector_CPF(CDF_DF_List_pu, PQ_BusArray, SolutionVector_V, NR_Type,  K_Vector, Current_Lambda)
    
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
    Compute_ToleranceSatisfaction_Ini_CPF(Tolerance, Corrections_Vector)

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
- 'K_Vector': K Vector of CPF.
- 'Lambda': A scalar Lambda value for CPF
'''
'''
# Output
- 'Tolerance_Satisfaction': A boolean value indicating tolerance satifaction
(True) and disatisfaction (False).
'''
"""
function Compute_ToleranceSatisfaction_Ini_CPF(CDF_DF_List_pu, Ybus, NR_Type, Tolerance, SolutionVector_NR, K_Vector, Current_Lambda)
        
    # Creating SolutionVector_V and SolutionVector_Delta
    SolutionVector_V, SolutionVector_Delta = Create_SolutionVector_VDelta_NR(CDF_DF_List_pu, SolutionVector_NR)

    # @show SolutionVector_V
    # @show SolutionVector_Delta

    #  Compute P-Q at Buses
    PQ_BusArray = Compute_PQ_BusArray(Ybus, SolutionVector_V, SolutionVector_Delta)

    #  Compute PQ Mistmatch (DEBUG add CPF version of Compute PQ mismatch vector)
    PQ_MismatchVector = Compute_PQ_MismatchVector_CPF(CDF_DF_List_pu, PQ_BusArray, SolutionVector_V, NR_Type,  K_Vector, Current_Lambda)
    
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
    PowerFlow_MainFunction_CPF(CDF_DF_List_pu, Ybus, NR_Type,
    Initial_SolutionVector_CPF, Tolerance, Tol_Num, K_Vector, Index_CPF)

Creates Ybus without taps for a power system network.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file coverted to pu: [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
- 'Ybus': Ybus of the Power Network.
- 'NR_Type': 1 -> Full Newton-Raphson, 2-> Decoupled Newton-Raphson,
3 -> Fast Decoupled Newton-Raphson.
- 'Initial_SolutionVector_CPF': Solution vector from Predictor of CPF.
- 'K_Vector': K Vector of CPF.
- 'Index_CPF': Index of the Continuation Parameter of CPF.
- 'Tolerance': Tolerance level for stopping criterion of Newton-Raphson Method.
- 'Tol_Num': Tolerance for being near zero.
'''
'''
# Output
- 'SolutionVector_NR': Final Solution of the NR method for CPF Corrector Step.
- 'PowerFlow_IterationTimeInfo_Array': An array containing Iteration Number and
Time in seconds for each iteration.
'''
"""
function PowerFlow_MainFunction_CPF(CDF_DF_List_pu, Ybus, NR_Type, Initial_SolutionVector_CPF, Tolerance, Tol_Num, K_Vector, Index_CPF)

    # Getting required data from CDF_DF_List
    BusDataCard_DF = CDF_DF_List_pu[2]

    # Number of Buses
    N_Bus = nrow(BusDataCard_DF)
    N_PQ_Bus = nrow(filter(row -> ((row.Type == 0) || (row.Type == 1)), BusDataCard_DF))
    N_PV_Bus = nrow(filter(row -> (row.Type == 2), BusDataCard_DF))
    N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

    N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

    # We have Y_Bus

    # We have Initial_SolutionVector_NR

    # Initializing Tolerance_Satisfaction
    Tolerance_Satisfaction = false

    # Initializing PowerFlow_IterationTimeInfo_Array
    PowerFlow_IterationTimeInfo_Array = zeros(1,2)

    # Intializing WhileLoop_Counter
    WhileLoop_Counter = 0

    # Initializations
    SolutionVector_NR = 0

    #Current_Lambda = SolutionVector_NR[end,1]

    # While Loop: For Newton-Raphson Method
    while ((!Tolerance_Satisfaction))

            # Incrementing WhileLoop_Counter
            WhileLoop_Counter = WhileLoop_Counter + 1

            # Starting Timer
            tick()

            # If Else Loop: For checking First Iteration
            if (WhileLoop_Counter == 1)

                    # Compute SolutionVector_V, SolutionVector_Delta
                    SolutionVector_V, SolutionVector_Delta = Create_SolutionVector_VDelta_CPF(CDF_DF_List_pu, Initial_SolutionVector_CPF)

                    # Initializing Solution Vector NR
                    SolutionVector_NR = Initial_SolutionVector_CPF

                    # Getting Current Lambda
                    Current_Lambda = SolutionVector_NR[end,1]

            else

                    # Compute SolutionVector_V, SolutionVector_Delta
                    SolutionVector_V, SolutionVector_Delta = Create_SolutionVector_VDelta_CPF(CDF_DF_List_pu, SolutionVector_NR)

                    # Increasing  PowerFlow_IterationTimeInfo_Array Size
                    PowerFlow_IterationTimeInfo_Array = vcat(PowerFlow_IterationTimeInfo_Array, zeros(1,2))

                    # Getting Current Lambda
                    Current_Lambda = SolutionVector_NR[end,1]

            end

            # Compute PQ_BusArray
            PQ_BusArray = Compute_PQ_BusArray(Ybus, SolutionVector_V, SolutionVector_Delta)

            # Compute PQ_MismatchVector
            PQ_MismatchVector = Compute_PQ_MismatchVector_CPF(CDF_DF_List_pu, PQ_BusArray, SolutionVector_V, NR_Type, K_Vector, Current_Lambda)

            # We have Solution Vector NR

            # Compute Continuation Power Flow Jacobian Corrector Step
            Jacobian_CPF_Correct = Create_Jacobian_CPF_Correct(CDF_DF_List_pu, Ybus, SolutionVector_NR, 1, K_Vector, Index_CPF)

            # Compute Correction_Vector_NR using LU Factorization
            Correction_Vector_NR = PLU_Solve(Jacobian_CPF_Correct, PQ_MismatchVector, Tol_Num)

            # Correcting Correction_Vector_NR
            Correction_Vector_NR = Compute_Corrected_CorrectionVector_CPF(CDF_DF_List_pu, Correction_Vector_NR, SolutionVector_NR, NR_Type)

            # Computing New SolutionVector_NR
            SolutionVector_NR  = SolutionVector_NR + Correction_Vector_NR

            # Compute Tolerance Satisfaction
            Tolerance_Satisfaction = Compute_ToleranceSatisfaction_CPF(CDF_DF_List_pu, Ybus, NR_Type, Tolerance, SolutionVector_NR,  K_Vector, Current_Lambda)
           
            # Stopping Timer
            IterationTime = tok()

            # Filling-up PowerFlow_IterationTimeInfo_Array
            PowerFlow_IterationTimeInfo_Array[WhileLoop_Counter,1:2] = [WhileLoop_Counter , IterationTime]

    end


    return SolutionVector_NR, PowerFlow_IterationTimeInfo_Array

end

"""
    PowerFlow_MainFunction_Phase_CPF(CDF_DF_List_pu, Ybus, NR_Type,
    Initial_SolutionVector_CPF, Tolerance, Tol_Num, K_Vector, Index_CPF, Phase_CPF)

Creates Ybus without taps for a power system network.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file coverted to pu: [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
- 'Ybus': Ybus of the Power Network.
- 'NR_Type': 1 -> Full Newton-Raphson, 2-> Decoupled Newton-Raphson,
3 -> Fast Decoupled Newton-Raphson.
- 'Initial_SolutionVector_CPF': Solution vector from Predictor of CPF.
- 'K_Vector': K Vector of CPF.
- 'Index_CPF': Index of the Continuation Parameter of CPF.
- 'Tolerance': Tolerance level for stopping criterion of Newton-Raphson Method.
- 'Tol_Num': Tolerance for being near zero.
'''
'''
# Output
- 'SolutionVector_NR': Final Solution of the NR method for CPF Corrector Step.
- 'PowerFlow_IterationTimeInfo_Array': An array containing Iteration Number and
Time in seconds for each iteration.
'''
"""
function PowerFlow_MainFunction_Phase_CPF(CDF_DF_List_pu, Ybus, NR_Type, Initial_SolutionVector_CPF, Tolerance, Tol_Num, K_Vector, Index_CPF, Phase_CPF)

    # Getting required data from CDF_DF_List
    BusDataCard_DF = CDF_DF_List_pu[2]

    # Number of Buses
    N_Bus = nrow(BusDataCard_DF)
    N_PQ_Bus = nrow(filter(row -> ((row.Type == 0) || (row.Type == 1)), BusDataCard_DF))
    N_PV_Bus = nrow(filter(row -> (row.Type == 2), BusDataCard_DF))
    N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

    N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

    # We have Y_Bus

    # We have Initial_SolutionVector_NR

    # Initializing Tolerance_Satisfaction
    Tolerance_Satisfaction = false

    # Initializing PowerFlow_IterationTimeInfo_Array
    PowerFlow_IterationTimeInfo_Array = zeros(1,2)

    # Intializing WhileLoop_Counter
    WhileLoop_Counter = 0

    # Initializations
    SolutionVector_NR = 0

    #Current_Lambda = SolutionVector_NR[end,1]

    # While Loop: For Newton-Raphson Method
    while ((!Tolerance_Satisfaction) && (WhileLoop_Counter < 25))

            # Incrementing WhileLoop_Counter
            WhileLoop_Counter = WhileLoop_Counter + 1

            # Starting Timer
            tick()

            # If Else Loop: For checking First Iteration
            if (WhileLoop_Counter == 1)

                    # Compute SolutionVector_V, SolutionVector_Delta
                    SolutionVector_V, SolutionVector_Delta = Create_SolutionVector_VDelta_CPF(CDF_DF_List_pu, Initial_SolutionVector_CPF)

                    # Initializing Solution Vector NR
                    SolutionVector_NR = Initial_SolutionVector_CPF

                    # Getting Current Lambda
                    Current_Lambda = SolutionVector_NR[end,1]

            else

                    # Compute SolutionVector_V, SolutionVector_Delta
                    SolutionVector_V, SolutionVector_Delta = Create_SolutionVector_VDelta_CPF(CDF_DF_List_pu, SolutionVector_NR)

                    # Increasing  PowerFlow_IterationTimeInfo_Array Size
                    PowerFlow_IterationTimeInfo_Array = vcat(PowerFlow_IterationTimeInfo_Array, zeros(1,2))

                    # Getting Current Lambda
                    Current_Lambda = SolutionVector_NR[end,1]

            end

            # Compute PQ_BusArray
            PQ_BusArray = Compute_PQ_BusArray(Ybus, SolutionVector_V, SolutionVector_Delta)

            # Compute PQ_MismatchVector
            PQ_MismatchVector = Compute_PQ_MismatchVector_Phase_CPF(CDF_DF_List_pu, PQ_BusArray, SolutionVector_V, NR_Type, K_Vector, Current_Lambda, Phase_CPF)

            # We have Solution Vector NR

            # Compute Continuation Power Flow Jacobian Corrector Step
            Jacobian_CPF_Correct = Create_Jacobian_Phase_CPF_Correct(CDF_DF_List_pu, Ybus, SolutionVector_NR, 1, K_Vector, Index_CPF, Phase_CPF)

            # Compute Correction_Vector_NR using LU Factorization
            Correction_Vector_NR = PLU_Solve(Jacobian_CPF_Correct, PQ_MismatchVector, Tol_Num)

            # Correcting Correction_Vector_NR
            Correction_Vector_NR = Compute_Corrected_CorrectionVector_CPF(CDF_DF_List_pu, Correction_Vector_NR, SolutionVector_NR, NR_Type)

            # Checking Phase of CPF
            if ((Phase_CPF) == 1 || (Phase_CPF == 3))

                # Add zero change to Lambda
                Correction_Vector_NR = vcat(Correction_Vector_NR,0)

            elseif (Phase_CPF == 2)

                # Do Nothing

            end

            # Computing New SolutionVector_NR
            SolutionVector_NR  = SolutionVector_NR + Correction_Vector_NR

            # Compute Tolerance Satisfaction
            Tolerance_Satisfaction = Compute_ToleranceSatisfaction_CPF(CDF_DF_List_pu, Ybus, NR_Type, Tolerance, SolutionVector_NR,  K_Vector, Current_Lambda)
           
            # Stopping Timer
            IterationTime = tok()

            # Filling-up PowerFlow_IterationTimeInfo_Array
            PowerFlow_IterationTimeInfo_Array[WhileLoop_Counter,1:2] = [WhileLoop_Counter , IterationTime]

    end

    if (WhileLoop_Counter >= 25)

        Corrector_Convergence_Indicator = false

    else

        Corrector_Convergence_Indicator = true

    end


    return SolutionVector_NR, PowerFlow_IterationTimeInfo_Array, Corrector_Convergence_Indicator

end

"""
    PowerFlow_MainFunction_Ini_CPF(CDF_DF_List_pu, Ybus, NR_Type,
    Initial_SolutionVector_CPF, Tolerance, Tol_Num, K_Vector, Index_CPF)

Creates Ybus without taps for a power system network.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file coverted to pu: [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
- 'Ybus': Ybus of the Power Network.
- 'NR_Type': 1 -> Full Newton-Raphson, 2-> Decoupled Newton-Raphson,
3 -> Fast Decoupled Newton-Raphson.
- 'Initial_SolutionVector_CPF': Solution vector from Predictor of CPF.
- 'K_Vector': K Vector of CPF.
- 'Index_CPF': Index of the Continuation Parameter of CPF.
- 'Tolerance': Tolerance level for stopping criterion of Newton-Raphson Method.
- 'Tol_Num': Tolerance for being near zero.
'''
'''
# Output
- 'SolutionVector_NR': Final Solution of the NR method for CPF Corrector Step.
- 'PowerFlow_IterationTimeInfo_Array': An array containing Iteration Number and
Time in seconds for each iteration.
'''
"""
function PowerFlow_MainFunction_Ini_CPF(CDF_DF_List_pu, Ybus, NR_Type, Tolerance, Tol_Num, K_Vector, Lambda_Ini)

    # Getting required data from CDF_DF_List
    BusDataCard_DF = CDF_DF_List_pu[2]

    # Number of Buses
    N_Bus = nrow(BusDataCard_DF)
    N_PQ_Bus = nrow(filter(row -> ((row.Type == 0) || (row.Type == 1)), BusDataCard_DF))
    N_PV_Bus = nrow(filter(row -> (row.Type == 2), BusDataCard_DF))
    N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

    N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

    # We have Y_Bus

    # We have Initial_SolutionVector_NR
    Initial_SolutionVector_CPF = Create_Initial_SolutionVector_NR(CDF_DF_List_pu)
    Current_Lambda =  Lambda_Ini

    # Initializing Tolerance_Satisfaction
    Tolerance_Satisfaction = false

    # Initializing PowerFlow_IterationTimeInfo_Array
    PowerFlow_IterationTimeInfo_Array = zeros(1,2)

    # Intializing WhileLoop_Counter
    WhileLoop_Counter = 0

    # Initializations
    SolutionVector_NR = 0

    #Current_Lambda = SolutionVector_NR[end,1]

    # While Loop: For Newton-Raphson Method
    while (!Tolerance_Satisfaction)

            # Incrementing WhileLoop_Counter
            WhileLoop_Counter = WhileLoop_Counter + 1
            @show WhileLoop_Counter

            # Starting Timer
            tick()

            # If Else Loop: For checking First Iteration
            if (WhileLoop_Counter == 1)

                    # Compute SolutionVector_V, SolutionVector_Delta
                    SolutionVector_V, SolutionVector_Delta = Create_SolutionVector_VDelta_NR(CDF_DF_List_pu, Initial_SolutionVector_CPF)

                    # Initializing Solution Vector NR
                    SolutionVector_NR = Initial_SolutionVector_CPF

                    # Getting Current Lambda
                    #Current_Lambda = SolutionVector_NR[end,1]

            else

                    # Compute SolutionVector_V, SolutionVector_Delta
                    SolutionVector_V, SolutionVector_Delta = Create_SolutionVector_VDelta_NR(CDF_DF_List_pu, SolutionVector_NR)

                    # Increasing  PowerFlow_IterationTimeInfo_Array Size
                    PowerFlow_IterationTimeInfo_Array = vcat(PowerFlow_IterationTimeInfo_Array, zeros(1,2))

                    # Getting Current Lambda
                    #Current_Lambda = SolutionVector_NR[end,1]

            end

            # Compute PQ_BusArray
            PQ_BusArray = Compute_PQ_BusArray(Ybus, SolutionVector_V, SolutionVector_Delta)

            # Compute PQ_MismatchVector
            PQ_MismatchVector = Compute_PQ_MismatchVector_Ini_CPF(CDF_DF_List_pu, PQ_BusArray, SolutionVector_V, NR_Type, K_Vector, Current_Lambda)

            # We have Solution Vector NR

            # Compute Continuation Power Flow Jacobian Corrector Step
            Jacobian_CPF_Correct = Create_Jacobian_NR(CDF_DF_List_pu, Ybus, SolutionVector_V, SolutionVector_Delta, PQ_BusArray, NR_Type, 1)

            # Compute Correction_Vector_NR using LU Factorization
            Correction_Vector_NR = PLU_Solve(Jacobian_CPF_Correct, PQ_MismatchVector, Tol_Num)

            # Correcting Correction_Vector_NR
            Correction_Vector_NR = Compute_Corrected_CorrectionVector(CDF_DF_List_pu, Correction_Vector_NR, SolutionVector_NR, NR_Type)

            # Computing New SolutionVector_NR
            SolutionVector_NR = SolutionVector_NR + Correction_Vector_NR

            # Compute Tolerance Satisfaction
            Tolerance_Satisfaction = Compute_ToleranceSatisfaction_Ini_CPF(CDF_DF_List_pu, Ybus, NR_Type, Tolerance, SolutionVector_NR,  K_Vector, Current_Lambda)
           
            # Stopping Timer
            IterationTime = tok()

            # Filling-up PowerFlow_IterationTimeInfo_Array
            PowerFlow_IterationTimeInfo_Array[WhileLoop_Counter,1:2] = [WhileLoop_Counter , IterationTime]

    end


    return SolutionVector_NR, PowerFlow_IterationTimeInfo_Array

end
