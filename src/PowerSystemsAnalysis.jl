module PowerSystemsAnalysis

# Importing External Modules
using DataFrames
using LinearAlgebra
using TickTock

# Export Functions from Included Files
export CDF_Parser,
       CDF_pu_Converter,
       Create_Ybus_WithoutTaps,
       Create_Ybus_WithTaps,
       Create_Initial_SolutionVector_NR,
       Create_SolutionVector_VDelta_NR,
       Create_SolutionVector_NR,
       Compute_PQ_BusArray,
       Compute_PQ_MismatchVector,
       PQ_PV_Bus_Check_Modify,
       Compute_LineFlows,
       Compute_ToleranceSatisfaction,
       Compute_Corrected_CorrectionVector,
       Create_Jacobian_NR,
       Create_Jacobian_CPF_Predict,
       Create_Jacobian_CPF_Correct,
       Compute_PLU,
       LU_ForwardSubstitution,
       LU_BackwardSubstitution,
       PLU_Solve,
       PowerFlow_MainFunction,
       Create_Initial_SolutionVector_CPF,
       Create_SolutionVector_VDelta_CPF,
       Create_SolutionVector_CPF,
       Create_KVector_CPF,
       Compute_Tangent_Vector_CPF,
       Check_CriticalPoint_CPF,
       Choose_ContinuationParameter_CPF,
       Compute_PredictVector_CPF,
       Compute_PQ_MismatchVector_CPF,
       Compute_Corrected_CorrectionVector_CPF,
       PowerFlow_MainFunction_CPF,
       CDF_AddMeasurements_SE




# Main Functions

"""
    PowerFlow_MainFunction(CDF_FilePath, Ybus_Taps_Indicator, NR_Type, Tolerance, Tol_Num))

Creates Ybus without taps for a power system network.

'''
# Arguments
- 'CDF_FilePath': File path to the IEEE CDF text file.
- 'Ybus_Taps_Indicator': 1 - Ybus with no taps, 2 - Ybus with Taps
- 'NR_Type': 1 -> Full Newton-Raphson, 2-> Decoupled Newton-Raphson,
3 -> Fast Decoupled Newton-Raphson
- 'Tolerance': Tolerance level for stopping criterion of Newton-Raphson Method.
- 'Tol_Num': Tolerance for being near zero
'''
'''
# Output
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF]. Changed with changes in Line flows
- 'LineFlow_Array': An array (N_Lines*2) for P and Q line flows ordered
according to Branch Data Card.
- 'PowerFlow_IterationTimeInfo_Array': An array containing Iteration Number and
Time in seconds for each iteration.
'''
"""
function PowerFlow_MainFunction(CDF_FilePath, Ybus_Taps_Indicator, NR_Type, Tolerance, Tol_Num)

    # Reading IEEE CDF File
    CDF_DF_List = CDF_Parser(CDF_FilePath)

    # Converting CDF DataFrame to PU
    CDF_DF_List_pu = CDF_pu_Converter(CDF_DF_List)

    # Getting required data from CDF_DF_List
    BusDataCard_DF = CDF_DF_List_pu[2]

    # Number of Buses
    N_Bus = nrow(BusDataCard_DF)
    N_PQ_Bus = nrow(filter(row -> ((row.Type == 0) || (row.Type == 1)), BusDataCard_DF))
    N_PV_Bus = nrow(filter(row -> (row.Type == 2), BusDataCard_DF))
    N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

    N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

    # Create Ybus
    if (Ybus_Taps_Indicator == 1) # Without Taps

            Ybus = Create_Ybus_WithoutTaps(CDF_DF_List_pu)

    elseif (Ybus_Taps_Indicator == 2) # With Taps

            Ybus_WithoutTaps = Create_Ybus_WithoutTaps(CDF_DF_List_pu)

            Ybus = Create_Ybus_WithTaps(Ybus_WithoutTaps,CDF_DF_List_pu)

    end

    # Creating Initial Solution Vector for Power Flow
    Initial_SolutionVector_NR = Create_Initial_SolutionVector_NR(CDF_DF_List_pu)

    # Initializing Tolerance_Satisfaction
    Tolerance_Satisfaction = false

    # Initializing PowerFlow_IterationTimeInfo_Array
    PowerFlow_IterationTimeInfo_Array = zeros(1,2)

    # Intializing WhileLoop_Counter
    WhileLoop_Counter = 0

    # While Loop: For Newton-Raphson Method
    while (!Tolerance_Satisfaction)

            # Incrementing WhileLoop_Counter
            WhileLoop_Counter = WhileLoop_Counter + 1

            # Starting Timer
            tick()

            # If Else Loop: For checking First Iteration
            if (WhileLoop_Counter == 1)

                    # Compute Mismatch
                    CDF_DF_List_pu, Ybus, PQ_BusArray, PQ_MismatchVector, SolutionVector_V, SolutionVector_Delta = PQ_PV_Bus_Check_Modify(CDF_DF_List_pu, Ybus, Ybus_Type, NR_Type, Initial_SolutionVector_NR)

            else

                    # Compute Mismatch
                    CDF_DF_List_pu, Ybus, PQ_BusArray, PQ_MismatchVector, SolutionVector_V, SolutionVector_Delta = PQ_PV_Bus_Check_Modify(CDF_DF_List_pu, Ybus, Ybus_Type, NR_Type, SolutionVector_NR)

                    # Increasing  PowerFlow_IterationTimeInfo_Array Size
                    PowerFlow_IterationTimeInfo_Array = vcat(PowerFlow_IterationTimeInfo_Array, zeros(1,2))

            end

            # Create Solution Vector NR
            SolutionVector_NR = Create_SolutionVector_NR(CDF_DF_List_pu, SolutionVector_V, SolutionVector_Delta)

            # Compute Power Flow Jacobian
            Jacobian_NR = Create_Jacobian_NR(CDF_DF_List_pu, Ybus, SolutionVector_V, SolutionVector_Delta, PQ_BusArray, NR_Type, 0)

            # Compute Correction_Vector_NR using LU Factorization
            if (NR_Type == 1) # Full NR

                    Correction_Vector_NR = PLU_Solve(Jacobian_NR, PQ_MismatchVector, Tol_Num)

                    # Correcting Correction_Vector_NR
                    Correction_Vector_NR = Compute_Corrected_CorrectionVector(CDF_DF_List_pu, Correction_Vector_NR, SolutionVector_NR, NR_Type)

            elseif ((NR_Type == 2) || (NR_Type == 3)) # Decoupled/Fast Decoupled NR

                    Correction_Vector_NR_Delta = PLU_Solve(Jacobian_NR[1], PQ_MismatchVector[1:(N_PQ_Bus+N_PV_Bus),1], Tol_Num)

                    Correction_Vector_NR_V = PLU_Solve(Jacobian_NR[2], PQ_MismatchVector[(N_PQ_Bus+N_PV_Bus+1):end,1], Tol_Num)

                    # Creating Correction_Vector_NR
                    Correction_Vector_NR = vcat(Correction_Vector_NR_Delta, Correction_Vector_NR_V)

                    # Correcting Corrections Vector
                    if (NR_Type == 2) # Decoupled NR

                            Correction_Vector_NR = Compute_Corrected_CorrectionVector(CDF_DF_List_pu, Correction_Vector_NR, SolutionVector_NR, NR_Type)

                    elseif (NR_Type == 3) # Fast Decoupled NR

                            Correction_Vector_NR =  Correction_Vector_NR

                    end

            end

            # Compute Tolerance Satisfaction
            Tolerance_Satisfaction = Compute_ToleranceSatisfaction(Tolerance, Correction_Vector_NR)

            # Computing New SolutionVector_NR
            SolutionVector_NR = SolutionVector_NR + Correction_Vector_NR

            # Stopping Timer
            IterationTime = tok()

            # Filling-up PowerFlow_IterationTimeInfo_Array
            PowerFlow_IterationTimeInfo_Array[WhileLoop_Counter,1:2] = [WhileLoop_Counter , IterationTime]

    end


    # Compute Line Flows
    CDF_DF_List_pu, LineFlow_Array = Compute_LineFlows(CDF_DF_List_pu, Ybus, SolutionVector_V, SolutionVector_Delta)

    # ReOrdering BusDataCard_DF: PQ->PV->Slack
    sort!(BusDataCard_DF, [order(:Type_Original)])

    return CDF_DF_List_pu, LineFlow_Array, PowerFlow_IterationTimeInfo_Array

end

"""
    ContinuationPowerFlow_MainFunction(CDF_FilePath, Ybus_Taps_Indicator, NR_Type, Tol_Num))

Creates Ybus without taps for a power system network.

'''
# Arguments
- 'CDF_FilePath': File path to the IEEE CDF text file.
- 'PQ_V_Curve_Tuple': A Tuple containing (BusNumber,P/Q to be parametrized)
- 'Ybus_Taps_Indicator': 1 - Ybus with no taps, 2 - Ybus with Taps
- 'StepSize_CPF': Step size for Continuation Power Flow Predictor Step
- 'Tolerance': Tolerance level for stopping criterion of Newton-Raphson Method.
- 'Tol_Num': Tolerance for being near zero
- 'PostCriticalPoint_Counter_Input': Number of iterations to be computed after
finding critical point.
'''
'''
# Output
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF]. Changed with changes in Line flows
- 'LineFlow_Array': An array (N_Lines*2) for P and Q line flows ordered
according to Branch Data Card.
- 'PowerFlow_IterationTimeInfo_Array': An array containing Iteration Number and
Time in seconds for each iteration.
'''
"""
function ContinuationPowerFlow_MainFunction(CDF_FilePath, PQ_V_Curve_Tuple, Ybus_Taps_Indicator, StepSize_CPF, Tolerance_NR, Tol_Num, PostCriticalPoint_Counter_Input)

    # Solving Initial Power Flow Problem
   CDF_DF_List_pu, LineFlow_Array, PowerFlow_IterationTimeInfo_Array =  PowerFlow_MainFunction(CDF_FilePath, Ybus_Taps_Indicator, 1, Tolerance_NR, Tol_Num)

   # Create Ybus
   if (Ybus_Taps_Indicator == 1) # Without Taps

           Ybus = Create_Ybus_WithoutTaps(CDF_DF_List_pu)

   elseif (Ybus_Taps_Indicator == 2) # With Taps

           Ybus_WithoutTaps = Create_Ybus_WithoutTaps(CDF_DF_List_pu)

           Ybus = Create_Ybus_WithTaps(Ybus_WithoutTaps,CDF_DF_List_pu)

   end

   # Creating Initial Solution Vector
   Initial_SolutionVector_CPF = Create_Initial_SolutionVector_CPF(CDF_DF_List_pu)

   # Creating K_Vector
   K_Vector = Create_KVector_CPF(CDF_DF_List_pu, PQ_V_Curve_Tuple)

   # Initializing Predictor/Corrector/Tangent Vector History
   Predictor_Vector_History = zeros(size(K_Vector)[1]+1,1)
   Corrector_Vector_History = zeros(size(K_Vector)[1]+1,1)
   Tangent_Vector_History = zeros(size(K_Vector)[1]+1,1)

   # Initializing ContinuationPowerFlow_IterationTimeInfo_Array
   PowerFlow_IterationTimeInfo_Array = zeros(1,2)

   # Initializing PostCriticalPoint_Counter
   PostCriticalPoint_Counter = 0

   # Intializing WhileLoop_Counter
   WhileLoop_Counter = 0

   # While Loop: For computing Continuation Power Flow
   while (PostCriticalPoint_Counter <= PostCriticalPoint_Counter_Input)

           # Incrementing WhileLoop_Counter
           WhileLoop_Counter = WhileLoop_Counter + 1

           # Starting Timer
           tick()

           # If Else Loop: For checking first iteration
           if (WhileLoop_Counter == 1)

                   # Getting Solution Vector comes from Load Flow Solution
                   SolutionVector_CPF = Initial_SolutionVector_CPF

                   # Getting Index k (Choosing Lambda)
                   Index_CPF = size(K_Vector)[1]+1

                   # Initializing Predict_b_Value_CPF
                   Predict_b_Value_CPF = 1

           else
                   # Getting Solution Vector From Corrector
                   SolutionVector_CPF = Corrector_Vector_History[:,WhileLoop_Counter-1]

                   # Getting Index k from previous iteration
                   Index_CPF = Index_CPF_New

           end

           # Compute Continuation Power Flow Jacobian Predictor Step
           Jacobian_CPF_Predict = Create_Jacobian_CPF_Predict(CDF_DF_List_pu, Ybus, SolutionVector_CPF, 1, K_Vector, Index_CPF)

           # Compute Tangent Vector
           Tangent_Vector = Compute_Tangent_Vector_CPF(Jacobian_CPF_Predict, Tol_Num, Predict_b_Value_CPF, Index_CPF)

           # Check if Critical Point is passed
           PostCriticalPoint_Counter = Check_CriticalPoint_CPF(Tangent_Vector, PostCriticalPoint_Counter)

           # Choose Continuation Parameter
           Index_CPF_New, Predict_b_Value_CPF = Choose_ContinuationParameter_CPF(Tangent_Vector)

           # Predict Solution
           CPF_Predictor_Vector = Compute_PredictVector_CPF(SolutionVector_CPF, Tangent_Vector, StepSize_CPF)

           # Compute Continuation Power Flow Jacobian Corrector Step
           Jacobian_CPF_Correct = Create_Jacobian_CPF_Correct(CDF_DF_List_pu, Ybus, CPF_Predictor_Vector, 1, K_Vector, Index_CPF)

           # Update/Correct Solution
           CPF_Corrector_Vector, PowerFlow_IterationTimeInfo_Array_Corrector = PowerFlow_MainFunction_CPF(CDF_DF_List_pu, Ybus, NR_Type, Initial_SolutionVector_CPF, Tolerance, Tol_Num, K_Vector, Index_CPF)

           # Update Predictor/Corrector Vector History
           Predictor_Vector_History = hcat(Predictor_Vector_History, CPF_Predictor_Vector)
           Corrector_Vector_History = hcat(Corrector_Vector_History, CPF_Corrector_Vector)
           Tangent_Vector_History = hcat(Tangent_Vector_History, Tangent_Vector)

           # Stopping Timer
           IterationTime = tok()

           # Filling-up PowerFlow_IterationTimeInfo_Array
           PowerFlow_IterationTimeInfo_Array[WhileLoop_Counter,1:2] = [WhileLoop_Counter , IterationTime]

   end

   ## Add Plotting Code ##

   return Predictor_Vector_History, Corrector_Vector_History, Tangent_Vector_History, PowerFlow_IterationTimeInfo_Array

end

"""
    PowerSystem_StateEstimation_MainFunction(CDF_FilePath, Ybus_Taps_Indicator, NR_Type, Tolerance, Tol_Num Tol_Num))

Creates Ybus without taps for a power system network.

'''
# Arguments
- 'CDF_FilePath': File path to the IEEE CDF text file.
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
- 'Ybus_Taps_Indicator': 1 - Ybus with no taps, 2 - Ybus with Taps
- 'StepSize_CPF': Step size for Continuation Power Flow Predictor Step
- 'Tolerance': Tolerance level for stopping criterion of Newton-Raphson Method.
- 'Tol_Num': Tolerance for being near zero
- 'PostCriticalPoint_Counter_Input': Number of iterations to be computed after
finding critical point.
'''
'''
# Output
- '':
'''
"""
function PowerSystem_StateEstimation_MainFunction(CDF_FilePath, Bus_Measurement_Array, Branch_Measurement_Array, Ybus_Taps_Indicator, StepSize_CPF, Tolerance_NR, Tol_Num, PostCriticalPoint_Counter_Input)

        # Reading IEEE CDF File
        CDF_DF_List = CDF_Parser(CDF_FilePath)

        # Add Bus_Measurement_Array, Branch_Measurement_Array to CDF_DF_List
        CDF_DF_List = CDF_AddMeasurements_SE(CDF_DF_List, Bus_Measurement_Array, Branch_Measurement_Array)

        # Converting CDF DataFrame to PU
        CDF_DF_List_pu = CDF_pu_Converter_SE(CDF_DF_List)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]
        BranchDataCard_DF = CDF_DF_List[3]

        # Number of Buses and Branches
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type == 0) || (row.Type == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

        N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

        N_Branch = nrow(BranchDataCard_DF)

        # Create Ybus
        if (Ybus_Taps_Indicator == 1) # Without Taps

                Ybus = Create_Ybus_WithoutTaps(CDF_DF_List_pu)

        elseif (Ybus_Taps_Indicator == 2) # With Taps

                Ybus_WithoutTaps = Create_Ybus_WithoutTaps(CDF_DF_List_pu)

                Ybus = Create_Ybus_WithTaps(Ybus_WithoutTaps,CDF_DF_List_pu)

        end

end

# Including component files
include("IEEE_CDF_Parser.jl")

include("Ybus_Builder.jl")

include("BasicPowerFlow_Functions.jl")

include("Jacobian_Builder.jl")

include("LU_Factorization.jl")

include("BasicContinuationPowerFlow_Functions.jl")

include("BasicStateEstimation_Functions.jl")

end
