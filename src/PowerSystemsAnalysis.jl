module PowerSystemsAnalysis

# Importing External Modules
using DataFrames
using LinearAlgebra
using Random
using Distributions
using TickTock
using CSV
using Plots
using LaTeXStrings

# Export Functions from Included Files
export CDF_Parser,
       CDF_pu_Converter,
       Create_Ybus_WithoutTaps,
       Create_Ybus_WithTaps,
       initializeVectors_pu,
       sortMatrixByBusTypes,
       extractSystemName,
       createFolderIfNotExisting,
       ybusGenerator ,
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
       Compute_Corrected_Matrix,
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
    PowerFlow_MainFunction(CDF_FilePath, Ybus_Taps_Indicator, NR_Type, Tolerance, Tol_Num)

Creates Ybus without taps for a power system network.

'''
# Arguments
- 'CDF_FilePath': File path to the IEEE CDF text file.
- 'Ybus_Taps_Indicator': False - Ybus with no taps, True - Ybus with Taps
- 'NR_Type': 1 -> Full Newton-Raphson, 2-> Decoupled Newton-Raphson,
3 -> Fast Decoupled Newton-Raphson
- 'Tolerance': Tolerance level for stopping criterion of Newton-Raphson Method.
- 'Tol_Num': Tolerance for being near zero
- 'SortValue': True -> Sort CDF File according to Bus Type PQ->PV->Slack, False -> do not sort
- 'BusSwitching': True -> Bus Switching Employed, False -> Bus Switching not Employed
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
function PowerFlow_MainFunction(CDF_FilePath; Ybus_Taps_Indicator=false, NR_Type=1, Tolerance=0.001, Tol_Num=0, SortValue=true, BusSwitching=false)

    # Reading IEEE CDF File
    CDF_DF_List = CDF_Parser(CDF_FilePath, SortValue)

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
    if (Ybus_Taps_Indicator == false) # Without Taps

            Ybus = Create_Ybus_WithoutTaps(CDF_DF_List_pu)

    elseif (Ybus_Taps_Indicator == true) # With Taps

            Ybus_WithoutTaps = Create_Ybus_WithoutTaps(CDF_DF_List_pu)

            Ybus = Create_Ybus_WithTaps(Ybus_WithoutTaps,CDF_DF_List_pu)

    end

    # Creating Initial Solution Vector for Power Flow
    Initial_SolutionVector_NR = Create_Initial_SolutionVector_NR(CDF_DF_List_pu)

    # Initializing SolutionVector_NR - For getting reference outside While Loop
    SolutionVector_NR = 0

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
            
            @show WhileLoop_Counter

            # Starting Timer
            TickTock.tick()

            # If Else Loop: For checking First Iteration
            if (WhileLoop_Counter == 1)

                    # Compute Mismatch
                    CDF_DF_List_pu, Ybus, PQ_BusArray, PQ_MismatchVector, SolutionVector_V, SolutionVector_Delta = PQ_PV_Bus_Check_Modify(CDF_DF_List_pu, Ybus, Ybus_Taps_Indicator, NR_Type, Initial_SolutionVector_NR, BusSwitching)

            else

                    # Compute Mismatch
                    CDF_DF_List_pu, Ybus, PQ_BusArray, PQ_MismatchVector, SolutionVector_V, SolutionVector_Delta = PQ_PV_Bus_Check_Modify(CDF_DF_List_pu, Ybus, Ybus_Taps_Indicator, NR_Type, SolutionVector_NR, BusSwitching)

                    # Increasing  PowerFlow_IterationTimeInfo_Array Size
                    PowerFlow_IterationTimeInfo_Array = vcat(PowerFlow_IterationTimeInfo_Array, zeros(1,2))

            end

            # Create Solution Vector NR
            SolutionVector_NR = Create_SolutionVector_NR(CDF_DF_List_pu, SolutionVector_V, SolutionVector_Delta)

            # Compute Power Flow Jacobian
            Jacobian_NR = Create_Jacobian_NR(CDF_DF_List_pu, Ybus, SolutionVector_V, SolutionVector_Delta, PQ_BusArray, NR_Type, 0)

            # Compute Correction_Vector_NR using LU Factorization
            if (NR_Type == 1) # Full NR

                    Jacobian_NR_Corrected = Compute_Corrected_Matrix(Jacobian_NR)

                    Correction_Vector_NR = PLU_Solve(Jacobian_NR_Corrected, PQ_MismatchVector, Tol_Num)

                    Correction_Vector_NR_1 = copy(Correction_Vector_NR)

                    # Correcting Correction_Vector_NR
                    Correction_Vector_NR = Compute_Corrected_CorrectionVector(CDF_DF_List_pu, Correction_Vector_NR, SolutionVector_NR, NR_Type)

            elseif ((NR_Type == 2) || (NR_Type == 3)) # Decoupled/Fast Decoupled NR

                    Jacobian_NR1_Corrected = Compute_Corrected_Matrix(Jacobian_NR[1])

                    Correction_Vector_NR_Delta = PLU_Solve(Jacobian_NR1_Corrected, PQ_MismatchVector[1:(N_PQ_Bus+N_PV_Bus),1], Tol_Num)

                    Jacobian_NR2_Corrected = Compute_Corrected_Matrix(Jacobian_NR[2])

                    Correction_Vector_NR_V = PLU_Solve(Jacobian_NR2_Corrected, PQ_MismatchVector[(N_PQ_Bus+N_PV_Bus+1):end,1], Tol_Num)

                    # Creating Correction_Vector_NR
                    Correction_Vector_NR = vcat(Correction_Vector_NR_Delta, Correction_Vector_NR_V)

                    #
                    Correction_Vector_NR_1 = copy(Correction_Vector_NR)

                    # Correcting Corrections Vector
                    if (NR_Type == 2) # Decoupled NR

                        Correction_Vector_NR = Compute_Corrected_CorrectionVector(CDF_DF_List_pu, Correction_Vector_NR, SolutionVector_NR, NR_Type)

                    elseif (NR_Type == 3) # Fast Decoupled NR
                        
                        Correction_Vector_NR = Compute_Corrected_CorrectionVector(CDF_DF_List_pu, Correction_Vector_NR, SolutionVector_NR, NR_Type)

                    end

            end

            @show Correction_Vector_NR

            # Computing New SolutionVector_NR
            SolutionVector_NR = SolutionVector_NR + Correction_Vector_NR

            # Compute Tolerance Satisfaction
            Tolerance_Satisfaction = Compute_ToleranceSatisfaction(CDF_DF_List_pu, Ybus, NR_Type, Tolerance, SolutionVector_NR)

            # Stopping Timer
            IterationTime = TickTock.tok()

            # Filling-up PowerFlow_IterationTimeInfo_Array
            PowerFlow_IterationTimeInfo_Array[WhileLoop_Counter,1:2] = [WhileLoop_Counter , IterationTime]

    end

    # Compute Line Flows
    CDF_DF_List_pu, LineFlow_Array = Compute_LineFlows(CDF_DF_List_pu, Ybus, SolutionVector_NR)

    # ReOrdering BusDataCard_DF: PQ->PV->Slack
    BusDataCard_DF = CDF_DF_List_pu[2]

    sort!(BusDataCard_DF, [order(:Type_Original)])

    CDF_DF_List_pu[2] = BusDataCard_DF

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
- 'SortValue': 1 -> Sort CDF File according to Bus Type PQ->PV->Slack, any other value -> do not sort
- 'BusSwitching': 1 -> Bus Switching Employed, any other -> Bus Switching not Employed
'''
'''
# Output
- 'Predictor_Vector_History': An array consiisting of each Predictor Vector computed per iteration of Continuation Power Flow.
- 'Corrector_Vector_History': An array consiisting of each Corrector Vector computed per iteration of Continuation Power Flow.
- 'Tangent_Vector_History': An array consiisting of each Tangent Vector computed per iteration of Continuation Power Flow.
- 'PowerFlow_IterationTimeInfo_Array': An array containing Iteration Number and
Time in seconds for each iteration.
'''
"""
function ContinuationPowerFlow_MainFunction(CDF_FilePath, PQ_V_Curve_Tuple, Ybus_Taps_Indicator, StepSize_Vector_CPF, Tolerance_NR, Tol_Num, PostCriticalPoint_Counter_Input, SortValue, BusSwitching)
    
    # Reading IEEE CDF File
    CDF_DF_List = CDF_Parser(CDF_FilePath, SortValue)

    # Converting CDF DataFrame to PU
    CDF_DF_List_pu = CDF_pu_Converter(CDF_DF_List)

   # Create Ybus
   if (Ybus_Taps_Indicator == false) # Without Taps

           Ybus = Create_Ybus_WithoutTaps(CDF_DF_List_pu)

   elseif (Ybus_Taps_Indicator == true) # With Taps

           Ybus_WithoutTaps = Create_Ybus_WithoutTaps(CDF_DF_List_pu)

           Ybus = Create_Ybus_WithTaps(Ybus_WithoutTaps,CDF_DF_List_pu)

   end

   # Creating K_Vector
   K_Vector = Create_KVector_CPF(CDF_DF_List_pu, PQ_V_Curve_Tuple)

   # Creating Initial Lambda
   Lambda_Ini = 0

   # Solving Initial Power Flow Problem
   Initial_SolutionVector, PowerFlow_IterationTimeInfo_Array =  PowerFlow_MainFunction_Ini_CPF(CDF_DF_List_pu, Ybus, 1, Tolerance_NR, Tol_Num, K_Vector, Lambda_Ini)

   # Creating Initial Solution Vector
   #Initial_SolutionVector_CPF = Create_Initial_SolutionVector_CPF(CDF_DF_List_pu)
   Initial_SolutionVector_CPF = vcat(Initial_SolutionVector, Lambda_Ini)   

   # Initializing Predictor/Corrector/Tangent Vector History
   Predictor_Vector_History = zeros(size(K_Vector)[1]+1,1)
   Corrector_Vector_History = zeros(size(K_Vector)[1]+1,1)
   Tangent_Vector_History = zeros(size(K_Vector)[1]+1,1)

   # Updating Corrector_Vector_History
   Corrector_Vector_History[1:size(Corrector_Vector_History)[1], 1] = Initial_SolutionVector_CPF

   # Initializing ContinuationPowerFlow_IterationTimeInfo_Array
   ContinuationPowerFlow_IterationTimeInfo_Array = zeros(1,2)

   # Initializing PostCriticalPoint_Counter
   PostCriticalPoint_Counter = 0

   # Intializing WhileLoop_Counter
   WhileLoop_Counter = 0

   # Initializing
   Index_CPF = 0
   Index_CPF_New = 0
   Predict_b_Value_CPF = 0
   Predict_b_Value_CPF_New = 0

   # While Loop: For computing Continuation Power Flow
   while (PostCriticalPoint_Counter <= PostCriticalPoint_Counter_Input)

           # Incrementing WhileLoop_Counter
           WhileLoop_Counter = WhileLoop_Counter + 1

           # Starting Timer
           TickTock.tick()

           # If Else Loop: For checking first iteration
           if (WhileLoop_Counter == 1)

                   # Getting Solution Vector comes from Load Flow Solution
                   SolutionVector_CPF = Corrector_Vector_History[:,WhileLoop_Counter]

                   # Getting Index k (Choosing Lambda)
                   Index_CPF = size(K_Vector)[1]+1

                   # Initializing Predict_b_Value_CPF
                   Predict_b_Value_CPF = 1

           else
                   # Getting Solution Vector From Corrector
                   SolutionVector_CPF = Corrector_Vector_History[:,WhileLoop_Counter]

                   # Getting Index k from previous iteration
                   Index_CPF = Index_CPF_New

                   # Initializing Predict_b_Value_CPF with Predict_b_Value_CPF_New
                   Predict_b_Value_CPF = Predict_b_Value_CPF_New

                   # Increasing  ContinuationPowerFlow_IterationTimeInfo_Array Size
                   ContinuationPowerFlow_IterationTimeInfo_Array = vcat(ContinuationPowerFlow_IterationTimeInfo_Array, zeros(1,2))

           end

           # Choosing correct step size based on current continuation parameterization
           if (Index_CPF == size(SolutionVector_CPF)[1])  # Lambda is the Continuation Parameter

                StepSize_CPF = StepSize_Vector_CPF[1]

           else # Lambda is not the Continuation Parameter

                StepSize_CPF = StepSize_Vector_CPF[2]

           end

           # Compute Continuation Power Flow Jacobian Predictor Step (Debug remove multiplication by voltage)
           Jacobian_CPF_Predict = Create_Jacobian_CPF_Predict(CDF_DF_List_pu, Ybus, SolutionVector_CPF, 1, K_Vector, Index_CPF)

           # Compute Tangent Vector
           Tangent_Vector = Compute_Tangent_Vector_CPF(Jacobian_CPF_Predict, Tol_Num, Predict_b_Value_CPF, Index_CPF)

           Lambda_T = Tangent_Vector[end,1]
           @show Lambda_T

           # Check if Critical Point is passed
           PostCriticalPoint_Counter = Check_CriticalPoint_CPF(Tangent_Vector, PostCriticalPoint_Counter)

           # Choose Continuation Parameter
           Index_CPF_New, Predict_b_Value_CPF_New = Choose_ContinuationParameter_CPF(Tangent_Vector)

           # Correcting Tangent Vector for Delts - Rad to Degree
           Tangent_Vector_Corrected = Compute_Corrected_TangentVector_CPF(CDF_DF_List_pu, Tangent_Vector)

           # Predict Solution
           CPF_Predictor_Vector = Compute_PredictVector_CPF(SolutionVector_CPF, Tangent_Vector_Corrected, StepSize_CPF)

           # Compute Continuation Power Flow Jacobian Corrector Step
           # Jacobian_CPF_Correct = Create_Jacobian_CPF_Correct(CDF_DF_List_pu, Ybus, CPF_Predictor_Vector, 1, K_Vector, Index_CPF)

           # Update/Correct Solution
           CPF_Corrector_Vector, PowerFlow_IterationTimeInfo_Array_Corrector = PowerFlow_MainFunction_CPF(CDF_DF_List_pu, Ybus, 1, CPF_Predictor_Vector, Tolerance_NR, Tol_Num, K_Vector, Index_CPF)

           Lambda_C = CPF_Corrector_Vector[end,1]
           @show Lambda_C

           # Update Predictor/Corrector Vector History
           Predictor_Vector_History = hcat(Predictor_Vector_History, CPF_Predictor_Vector)
           Corrector_Vector_History = hcat(Corrector_Vector_History, CPF_Corrector_Vector)
           Tangent_Vector_History = hcat(Tangent_Vector_History, Tangent_Vector)

           # Stopping Timer
           IterationTime = TickTock.tok()

           # Filling-up ContinuationPowerFlow_IterationTimeInfo_Array
           ContinuationPowerFlow_IterationTimeInfo_Array[WhileLoop_Counter,1:2] = [WhileLoop_Counter , IterationTime]

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
- 'Measurement_Error_Variance' - Vector of measurement error variance in order V -> P -> Q.
- 'Bad_Bus_Measurement_Input': Array of format -  [[Bus_Num, V_scale, P_scale, Q_scale], []] when bad data present, 
otherwise 'nothing'; the measurements are scaled through scale variables to make bad data.
- 'Bad_Branch_Measurement_Input': Array of format - [[Bus_Num_i, Bus_Num_j, P_+_scale, P_-_scale, Q_+_scale, Q_-_scale], []] when bad data present, 
otherwise 'nothing'; the measurements are scaled through scale variables to make bad data.
- 'alpha': Alpha value for Chi-Square Distribution computation for Bad Data detection
- 'Ybus_Taps_Indicator': 1 - Ybus with no taps, 2 - Ybus with Taps
- 'Tolerance': Tolerance level for stopping criterion of Newton-Raphson Method.
- 'Tol_Num': Tolerance for being near zero
- 'SortValue': 1 -> Sort CDF File according to Bus Type PQ->PV->Slack, any other value -> do not sort
- 'BusSwitching': 1 -> Bus Switching Employed, any other -> Bus Switching not Employed
'''
'''
# Output
- '':
'''
"""
function PowerSystem_StateEstimation_MainFunction(CDF_FilePath, Measurement_Error_Variance, Bad_Bus_Measurement_Input, Bad_Branch_Measurement_Input, Tolerance_SE, alpha, Ybus_Taps_Indicator, Tolerance_NR, Tol_Num, SortValue, BusSwitching)

        # Solving Initial Power Flow Problem
        CDF_DF_List_pu, LineFlow_Array, PowerFlow_IterationTimeInfo_Array =  PowerFlow_MainFunction(CDF_FilePath; Ybus_Taps_Indicator=Ybus_Taps_Indicator, NR_Type=1, Tolerance=Tolerance_NR, Tol_Num=Tol_Num, SortValue=SortValue, BusSwitching=BusSwitching)

        # Create Ybus
        if (Ybus_Taps_Indicator == false) # Without Taps

                Ybus = Create_Ybus_WithoutTaps(CDF_DF_List_pu)

        elseif (Ybus_Taps_Indicator == true) # With Taps

                Ybus_WithoutTaps = Create_Ybus_WithoutTaps(CDF_DF_List_pu)

                Ybus = Create_Ybus_WithTaps(Ybus_WithoutTaps,CDF_DF_List_pu)

        end

        # Creating Initial Solution Vector
        Initial_SolutionVector_SE = Create_Initial_SolutionVector1_SE(CDF_DF_List_pu)

        # Getting Solution Vectors divided in Voltages and Deltas
        SolutionVector_V_Ini, SolutionVector_Delta_Ini = Create_SolutionVector_VDelta_SE(CDF_DF_List_pu, Initial_SolutionVector_SE)

        # Compute Bus Injections
        Bus_PQ_Array = Compute_PQ_BusArray(Ybus, SolutionVector_V_Ini, SolutionVector_Delta_Ini)        
        
        # Create Bus Branch Measurement Arrays
        Bus_Measurement_Array, Branch_Measurement_Array = Create_BusBranch_Measurement_Arrays_SE(CDF_DF_List_pu, LineFlow_Array, Bus_PQ_Array, Measurement_Error_Variance)

        # Create Bad Bus/Branch Measurement Arrays
        Bad_Bus_Measurement_Array, Bad_Branch_Measurement_Array = Create_Bad_BusBranch_Measurement_Arrays_SE(Bus_Measurement_Array, Branch_Measurement_Array, Bad_Bus_Measurement_Input, Bad_Branch_Measurement_Input)

        # Add Bus_Measurement_Array, Branch_Measurement_Array to CDF_DF_List
        CDF_DF_List_pu = CDF_AddMeasurements_SE(CDF_DF_List_pu, Bad_Bus_Measurement_Array, Bad_Branch_Measurement_Array)    
        
        # Creating Incidence Matrix (A) and Initial Detected Bad Data Vector
        IncidenceMatrix_A, Detected_BadData_Vector_Ini = Compute_AMat_BadDataVec_SE(CDF_DF_List_pu)

        # Initializing Bad Data Debug_Indicator
        Bad_Data_Indicator = true

        Detected_BadData_Vector = Detected_BadData_Vector_Ini 

        # Initializing StateEstimation_IterationTimeInfo_Array
        StateEstimation_IterationTimeInfo_Array = zeros(1,2)

        # Initializing While Loop Counter
        WhileLoop_Counter = 0

        # While Loop: For State Estimation and Bad Data Detection
        while (Bad_Data_Indicator == true)

                # Incrementing While Loop Counter
                WhileLoop_Counter = WhileLoop_Counter + 1

                # Starting Timer
                TickTock.tick()

                # If Else Loop: For checking first iteration
                if (WhileLoop_Counter == 1)

                        # Getting Solution Vector comes from Load Flow Solution
                        SolutionVector_SE = Initial_SolutionVector_SE

                else
                        # Getting Solution Vector from Load Flow Solution
                        SolutionVector_SE = Initial_SolutionVector_SE

                        # Increasing  ContinuationPowerFlow_IterationTimeInfo_Array Size
                        StateEstimation_IterationTimeInfo_Array = vcat(StateEstimation_IterationTimeInfo_Array, zeros(1,2))

                end

                # Perform State Estimation
                State_Estimate, CDF_DF_List_pu, R_Inv_Matrix, Z_Measured_Vector  = Compute_StateEstimation_SE(CDF_DF_List_pu, SolutionVector_SE, Ybus, IncidenceMatrix_A, Detected_BadData_Vector, Tolerance_SE)

                # Perform Bad Data Detection
                Bad_Data_Indicator, Detected_BadData_Vector = Compute_Bad_Data_Detection_SE(CDF_DF_List_pu, State_Estimate, R_Inv_Matrix, IncidenceMatrix_A, Z_Measured_Vector, Ybus, Detected_BadData_Vector, alpha)

                # Stopping Timer
                IterationTime = TickTock.tok()

                # Filling-up ContinuationPowerFlow_IterationTimeInfo_Array
                StateEstimation_IterationTimeInfo_Array[WhileLoop_Counter,1:2] = [WhileLoop_Counter , IterationTime]


        end

        # return State_Estimate, CDF_DF_List_pu, StateEstimation_IterationTimeInfo_Array  
 
        return nothing

end

# Including component files
include("IEEE_CDF_Parser.jl")

include("Ybus_Builder.jl")

include("BasicPowerFlow_Functions.jl")

include("Jacobian_Builder.jl")

include("LU_Factorization.jl")

include("BasicContinuationPowerFlow_Functions.jl")

include("BasicStateEstimation_Functions.jl")

include("Helper_Functions.jl")


end

