# BasicOPF_Functions.jl

""" 
Compute_EconomicDispatch_A_Matrix(GeneratorCostCurve_Array)

Computes A matrix for economic dispatch problem.

'''     
# Arguments 

'''
'''
# Output

'''
"""
function Compute_EconomicDispatch_A_Matrix(GeneratorCostCurve_Array)

    # Getting the length of GeneratorCostCurve_Array
    Len_GeneratorCostCurve_Array = length(GeneratorCostCurve_Array)

    # Initializing the EconomicDispatch_A_Matrix
    EconomicDispatch_A_Matrix = zeros(Len_GeneratorCostCurve_Array+1, Len_GeneratorCostCurve_Array+1)

    # Adding the -1 end column 
    EconomicDispatch_A_Matrix[1:end-1,end] = -1 * ones(Len_GeneratorCostCurve_Array, 1)

    # Adding the 1 end row
    EconomicDispatch_A_Matrix[end,1:end-1] = ones(1, Len_GeneratorCostCurve_Array)

    # Adding the diagonal elements
    for ii in 1:Len_GeneratorCostCurve_Array  # Through the rows of EconomicDispatch_A_Matrix

        for jj in 1:Len_GeneratorCostCurve_Array  # Through the columns of EconomicDispatch_A_Matrix

            if (ii == jj) # Diagonal Element

                EconomicDispatch_A_Matrix[ii,jj] =  2 * GeneratorCostCurve_Array[ii][1]

            end

        end

    end

    return EconomicDispatch_A_Matrix 

end

""" 
Compute_EconomicDispatch_b_Vector(GeneratorCostCurve_Array, Load_Demand)

Computes b vector for economic dispatch problem.

'''     
# Arguments 

'''
'''
# Output

'''
"""
function Compute_EconomicDispatch_b_Vector(GeneratorCostCurve_Array, Load_Demand)

    # Getting the length of GeneratorCostCurve_Array
    Len_GeneratorCostCurve_Array = length(GeneratorCostCurve_Array)

    # Initializing the EconomicDispatch_b_Vector
    EconomicDispatch_b_Vector = zeros(Len_GeneratorCostCurve_Array+1, 1)

    # Adding the load demand at the end row
    EconomicDispatch_b_Vector[end,1] = Load_Demand

    # Adding the rest of the elements
    for ii in 1:Len_GeneratorCostCurve_Array

        EconomicDispatch_b_Vector[ii,1] = -GeneratorCostCurve_Array[ii][2]

    end

    return EconomicDispatch_b_Vector

end

"""
    Create_SolutionVector_OPF(CDF_DF_List_pu)

Creates Solution Vector for the Power System Optimal Power Flow.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
'''
'''
# Output
- 'Initial_SolutionVector_x': Power Flow Solution Voltage and Angle at each
bus ordered according to bus type: PQ->PV.
'''
"""
function Create_SolutionVector_x_OPF(CDF_DF_List_pu)

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]

        # Number of Buses
        N_Bus = nrow(BusDataCard_DF)
        N_PQ_Bus = nrow(filter(row -> ((row.Type_Original == 0) || (row.Type_Original == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type_Original == 2), BusDataCard_DF))
        N_Slack_Bus = nrow(filter(row -> (row.Type_Original == 3), BusDataCard_DF))

        # Creating SolutionVector_OPF Length
        SolutionVector_OPF_Len = (2*N_PQ_Bus)+(N_PV_Bus)

        # Initializing SolutionVector_OPF
        SolutionVector_OPF = zeros(SolutionVector_OPF_Len, 1)

        # Creating SolutionVector_OPF Delta part for both PQ and PV Buses
        SolutionVector_OPF[1:(N_PQ_Bus+N_PV_Bus),1] =  BusDataCard_DF.Final_A_deg[1:(N_PQ_Bus+N_PV_Bus),1]

        # Creating SolutionVector_OPF V part for PQ Buses
        SolutionVector_OPF[(N_PQ_Bus+N_PV_Bus)+1:end,1] =  BusDataCard_DF.Final_V_pu[1:N_PQ_Bus,1]

        return SolutionVector_OPF

end

""" 
Get_BaseMVA_OPF(CDF_DF_List_pu)

Gets system base MVA from CDF file.

'''     
# Arguments 

'''
'''
# Output

'''
"""
function Get_BaseMVA_OPF(CDF_DF_List_pu)

    # Get Title Card from CDF_DF_List_pu
    TitleCard_DF = CDF_DF_List_pu[1]

    # Get Base MVA from Title Card
    Base_MVA = TitleCard_DF.MVA_Base[1]

    return Base_MVA 

end

""" 
Get_Generator_BusNum_CostCurve_Arrays_OPF(CDF_DF_List_pu, Generator_BusNum_CostCurve_Array)

Gets system base MVA from CDF file.

'''     
# Arguments 

'''
'''
# Output

'''
"""
function Get_Generator_BusNum_CostCurve_Arrays_OPF(CDF_DF_List_pu, Generator_BusNum_CostCurve_Array)

    # Getting required data from CDF_DF_List
    BusDataCard_DF = CDF_DF_List_pu[2]

    # Number of Buses
    N_Bus = nrow(BusDataCard_DF)

    # Getting length of Generator_BusNum_CostCurve_Array
    Len_Generator_BusNum_CostCurve_Array = length(Generator_BusNum_CostCurve_Array)

    # Initializing Generator_CostCurve_Matrix
    Generator_CostCurve_Matrix = zeros(Len_Generator_BusNum_CostCurve_Array, 4)

    # Initializing u_Index_Vector
    u_Index_Vector = zeros(Len_Generator_BusNum_CostCurve_Array,1)

    # Filling up Generator_CostCurve_Array and u_Index_Vector
    for ii in 1:Len_Generator_BusNum_CostCurve_Array  # For each element in Generator_BusNum_CostCurve_Array

        for jj in 1:N_Bus  # For each bus in BusDataCard_DF

            if (Int64(Generator_BusNum_CostCurve_Array[ii][1]) == BusDataCard_DF.Bus_Num[jj])

                if (jj == N_Bus)  # Slack Bus

                    # Updating Generator_CostCurve_Matrix
                    Generator_CostCurve_Matrix[ii,1] = 1
                    Generator_CostCurve_Matrix[ii,2] = Generator_BusNum_CostCurve_Array[ii][2]
                    Generator_CostCurve_Matrix[ii,3] = Generator_BusNum_CostCurve_Array[ii][3]
                    Generator_CostCurve_Matrix[ii,4] = Generator_BusNum_CostCurve_Array[ii][4]

                else  # PV Bus

                    # Updating Generator_CostCurve_Matrix
                    Generator_CostCurve_Matrix[ii,1] = jj+1
                    Generator_CostCurve_Matrix[ii,2] = Generator_BusNum_CostCurve_Array[ii][2]
                    Generator_CostCurve_Matrix[ii,3] = Generator_BusNum_CostCurve_Array[ii][3]
                    Generator_CostCurve_Matrix[ii,4] = Generator_BusNum_CostCurve_Array[ii][4]

                end

            end

        end

    end

    # Sort Generator_CostCurve_Matrix w.r.t Column 1 - Bus Order 
    Generator_CostCurve_Matrix[sortperm(Generator_CostCurve_Matrix[:, 1]), :]

    # Get u_Index_Vector
    u_Index_Vector[:,1] = Generator_CostCurve_Matrix[:,1]

    # Getting u_P_Index_Vector
    if (u_Index_Vector[1,1] == 1.0)  # Slack Bus Present

        u_P_Index_Vector = u_Index_Vector[2:end, 1]

    else # Slack Bus absent

        u_P_Index_Vector = u_Index_Vector

    end
    
    # Get Just Generator Cost parameters
    Generator_CostCurve_Matrix_New = Generator_CostCurve_Matrix[:,2:end]

    return Generator_CostCurve_Matrix_New, u_Index_Vector, u_P_Index_Vector

end

""" 
Get_Line_Bus_Index_PLimit_Arrays_OPF(CDF_DF_List_pu, Line_PowerFlow_Limit_Array, Base_MVA)

Gets system base MVA from CDF file.

'''     
# Arguments 

'''
'''
# Output

'''
"""
function Get_Line_Bus_Index_PLimit_Arrays_OPF(CDF_DF_List_pu, Line_PowerFlow_Limit_Array, Base_MVA)

    if (Line_PowerFlow_Limit_Array == nothing)  # No Line Limits present

        # Initializing Line_Index_Vector, Line_Bus_Index_Matrix, and Line_P_Limit_Vector
        Line_Index_Vector = nothing
        Line_Bus_Index_Matrix = nothing   
        Line_P_Limit_Vector = nothing

    else  # Line Limits present

        # Getting required data from CDF_DF_List
        BusDataCard_DF = CDF_DF_List_pu[2]
        BranchDataCard_DF = CDF_DF_List_pu[3]

        # Number of Buses and Branches
        N_Bus = nrow(BusDataCard_DF)
        N_Branch = nrow(BranchDataCard_DF)

        # Getting size of Line_PowerFlow_Limit_Array
        Len_Line_PowerFlow_Limit_Array = size(Line_PowerFlow_Limit_Array)[1]

        # Initializing Line_Index_Vector, Line_Bus_Index_Matrix, and Line_P_Limit_Vector
        Line_Index_Vector = zeros(Len_Line_PowerFlow_Limit_Array,1)
        Line_Bus_Index_Matrix = zeros(Len_Line_PowerFlow_Limit_Array,2)    
        Line_P_Limit_Vector = zeros(Len_Line_PowerFlow_Limit_Array,1)

        # Updating Line_Index_Vector
        for ii in 1:Len_Line_PowerFlow_Limit_Array  # For each element in Line_PowerFlow_Limit_Array

            # Getting Bus i and Bus j numbers
            Bus_i = Line_PowerFlow_Limit_Array[ii][1]
            Bus_j = Line_PowerFlow_Limit_Array[ii][2]

            for jj in 1:N_Branch  # For each element in BranchDataCard_DF

                # Checking if Current Bus numbers match
                if (((Bus_i == BranchDataCard_DF.Tap_Bus_Num[jj]) && (Bus_j == BranchDataCard_DF.Z_Bus_Num[jj])) || ((Bus_j == BranchDataCard_DF.Tap_Bus_Num[jj]) && (Bus_i == BranchDataCard_DF.Z_Bus_Num[jj])))
                    
                    Line_Index_Vector[ii,1] = jj

                end

            end

        end

        # Updating Line_Bus_Index_Matrix
        for ii in 1:Len_Line_PowerFlow_Limit_Array  # # For each element in Line_PowerFlow_Limit_Array

            # Getting Bus i and Bus j numbers
            Bus_i = Line_PowerFlow_Limit_Array[ii][1]
            Bus_j = Line_PowerFlow_Limit_Array[ii][2]

            # Getting New Bus i and Bus j numbers (if the branch is numbered in reverse)
            Bus_i_New  = 0  # Initialization
            Bus_j_New  = 0  # Initialization
            for jj in 1:N_Branch  # For each element in BranchDataCard_DF

                # Checking if Current Bus numbers match
                if ((Bus_i == BranchDataCard_DF.Tap_Bus_Num[jj]) && (Bus_j == BranchDataCard_DF.Z_Bus_Num[jj])) 
                    
                    Bus_i_New = Bus_i
                    Bus_j_New = Bus_j

                elseif ((Bus_j == BranchDataCard_DF.Tap_Bus_Num[jj]) && (Bus_i == BranchDataCard_DF.Z_Bus_Num[jj]))

                    Bus_i_New = Bus_j
                    Bus_j_New = Bus_i

                end

            end

            # Getting Correct Bus i Index for updating Line_Bus_Index_Matrix based on New Bus i 
            for jj in 1:N_Bus  # For each element in BusDataCard_DF

                # Checking if Current Bus numbers match
                if (Bus_i_New == BusDataCard_DF.Bus_Num[jj])
                    
                    if (jj == N_Bus)  # Slack Bus

                        Line_Bus_Index_Matrix[ii,1] = 1

                    else # Not Slack Bus

                        Line_Bus_Index_Matrix[ii,1] = jj + 1

                    end

                end

            end

            # Getting Correct Bus i Index for updating Line_Bus_Index_Matrix
            for jj in 1:N_Bus  # For each element in BusDataCard_DF

                # Checking if Current Bus numbers match
                if (Bus_j_New == BusDataCard_DF.Bus_Num[jj])
                    
                    if (jj == N_Bus)  # Slack Bus

                        Line_Bus_Index_Matrix[ii,2] = 1

                    else # Not Slack Bus

                        Line_Bus_Index_Matrix[ii,2] = jj + 1

                    end

                end

            end

        end

        # Updating Line_P_Limit_Vector
        for ii in 1:Len_Line_PowerFlow_Limit_Array  # For each element in Line_PowerFlow_Limit_Array

            Line_P_Limit_Vector[ii,1] = Line_PowerFlow_Limit_Array[ii][3] / Base_MVA

        end

    end

    return Line_Index_Vector, Line_Bus_Index_Matrix, Line_P_Limit_Vector

end

""" 
Create_SolutionVector_p_OPF(CDF_DF_List_pu,, u_P_Index_Vector)

Gets system base MVA from CDF file.

'''     
# Arguments 

'''
'''
# Output

'''
"""
function Create_SolutionVector_p_OPF(CDF_DF_List_pu, u_P_Index_Vector)

    # Getting required data from CDF_DF_List
    BusDataCard_DF = CDF_DF_List_pu[2]

    # Number of Buses
    N_Bus = nrow(BusDataCard_DF)

    # Initialize SolutionVector_p and SolutionVector_p_Full
    SolutionVector_p = zeros(size(u_P_Index_Vector)[1], 1)
    
    SolutionVector_p_Full = zeros(N_Bus, 1)

    # Getting SolutionVector_p
    for ii in 1:size(u_P_Index_Vector)[1]  # For each element in u_P_Index_Vector

        for jj in 1:N_Bus  # For each Bus

            if (Int64(u_P_Index_Vector[ii, 1]) == jj+1)

                SolutionVector_p[ii,1] = BusDataCard_DF.Gen_MW[jj]
            
            end

        end

    end

    # Getting SolutionVector_p_Full
    for ii in 1:N_Bus  # For each Bus

        if (ii == N_Bus)  # Slack-Bus

            SolutionVector_p_Full[1,1] = BusDataCard_DF.Gen_MW[ii]

        else  # PQ-PV Buses

            SolutionVector_p_Full[ii+1,1] = BusDataCard_DF.Gen_MW[ii]

        end

    end
    
    return SolutionVector_p, SolutionVector_p_Full

end

""" 
Compute_P_MW(P_pu_Vector, Base_MVA)

Gets system base MVA from CDF file.

'''     
# Arguments 

'''
'''
# Output

'''
"""
function Compute_P_MW(P_pu_Vector, Base_MVA)

    P_MW_Vector = Base_MVA * P_pu_Vector

    return P_MW_Vector

end

"""
    PowerFlow_MainFunction_OPF(CDF_FilePath, Ybus_Taps_Indicator, NR_Type, Tolerance, Tol_Num)

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
function PowerFlow_MainFunction_OPF(CDF_DF_List_pu, Ybus, Ybus_Taps_Indicator, NR_Type, Tolerance, Tol_Num, SortValue, BusSwitching)    

    # Getting required data from CDF_DF_List
    BusDataCard_DF = CDF_DF_List_pu[2]

    # Number of Buses
    N_Bus = nrow(BusDataCard_DF)
    N_PQ_Bus = nrow(filter(row -> ((row.Type == 0) || (row.Type == 1)), BusDataCard_DF))
    N_PV_Bus = nrow(filter(row -> (row.Type == 2), BusDataCard_DF))
    N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

    N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

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

    # Compute Bus Injections

    # Creating SolutionVector_V and SolutionVector_Delta
    SolutionVector_V, SolutionVector_Delta = Create_SolutionVector_VDelta_NR(CDF_DF_List_pu, SolutionVector_NR)

    #  Compute P-Q at Buses
    PQ_BusArray = Compute_PQ_BusArray(Ybus, SolutionVector_V, SolutionVector_Delta)

    # Compute Line Flows
    CDF_DF_List_pu, LineFlow_Array = Compute_LineFlows(CDF_DF_List_pu, Ybus, SolutionVector_NR)

    # ReOrdering BusDataCard_DF: PQ->PV->Slack
    BusDataCard_DF = CDF_DF_List_pu[2]

    sort!(BusDataCard_DF, [order(:Type_Original)])

    CDF_DF_List_pu[2] = BusDataCard_DF

    return CDF_DF_List_pu, LineFlow_Array, PQ_BusArray, PowerFlow_IterationTimeInfo_Array

end

"""
    Compute_ToleranceSatisfaction_OPF(Tolerance_OPF, Corrections_C_Vector)

Computes if the correction vector satifies the tolerance.

'''
# Arguments
- '': 
'''
# Output
- '': 
'''
"""
function Compute_ToleranceSatisfaction_OPF(Tolerance_OPF, Corrections_C_Vector)
        
       
    # Initializing Tolerance Counter
    ToleranceCounter = 0

    # Computing Tolerance Counter
    for ii in 1: length(Corrections_C_Vector)

            if (abs(Corrections_C_Vector[ii]) > Tolerance_OPF)

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
Update_CDFDFListpu_uP_OPF(CDF_DF_List_pu, SolutionVector_p, u_P_Index_Vector)

Sets the P Generation in CDF DF Bus Data Card to values computed in optimal power flow iteration.

'''
# Arguments
- '': 
'''
# Output
- '': 
'''
"""
function Update_CDFDFListpu_uP_OPF(CDF_DF_List_pu, SolutionVector_p, u_P_Index_Vector)

    # Getting required data from CDF_DF_List
    BusDataCard_DF = CDF_DF_List_pu[2]

    # Updating BusDataCard_DF
    for ii in 1:size(u_P_Index_Vector)[1]  # For each element in u_P_Index_Vector

        BusDataCard_DF.Gen_MW[Int64(u_P_Index_Vector[ii,1])-1] = SolutionVector_p[ii,1]

    end

    return CDF_DF_List_pu

end

"""
Update_CDFDFListpu_SlackP_OPF(CDF_DF_List_pu, PQ_BusArray)

Sets the P Generation for Slack Bus in CDF DF Bus Data Card to values computed in converged Power Flow iteration.

'''
# Arguments
- '': 
'''
# Output
- '': 
'''
"""
function Update_CDFDFListpu_SlackP_OPF(CDF_DF_List_pu, PQ_BusArray)

    # Getting required data from CDF_DF_List
    BusDataCard_DF = CDF_DF_List_pu[2]
    
    # Number of Buses
    N_Bus = nrow(BusDataCard_DF)

    # Getting Slack Bus P injected
    Slack_P_Inj = PQ_BusArray[1,1]

    # Getting Slack Bus P demand
    Slack_P_Demand = BusDataCard_DF.Load_MW[N_Bus]

    # Update BusDataCard_DF Slack Bus P Generated
    BusDataCard_DF.Gen_MW[N_Bus] = Slack_P_Inj + Slack_P_Demand

    return CDF_DF_List_pu 

end

"""
Compute_GenerationCost_OPF(SolutionVector_p, SolutionVector_p_Full, Generator_CostCurve_Matrix_New, Base_MVA)

Sets the P Generation for Slack Bus in CDF DF Bus Data Card to values computed in converged Power Flow iteration.

'''
# Arguments
- '': 
'''
# Output
- '': 
'''
"""
function Compute_GenerationCost_OPF(SolutionVector_p, SolutionVector_p_Full, Generator_CostCurve_Matrix_New, Base_MVA)

    # Get size of SolutionVector_p
    Len_SolutionVector_p = size(SolutionVector_p)[1]

    # Get size of Generator_CostCurve_Matrix_New
    Len_Generator_CostCurve_Matrix_New = size(Generator_CostCurve_Matrix_New)[1]

    # Convert SolutionVector_p from pu to MW
    SolutionVector_p_MW = Base_MVA * SolutionVector_p

    # Convert SolutionVector_p_Full from pu to MW
    SolutionVector_p_Full_MW = Base_MVA * SolutionVector_p_Full

    # Initializing CostGeneration
    CostGeneration = 0

    # Compute Cost of Generation
    if (Len_SolutionVector_p != Len_Generator_CostCurve_Matrix_New)  # Slack Bus is a part of Cost Function

        for ii in 1:Len_Generator_CostCurve_Matrix_New  # For each element in Len_Generator_CostCurve_Matrix_New

            if (ii == 1)  # Firstirst iteration

                # Updating CostGeneration
                CostGeneration = CostGeneration + ((Generator_CostCurve_Matrix_New[ii,1] * SolutionVector_p_Full_MW[1,1]^(2)) + (Generator_CostCurve_Matrix_New[ii,2] * SolutionVector_p_Full_MW[1,1]) + (Generator_CostCurve_Matrix_New[ii,3]))

            else  # Not the first iteration

                # Updating CostGeneration
                CostGeneration = CostGeneration + ((Generator_CostCurve_Matrix_New[ii,1] * SolutionVector_p_MW[ii-1,1]^(2)) + (Generator_CostCurve_Matrix_New[ii,2] * SolutionVector_p_MW[ii-1,1]) + (Generator_CostCurve_Matrix_New[ii,3]))

            end
            
        end

    else  # Slack Bus is not part of cost function

        for ii in 1:Len_SolutionVector_p  # For each element in SolutionVector_p_MW

            # Updating CostGeneration
            CostGeneration = CostGeneration + ((Generator_CostCurve_Matrix_New[ii,1] * SolutionVector_p_MW[ii,1]^(2)) + (Generator_CostCurve_Matrix_New[ii,2] * SolutionVector_p_MW[ii,1]) + (Generator_CostCurve_Matrix_New[ii,3]))

        end

    end

    return CostGeneration

end

"""
Compute_LineLimit_Violation_Indicator_OPF(LineFlow_Array, Line_Index_Vector, Line_P_Limit_Vector)

Sets the P Generation for Slack Bus in CDF DF Bus Data Card to values computed in converged Power Flow iteration.

'''
# Arguments
- '': 
'''
# Output
- '': 
'''
"""
function Compute_LineLimit_Violation_Indicator_OPF(LineFlow_Array, Line_Index_Vector, Line_Bus_Index_Matrix,  Line_P_Limit_Vector)

    # Getting Length of Line_Index_Vector
    Len_Line_Index_Vector  = size(Line_Index_Vector)[1]

    # Initializing Line_LimitViolated_Index_Vector
    Line_LimitViolated_Index_Vector = zeros(1,1)

    # Initializing Line_LimitViolated_Bus_Index_Matrix
    Line_LimitViolated_Bus_Index_Matrix = zeros(1,2)

    # Initializing Line_LimitViolated_P_Limit_Vector
    Line_LimitViolated_P_Limit_Vector = zeros(1,1)

    # Initializing LineLimit_Violation_Counter
    LineLimit_Violation_Counter = 0

    # Checking for Line Limit Violations
    for ii in 1:Len_Line_Index_Vector  # For each element in Line_Index_Vector

        if (LineFlow_Array[Int64(Line_Index_Vector[ii,1]),1] > Line_P_Limit_Vector[ii,1])  # Line Limit is violated 

            # Updating LineLimit_Violation_Counter
            LineLimit_Violation_Counter = LineLimit_Violation_Counter + 1

            if (LineLimit_Violation_Counter == 1)  # First Iteration

                # Updating Line_Index_Vector_New
                Line_LimitViolated_Index_Vector[1,1] = Line_Index_Vector[ii,1]

                # Updating Line_LimitViolated_Bus_Index_Matrix
                Line_LimitViolated_Bus_Index_Matrix[1,:] = Line_Bus_Index_Matrix[ii,:]

                # Updating Line_LimitViolated_P_Limit_Vector
                Line_LimitViolated_P_Limit_Vector[1,1] = Line_P_Limit_Vector[ii,1]

            else  # Greater than first iteration

                # Updating Line_Index_Vector_New
                Line_LimitViolated_Index_Vector = vcat(Line_LimitViolated_Index_Vector, Line_Index_Vector[ii,1])

                # Updating Line_LimitViolated_Bus_Index_Matrix
                Line_LimitViolated_Bus_Index_Matrix = vcat(Line_LimitViolated_Bus_Index_Matrix, Line_Bus_Index_Matrix[ii,:])

                # Updating Line_LimitViolated_P_Limit_Vector
                Line_LimitViolated_P_Limit_Vector = vcat(Line_LimitViolated_P_Limit_Vector, Line_P_Limit_Vector[ii,1])

            end
            
        end

    end

    # Computing LineLimit_Violation_Indicator
    if (LineLimit_Violation_Counter > 0)  # Line Limit Violations present

        LineLimit_Violation_Indicator = true

    else  # Line Limit Violations not present

        LineLimit_Violation_Indicator = false

    end
    
    return LineLimit_Violation_Indicator, Line_LimitViolated_Index_Vector, Line_LimitViolated_Bus_Index_Matrix, Line_LimitViolated_P_Limit_Vector

end