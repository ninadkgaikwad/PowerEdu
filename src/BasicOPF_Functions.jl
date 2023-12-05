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
        SolutionVector_OPF[(N_PQ_Bus+N_PV_Bus)+1:end,1] =  BusDataCard_DF.Final_V_dpu[(N_PQ_Bus+N_PV_Bus)+1:end-1,1]

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

    return Base_MVA 

end