include("PowerSystemsAnalysis.jl")

using .PowerSystemsAnalysis

using DataFrames

using DelimitedFiles

folderInput = "data/"
folder_processedData = "processedData/";
systemName = "IEEE_14";
# systemName = "IEEE_30";

SortValue = 1

# createFolderIfNotExisting(systemName, folder_processedData)

fileType_CDFFile = ".txt";
filename_CDFFile = folderInput*systemName*"/"*systemName*"_Data"*fileType_CDFFile

# Which Code to Debug
Debug_Indicator = 4

if (Debug_Indicator == 1) # LU Factorization

    ## Debugging LU Factorization

    #A = [2 -1 -2; -4 6 3; -4 -2 8]

    #A = [1 2. 5; 3 8 0; -1 3 4]

    #b = [1; 1. ;6]

    A = [2 -6 -1; -3. -1 7; -8 1 -2]

    b = [1.; 0 ;0]

    Tol_Num = 0.0001


    #P, L, U = PowerSystemsAnalysis.Compute_PLU(A, Tol_Num)   

    x = PowerSystemsAnalysis.PLU_Solve(A, b, Tol_Num) 

elseif (Debug_Indicator == 2) # CDF Parser

    ## Debugging CDF Parser

    CDF_FilePath = "C:/Users/ninad/Dropbox (Personal)/NinadGaikwad_PhD/Gaikwad_Research/Gaikwad_Research_Work/PowerSystemsAnalysis/data/IEEE14/IEEE_14_Data.txt"
    CDF_DF_List = PowerSystemsAnalysis.CDF_Parser(CDF_FilePath)

    CDF_DF_List_pu = PowerSystemsAnalysis.CDF_pu_Converter(CDF_DF_List)

elseif (Debug_Indicator == 3) # YBus Ybus_Builder

    # Debugging YBus Builder
    CDF_FilePath = "C:/Users/ninad/Dropbox (Personal)/NinadGaikwad_PhD/Gaikwad_Research/Gaikwad_Research_Work/PowerSystemsAnalysis/data/IEEE14/IEEE_14_Data.txt"

    Ybus_Taps_Indicator = false

    # Reading IEEE CDF File
    CDF_DF_List = PowerSystemsAnalysis.CDF_Parser(CDF_FilePath)

    # Converting CDF DataFrame to PU
    CDF_DF_List_pu = PowerSystemsAnalysis.CDF_pu_Converter(CDF_DF_List)

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

            Ybus = PowerSystemsAnalysis.Create_Ybus_WithoutTaps(CDF_DF_List_pu)

    elseif (Ybus_Taps_Indicator == 2) # With Taps

            Ybus_WithoutTaps = PowerSystemsAnalysis.Create_Ybus_WithoutTaps(CDF_DF_List_pu)

            Ybus = PowerSystemsAnalysis.Create_Ybus_WithTaps(Ybus_WithoutTaps,CDF_DF_List_pu)

    end

elseif (Debug_Indicator == 4) # Full Power Flow

     # Debugging YBus Builder
     #CDF_FilePath = "D:/Sajjad_Work/Projects/Project_PowerEdu/PowerEdu/data/IEEE_14/IEEE_14_Data.txt"
    #  @show pwd()
    CDF_FilePath = "data/IEEE_14/IEEE_14_Data.txt"

     Ybus_Taps_Indicator_User = false

     NR_Type_User = 1

     Tolerance_User = 0.0001

     Tol_Num_User = 0

     SortValue_User = true
     
     BusSwitching_User = false
     
    #  PowerSystemsAnalysis.PowerFlow_MainFunction(CDF_FilePath, Ybus_Taps_Indicator, NR_Type, Tolerance, Tol_Num)
     PowerSystemsAnalysis.PowerFlow_MainFunction(CDF_FilePath,  Ybus_Taps_Indicator=Ybus_Taps_Indicator_User, NR_Type=NR_Type_User, Tolerance=Tolerance_User, Tol_Num=Tol_Num_User,SortValue=SortValue_User, BusSwitching=BusSwitching_User)

elseif (Debug_Indicator == 5) # Full Continuation Power Flow

    # Debugging YBus Builder
    #CDF_FilePath = "D:/Sajjad_Work/Projects/Project_PowerEdu/PowerEdu/data/IEEE_14/IEEE_14_Data.txt"
    #  @show pwd()
    CDF_FilePath = "data/IEEE_14/IEEE_14_Data.txt"

    PQ_V_Curve_Tuple = (10, 'P')

    Ybus_Taps_Indicator = false

    StepSize_Vector_CPF = [0.1 0.025]

    Tolerance_NR = 0.1

    Tol_Num = 0

    PostCriticalPoint_Counter_Input = 10

    SortValue = true
    
    BusSwitching = false
    
    #  PowerSystemsAnalysis.PowerFlow_MainFunction(CDF_FilePath, Ybus_Taps_Indicator, NR_Type, Tolerance, Tol_Num)
    PowerSystemsAnalysis.ContinuationPowerFlow_MainFunction(CDF_FilePath, PQ_V_Curve_Tuple, Ybus_Taps_Indicator, StepSize_Vector_CPF, Tolerance_NR, Tol_Num, PostCriticalPoint_Counter_Input, SortValue, BusSwitching)

elseif (Debug_Indicator == 6) # Full State Estimation


     # Debugging YBus Builder
     #CDF_FilePath = "D:/Sajjad_Work/Projects/Project_PowerEdu/PowerEdu/data/IEEE_14/IEEE_14_Data.txt"
    #  @show pwd()
    CDF_FilePath = "data/IEEE_14/IEEE_14_Data.txt"

    Measurement_Error_Variance = [0.015, 0.02, 0.02]

    Bad_Bus_Measurement_Input= nothing
    Bad_Branch_Measurement_Input = nothing

    # Bad_Bus_Measurement_Input = [6 1 1 1; 7 1 1 1]    
    
    # Bad_Branch_Measurement_Input = [1 2 1 1 1 1; 2 3 1 1 1 1]

    Tolerance_SE = 0.01

    alpha = 0.01

    Ybus_Taps_Indicator = false

    Tolerance_NR = 0.001

    Tol_Num = 0

    SortValue = true
    
    BusSwitching = false
     
    #  PowerSystemsAnalysis.PowerFlow_MainFunction(CDF_FilePath, Ybus_Taps_Indicator, NR_Type, Tolerance, Tol_Num)
     PowerSystemsAnalysis.PowerSystem_StateEstimation_MainFunction(CDF_FilePath, Measurement_Error_Variance, Bad_Bus_Measurement_Input, Bad_Branch_Measurement_Input, Tolerance_SE, alpha, Ybus_Taps_Indicator, Tolerance_NR, Tol_Num, SortValue, BusSwitching)

elseif (Debug_Indicator == 7)  # OPF: Economic Dispatch

    GeneratorCostCurve_Array = [[0.0625,1,0],[0.0125,1,0],[0.0250,1,0]] 
    
    Load_Demand  = 952

    PowerSystemsAnalysis.PowerSystem_EconomicDispatch_MainFunction(GeneratorCostCurve_Array, Load_Demand)

elseif (Debug_Indicator == 8)  # OPF: Full OPF



end