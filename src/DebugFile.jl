include("PowerSystemsAnalysis.jl")

using .PowerSystemsAnalysis

using DataFrames

using DelimitedFiles

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

    Ybus_Taps_Indicator = 2

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

elseif (Debug_Indicator == 4) # Create Initial Solution Vector

     # Debugging YBus Builder
     CDF_FilePath = "C:/Users/ninad/Dropbox (Personal)/NinadGaikwad_PhD/Gaikwad_Research/Gaikwad_Research_Work/PowerSystemsAnalysis/data/IEEE14/IEEE_14_Data.txt"

     Ybus_Taps_Indicator = 2
 
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

     # Creating Initial Solution Vector for Power Flow
    Initial_SolutionVector_NR = PowerSystemsAnalysis.Create_Initial_SolutionVector_NR(CDF_DF_List_pu)

end