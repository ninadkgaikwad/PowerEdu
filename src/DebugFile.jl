include("PowerSystemsAnalysis.jl")

using .PowerSystemsAnalysis

# Which Code to Debug
Debug_Indicator = 2

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

    CDF_FilePath = "C:/Users/ninad/Dropbox (Personal)/NinadGaikwad_PhD/Gaikwad_Research/Gaikwad_Research_Work/PowerSystemsAnalysis/data/IEEE_14_Data.txt"

    CDF_DF_List = PowerSystemsAnalysis.CDF_Parser(CDF_FilePath)

    CDF_DF_List_pu = CDF_pu_Converter(CDF_DF_List)

end