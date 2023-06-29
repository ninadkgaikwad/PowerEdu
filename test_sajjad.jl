using CSV
using DataFrames

include("src/IEEE_CDF_Parser.jl")

## Debugging CDF Parser

CDF_FilePath = "D:/Sajjad_Work/Projects/Project_PowerEdu/PowerEdu/data/IEEE_14/IEEE_14_Data.txt"

CDF_DF_List = CDF_Parser(CDF_FilePath)

CDF_DF_List_pu = CDF_pu_Converter(CDF_DF_List)