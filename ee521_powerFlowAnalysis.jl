using Pkg
Pkg.instantiate()
using CSV
using DelimitedFiles
using DataFrames


include("functions/Foo.jl")
include("src/IEEE_CDF_Parser.jl")


# folderInput = "rawData/";
folderInput = "data/"
# folder_processedData = "processedData/";
# systemName = "ieee14";
systemName = "IEEE_14"
fileType_CDFFile = ".txt"

filename_CDFFile = folderInput*systemName*"_Data"*fileType_CDFFile
CDF_DF_List = CDF_Parser(filename_CDFFile);
CDF_DF_List_pu = CDF_pu_Converter(CDF_DF_List)
busData_pu = CDF_DF_List_pu[2]
# vscodedisplay(busData_pu)
branchData_pu = CDF_DF_List_pu[3]


ybus, BMatrix, b, A, branchNames, E   = Foo.ybusGenerator(busData_pu, branchData_pu);
vscodedisplay(ybus)
ybus, BMatrix, b, A, branchNames, E   = Foo.ybusGenerator(busData, branchData, saveTables=true, systemName=systemName);
ybusByTypes= Foo.ybusGenerator(busData, branchData, sortBy="busTypes", verbose=true, saveTables=true, saveLocation=folder_processedData, systemName=systemName)   