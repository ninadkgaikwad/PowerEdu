using Pkg
Pkg.instantiate()
using CSV
using DataFrames


include("src/Ybus_Builder.jl")
include("src/IEEE_CDF_Parser.jl")


folderInput = "data/"
folder_processedData = "processedData/";
systemName = "IEEE_14"
fileType_CDFFile = ".txt"

filename_CDFFile = folderInput*systemName*"_Data"*fileType_CDFFile
CDF_DF_List = IEEE_CDF_Parser.CDF_Parser(filename_CDFFile);

CDF_DF_List_pu = IEEE_CDF_Parser.CDF_pu_Converter(CDF_DF_List);

systemName = Ybus_Builder.extractSystemName(CDF_DF_List)
busData_pu = CDF_DF_List_pu[2];
# vscodedisplay(busData_pu)
branchData_pu = CDF_DF_List_pu[3];


ybus, BMatrix, b, A, branchNames, E   = Ybus_Builder.ybusGenerator(CDF_DF_List_pu);
ybus, BMatrix, b, A, branchNames, E   = Ybus_Builder.ybusGenerator(CDF_DF_List_pu, saveTables=true);
ybusByTypes, = Ybus_Builder.ybusGenerator(CDF_DF_List_pu, sortBy="busTypes", saveTables=true, saveLocation=folder_processedData);