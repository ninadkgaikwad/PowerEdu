using Pkg
using CSV
using DataFrames


include("src/Ybus_Builder.jl")
include("src/IEEE_CDF_Parser.jl")
include("src/SparseTechniques_Functions.jl")

folderInput = "data/"
folder_processedData = "processedData/";
systemName = "IEEE_14";
# systemName = "IEEE_30";

# Ybus_Builder.createFolderIfNotExisting(systemName, folder_processedData)
createFolderIfNotExisting(systemName, folder_processedData)

fileType_CDFFile = ".txt";
filename_CDFFile = folderInput*systemName*"/"*systemName*"_Data"*fileType_CDFFile
CDF_DF_List = CDF_Parser(filename_CDFFile, saveTables = true, saveLocation = folder_processedData);

CDF_DF_List_pu = CDF_pu_Converter(CDF_DF_List, saveTables = true, saveLocation = folder_processedData);
CDF_DF_List_pu = CDF_pu_Converter(CDF_DF_List);

systemName = extractSystemName(CDF_DF_List)

busData_pu = CDF_DF_List_pu[2];
branchData_pu = CDF_DF_List_pu[3];

NYBus, nnzYBus = constructSparseYBus(CDF_DF_List_pu)
# vscodedisplay(NYBus)
vscodedisplay(nnzYBus)

ybus, BMatrix, b, A, branchNames, E   = ybusGenerator(CDF_DF_List_pu, saveTables=true);
vscodedisplay(ybus)
ybus_ByTypes, BMatrix_ByTypes, b_ByTypes, A_ByTypes, branchNames_ByTypes, E_ByTypes = ybusGenerator(CDF_DF_List_pu, sortBy="busTypes", saveTables=true, saveLocation=folder_processedData);