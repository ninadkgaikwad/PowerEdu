using Revise
using Pkg
using CSV
using DataFrames
using LinearAlgebra
using Cthulhu
using BenchmarkTools
using Plots
using Statistics

include("src/Helper_Functions.jl");
include("src/Ybus_Builder.jl");
include("src/IEEE_CDF_Parser.jl");
include("src/SparseTechniques_Functions.jl");
include("src/Jacobian_Builder.jl");
include("src/LU_Factorization.jl");

folderInput = "data/";
folder_processedData = "processedData/";

# systemName = "IEEE_14";
# systemName = "IEEE_118"
systemName = "IEEE_30";
# systemName = "IEEE_57"




createFolderIfNotExisting(systemName, folder_processedData)

fileType_CDFFile = ".txt";
filename_CDFFile = folderInput*systemName*"/"*systemName*"_Data"*fileType_CDFFile
CDF_DF_List = CDF_Parser(filename_CDFFile);
CDF_DF_List_pu = CDF_pu_Converter(CDF_DF_List);

results = solveForPowerFlow_Sparse(CDF_DF_List_pu, verbose=true)


