using Revise
using Pkg
using CSV
using DataFrames
using LinearAlgebra
using Cthulhu
using BenchmarkTools
using Plots
using Statistics
using LaTeXStrings
using Symbolics
using ForwarDiff
using NLsolve

include("src/Helper_Functions.jl");
include("src/Ybus_Builder.jl");
include("src/IEEE_CDF_Parser.jl");
include("src/SparseTechniques_Functions.jl");
include("src/Jacobian_Builder.jl");
include("src/LU_Factorization.jl");
include("src/BasicPowerFlow_Functions.jl")

folderInput = "data/";
folder_processedData = "processedData/";

systemName = "IEEE_14";
# systemName = "IEEE_118"
# systemName = "IEEE_30";
# systemName = "IEEE_57"

createFolderIfNotExisting(systemName, folder_processedData)

fileType_CDFFile = ".txt";
filename_CDFFile = folderInput*systemName*"/"*systemName*"_Data"*fileType_CDFFile
CDF_DF_List = CDF_Parser(filename_CDFFile);
dfpu = CDF_pu_Converter(CDF_DF_List);

results = solveForPowerFlow_Sparse(dfpu, verbose=false)

plotBuswiseDifferences(dfpu, results, savePlots=false)

solutions = solveForEconomicDispatch(dfpu, x, f, h)

# @variables P₁ P₂ λ₁;
# x = [P₁, P₂, λ₁];
# f₁ = 0.004P₁^2 + 8P₁;
# f₂ = 0.0048P₂^2 + 6.4P₂;
h₁ = P₁ + P₂ - Pₗ;
# f = f₁ + f₂ - λ₁*h₁;


# P₁′, P₂′, λ₁′ = solutions;


