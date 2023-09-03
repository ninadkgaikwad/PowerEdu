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
using ForwardDiff
using NLsolve

include("src/Helper_Functions.jl");
include("src/Ybus_Builder.jl");
include("src/IEEE_CDF_Parser.jl");
include("src/SparseTechniques_Functions.jl");
include("src/Jacobian_Builder.jl");
include("src/LU_Factorization.jl");
include("src/BasicPowerFlow_Functions.jl");
include("src/OptimalPowerFlow_Functions.jl");

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

Pₗ₁ = 259;
@variables P₁ P₂;
x = [P₁, P₂];
f₁ = 0.004P₁^2 + 8P₁;
f₂ = 0.0048P₂^2 + 6.4P₂;
h₁ = P₁ + P₂ - Pₗ₁;
# h₂ = h₁;
h = [h₁];
# h = [h₁, h₂];
f = f₁ + f₂;
solutions1 = solveForEconomicDispatch(dfpu, x, f, h, verbose=false);
P₁′, P₂′, λ₁′ = solutions1;

MVA_B = getBaseMVA(dfpu);
Pₗ₂ = sum(results.P[1:2])*MVA_B;
h₂ = P₁ + P₂ - Pₗ₂;
h = [h₂];
solutions2 = solveForEconomicDispatch(dfpu, x, f, h, verbose=false);    
P₁′, P₂′, λ₁′ = solutions2;


