using Revise
using Pkg
using CSV
using DataFrames
using LinearAlgebra
using Cthulhu
using BenchmarkTools
using Plots

include("src/Helper_Functions.jl");
include("src/Ybus_Builder.jl");
include("src/IEEE_CDF_Parser.jl");
include("src/SparseTechniques_Functions.jl");
include("src/Jacobian_Builder.jl");

systemName = "IEEE_14";
# systemName = "IEEE_118"
# systemName = "IEEE_30";
# systemName = "IEEE_57"

folderInput = "data/";
folder_processedData = "processedData/";


createFolderIfNotExisting(systemName, folder_processedData)

fileType_CDFFile = ".txt";
filename_CDFFile = folderInput*systemName*"/"*systemName*"_Data"*fileType_CDFFile
CDF_DF_List = CDF_Parser(filename_CDFFile);
CDF_DF_List_pu = CDF_pu_Converter(CDF_DF_List);

systemName = extractSystemName(CDF_DF_List)

busData_pu = CDF_DF_List_pu[2];
branchData_pu = CDF_DF_List_pu[3];

sparYBus = constructSparseYBus(CDF_DF_List_pu);
ybusSpar2Full = spar2Full(sparYBus, readMethod="col-wise");

ybusResults = ybusGenerator(CDF_DF_List_pu, saveTables=true);
ybusDense = ybusResults.ybus;
E = ybusResults.E ;

@test ybusSpar2Full ≈ ybusDense

powSysData = initializeVectors_pu(CDF_DF_List_pu);
PSpecified = powSysData.PSpecified;
QSpecified = powSysData.QSpecified;
V = powSysData.V;
delta = powSysData.delta;
lSlack = powSysData.listOfSlackBuses;
lPV = powSysData.listOfPVBuses;
lPQ = powSysData.listOfPQBuses;
lNonSlack = powSysData.listOfNonSlackBuses;
nPV = powSysData.nPV
nPQ = powSysData.nPQ
nSlack = powSysData.nSlack
nNonSlack = nPV+nPQ
N = nSlack + nNonSlack
# ATestDense = [10 1 0 3 0 0 0 5 0 0
# 2 9 0 0 0 0 5 0 0 2;
# 0 0 21 5 7 0 0 0 0 4;
# 4 0 1 18 8 0 0 0 0 0;
# 0 0 4 7 25 4 1 0 0 2;
# 0 0 0 0 3 14 9 0 0 0;
# 0 1 4 0 2 3 12 1 1 0;
# 1 0 5 0 0 0 5 10 0 0;
# 0 0 0 0 0 0 6 0 20 0;
# 0 2 3 0 4 0 0 0 0 35];

# ADense = [1 3 4 8; 2 1 2 3; 4 3 5 8; 9 2 7 4];
# A = sparmat(ADense);
# qluA = sparLU(A);
# QA = qluA.Q;
# α = qluA.α;
# QASpar2Full = spar2Full(QA);
# b = ones(Float64, 4);
# y, β▶ = sparForwardSolve(QA, b, verbose=false);
# x, β◀ = sparBackwardSolve(QA, y, verbose=false);



# x, numOperations, α, β = solveUsingSparseLU(A, b)

ΔP, ΔQ = computeMismatchesViaSparseYBus(PSpecified, QSpecified, V, delta, sparYBus);
mismatch = vcat(ΔP[lNonSlack], ΔQ[lPQ])

P = PSpecified - ΔP;
Q = QSpecified - ΔQ;

sparJ = constructSparseJacobian(CDF_DF_List_pu, P, Q, V, delta, sparYBus);
# JSpar2Full = spar2Full(sparJ);

qluJ = sparLU(sparJ);
QJ, LJ, UJ, αJ, fill_insJ = qluJ.Q, qluJ.L, qluJ.U, qluJ.α, qluJ.fills;

# QJSpar2Full = spar2Full(QJ);


correction = solveUsingSparseLU(QJ, mismatch).x
Δδ = correction[1:nNonSlack]
ΔVbyV = correction[N:N+nPQ-1]
delta[lNonSlack] = delta[lNonSlack] + Δδ
V[lPQ] .*= (1 .+ (ΔVbyV))