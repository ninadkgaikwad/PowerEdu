using Revise

using Pkg
using CSV
using DataFrames
using LinearAlgebra

include("src/Helper_Functions.jl")
include("src/Ybus_Builder.jl")
include("src/IEEE_CDF_Parser.jl")
include("src/SparseTechniques_Functions.jl")
include("src/Jacobian_Builder.jl")
folderInput = "data/"
folder_processedData = "processedData/";
# systemName = "IEEE_14";
# systemName = "IEEE_118"
systemName = "IEEE_30";

createFolderIfNotExisting(systemName, folder_processedData)

fileType_CDFFile = ".txt";
filename_CDFFile = folderInput*systemName*"/"*systemName*"_Data"*fileType_CDFFile
CDF_DF_List = CDF_Parser(filename_CDFFile, saveTables = true, saveLocation = folder_processedData);
busData = CDF_DF_List[2];
# vscodedisplay(busData)
# CDF_DF_List_pu = CDF_pu_Converter(CDF_DF_List, saveTables = true, saveLocation = folder_processedData);
CDF_DF_List_pu = CDF_pu_Converter(CDF_DF_List);

systemName = extractSystemName(CDF_DF_List)

busData_pu = CDF_DF_List_pu[2];
# vscodedisplay(busData_pu)
branchData_pu = CDF_DF_List_pu[3];

sparYBus = constructSparseYBus(CDF_DF_List_pu);
NYBus = sparYBus.NVec; MYBus = sparYBus.MVec; nnzYBus = sparYBus.nnzVec;

ybusFull = spar2Full(sparYBus, readMethod="col-wise")
# vscodedisplay(ybusFull)
ybus, BMatrix, b, A, branchNames, E   = ybusGenerator(CDF_DF_List_pu, saveTables=false);
# vscodedisplay(ybus)
# ybus_ByTypes, BMatrix_ByTypes, b_ByTypes, A_ByTypes, branchNames_ByTypes, E_ByTypes = ybusGenerator(CDF_DF_List_pu, sortBy="busTypes", saveTables=true, saveLocation=folder_processedData);

PSpecified, QSpecified, V, delta, listOfSlackBuses, listOfPVBuses, 
listOfPQBuses, listOfNonSlackBuses, nSlack, nPV, nPQ = initializeVectors_pu(CDF_DF_List_pu);

deltaP, deltaQ = computeMismatchesViaSparseYBus(PSpecified, QSpecified, V, delta, sparYBus);

P = PSpecified - deltaP;
Q = QSpecified - deltaQ;

sparJ = constructSparseJacobian(CDF_DF_List_pu, P, Q, V, delta, sparYBus);
JFull = real.(spar2Full(sparJ))::Matrix{Float64};

# A = [1 3 4 8; 2 1 2 3; 4 3 5 8; 9 2 7 4];

JRegular = constructJacobian(CDF_DF_List_pu, P, Q, V, delta, ybus, E=E);
# @vscodedisplay(JRegular)
@test JFull ≈ JRegular atol=1e-3

# ATestFull = [10 1 0 3 0 0 0 5 0 0
# 2 9 0 0 0 0 5 0 0 2;
# 0 0 21 5 7 0 0 0 0 4;
# 4 0 1 18 8 0 0 0 0 0;
# 0 0 4 7 25 4 1 0 0 2;
# 0 0 0 0 3 14 9 0 0 0;
# 0 1 4 0 2 3 12 1 1 0;
# 1 0 5 0 0 0 5 10 0 0;
# 0 0 0 0 0 0 6 0 20 0;
# 0 2 3 0 4 0 0 0 0 35];

qluJ = sparLU(sparJ);
QJ, LJ, UK, αJ, fill_insJ = qluJ.Q, qluJ.L, qluJ.U, qluJ.α, qluJ.fills;

αJ
fill_insJ
QJFull = real.(spar2Full(QJ))
# vscodedisplay(QJFull)
LJFull = real.(spar2Full(LJ));
UJFull = real.(spar2Full(UJ));

@test JFull ≈ JRegular atol=1e-3
@test LJFull*UJFull ≈ JRegular atol=1e-3
@test LJFull*UJFull ≈ JFull atol=1e-3