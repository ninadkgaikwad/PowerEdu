# systemName = extractSystemName(CDF_DF_List)

# busData_pu = CDF_DF_List_pu[2];
# branchData_pu = CDF_DF_List_pu[3];

# sparYBus = constructSparseYBus(CDF_DF_List_pu);
# ybusSpar2Full = spar2Full(sparYBus, readMethod="col-wise");

# ybusResults = ybusGenerator(CDF_DF_List_pu, saveTables=false);
# ybusDense = ybusResults.ybus;
# E = ybusResults.E ;

# @test ybusSpar2Full ≈ ybusDense

# powSysData = initializeVectors_pu(CDF_DF_List_pu);
# PSpecified = powSysData.PSpecified;
# QSpecified = powSysData.QSpecified;
# V = powSysData.V;
# delta = powSysData.delta;
# lSlack = powSysData.listOfSlackBuses;
# lPV = powSysData.listOfPVBuses;
# lPQ = powSysData.listOfPQBuses;
# lNonSlack = powSysData.listOfNonSlackBuses;
# nPV = powSysData.nPV;
# nPQ = powSysData.nPQ;
# nSlack = powSysData.nSlack;
# nNonSlack = nPV+nPQ;
# N = nSlack + nNonSlack;
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

# ADense = [2 4 2; 1 1 2; -1 0 2]
# A = sparmat(ADense);
# QA = sparLU(A).Q;
# QASpar2Full = spar2Full(QA)

# ADense = [1 3 4 8; 2 1 2 3; 4 3 5 8; 9 2 7 4]

# # ADense = [3 -7 -2 2; -3 5 1 0; 6 -4 0 -5; -9 5 -5 12]
# # ADense = [1 -1 2; 3 -1 7; 2 -4 5]
# A = sparmat(ADense);
# qluA = sparLU(A, method="Crout", verbose=true);
# LASpar2Full = spar2Full(qluA.L)
# UASpar2Full = spar2Full(qluA.U)
# LASpar2Full*UASpar2Full - ADense
# QA = qluA.Q;
# QASpar2Full = spar2Full(QA)

# P, L, U = Compute_PLU(ADense, 1e-3);
# L
# U
# Q = L + U - I
# x, numOperations, α, β = solveUsingSparseLU(A, b)



# JSpar2Full = spar2Full(sparJ);



# QJSpar2Full = spar2Full(QJ);
# vscodedisplay(QJSpar2Full)