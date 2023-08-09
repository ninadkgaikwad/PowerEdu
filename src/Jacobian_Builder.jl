# Jacobian_Builder.jl

"""
    Create_Jacobian_NR(CDF_DF_List_pu, Ybus, SolutionVector_V,
    SolutionVector_Delta, PQ_BusArray)

Creates Jacobian Matrix for power system network for full Newton-Raphson method.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
- 'Ybus': A complex array of Ybus elements ordered according to bus
type: Slack->PQ->PV.
- 'SolutionVector_V': Voltage at each bus ordered according to bus
type: Slack->PQ->PV.
- 'SolutionVector_Delta': Angle at each bus ordered according to bus
type: Slack->PQ->PV.
- 'PQ_BusArray': An array (N*2) for P and Q vectors ordered according to bus
type: Slack->PQ->PV.
- 'NR_Type': 1 -> Full Newton-Raphson, 2-> Decoupled Newton-Raphson,
3 -> Fast Decoupled Newton-Raphson
- 'ContinuationPowerFlow_Indicator': 0 - Not Continuation Power Flow,
1 - Yes Continuation Power Flow
'''
'''
# Output
- 'Jacobian_Matrix': Jacobian Matrix for full Newton-Raphson method: PQ->PV.
'''
"""
function Create_Jacobian_NR(CDF_DF_List_pu, Ybus, SolutionVector_V, SolutionVector_Delta, PQ_BusArray, NR_Type, ContinuationPowerFlow_Indicator)

    # Getting required data from CDF_DF_List
    BusDataCard_DF = CDF_DF_List_pu[2]

    # Number of Buses
    N_Bus = nrow(BusDataCard_DF)
    if (ContinuationPowerFlow_Indicator == 0)

        N_PQ_Bus = nrow(filter(row -> ((row.Type == 0) || (row.Type == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type == 2), BusDataCard_DF))

    elseif (ContinuationPowerFlow_Indicator == 1)

        N_PQ_Bus = nrow(filter(row -> ((row.Type_Original == 0) || (row.Type_Original == 1)), BusDataCard_DF))
        N_PV_Bus = nrow(filter(row -> (row.Type_Original == 2), BusDataCard_DF))

    end

    # Creating Jacobian Matrix based on NR_Type
    if (NR_Type == 1) # Full Newton-Raphson

        # Computing Sizes for Jacobian Submatrices
        J_11_RowNum = N_PQ_Bus+N_PV_Bus
        J_11_ColNum = N_PQ_Bus+N_PV_Bus

        J_12_RowNum = N_PQ_Bus+N_PV_Bus
        J_12_ColNum = N_PQ_Bus

        J_21_RowNum = N_PQ_Bus
        J_21_ColNum = N_PQ_Bus+N_PV_Bus

        J_22_RowNum = N_PQ_Bus
        J_22_ColNum = N_PQ_Bus

        # Initializing the Jacobian Submatrices
        J_11 = zeros(J_11_RowNum,J_11_ColNum)
        J_12 = zeros(J_12_RowNum,J_12_ColNum)
        J_21 = zeros(J_21_RowNum,J_21_ColNum)
        J_22 = zeros(J_22_RowNum,J_22_ColNum)

        # Computing J_11
        for ii in 1:J_11_RowNum # Through Rows

            for jj in 1:J_11_ColNum # Through Columns

                if (ii == jj) # Diagonal Term

                    J_11[ii,jj] = -(PQ_BusArray[ii+1,2]) - (SolutionVector_V[ii+1]^(2)*imag(Ybus[ii+1,ii+1]))

                else # Off-Diagonal Term

                    J_11[ii,jj] = -(SolutionVector_V[ii+1]*SolutionVector_V[jj+1]*abs(Ybus[ii+1,jj+1]) *sin(angle(Ybus[ii+1,jj+1]) + deg2rad(SolutionVector_Delta[jj+1]) - deg2rad(SolutionVector_Delta[ii+1])))

                end

            end

        end

        # Computing J_12
        for ii in 1:J_12_RowNum # Through Rows

            for jj in 1:J_12_ColNum # Through Columns

                if (ii == jj) # Diagonal Term

                    J_12[ii,jj] = (PQ_BusArray[ii+1,1]) + (SolutionVector_V[ii+1]^(2)*real(Ybus[ii+1,ii+1]))

                else # Off-Diagonal Term

                    J_12[ii,jj] = (SolutionVector_V[ii+1]*SolutionVector_V[jj+1]*abs(Ybus[ii+1,jj+1]) *cos(angle(Ybus[ii+1,jj+1]) + deg2rad(SolutionVector_Delta[jj+1]) - deg2rad(SolutionVector_Delta[ii+1])))

                end

            end

        end

        # Computing J_21
        for ii in 1:J_21_RowNum # Through Rows

            for jj in 1:J_21_ColNum # Through Columns

                if (ii == jj) # Diagonal Term

                    J_21[ii,jj] = (PQ_BusArray[ii+1,1]) - (SolutionVector_V[ii+1]^(2)*real(Ybus[ii+1,ii+1]))

                else # Off-Diagonal Term

                    J_21[ii,jj] = -(SolutionVector_V[ii+1]*SolutionVector_V[jj+1]*abs(Ybus[ii+1,jj+1]) *cos(angle(Ybus[ii+1,jj+1]) + deg2rad(SolutionVector_Delta[jj+1]) - deg2rad(SolutionVector_Delta[ii+1])))

                end

            end

        end

        # Computing J_22
        for ii in 1:J_22_RowNum # Through Rows

            for jj in 1:J_22_ColNum # Through Columns

                if (ii == jj) # Diagonal Term

                    J_22[ii,jj] = (PQ_BusArray[ii+1,2]) - (SolutionVector_V[ii+1]^(2)*imag(Ybus[ii+1,ii+1]))

                else # Off-Diagonal Term

                    J_22[ii,jj] = -(SolutionVector_V[ii+1]*SolutionVector_V[jj+1]*abs(Ybus[ii+1,jj+1]) *sin(angle(Ybus[ii+1,jj+1]) + deg2rad(SolutionVector_Delta[jj+1]) - deg2rad(SolutionVector_Delta[ii+1])))

                end

            end

        end

        # Creating Jacobian_NR
        Jacobian_NR1 = hcat(J_11,J_12)
        Jacobian_NR2 = hcat(J_21,J_22)

        Jacobian_NR = vcat(Jacobian_NR1,Jacobian_NR2)

    elseif (NR_Type == 2) # Decoupled Newton-Raphson

            # Computing Sizes for Jacobian Submatrices
            J_11_RowNum = N_PQ_Bus+N_PV_Bus
            J_11_ColNum = N_PQ_Bus+N_PV_Bus

            J_22_RowNum = N_PQ_Bus
            J_22_ColNum = N_PQ_Bus

            # Initializing the Jacobian Submatrices
            J_11 = zeros(J_11_RowNum,J_11_ColNum)
            J_22 = zeros(J_22_RowNum,J_22_ColNum)

            # Computing J_11
            for ii in 1:J_11_RowNum # Through Rows

                for jj in 1:J_11_ColNum # Through Columns

                    if (ii == jj) # Diagonal Term

                        J_11[ii,jj] = -(PQ_BusArray[ii+1,2]) - (SolutionVector_V[ii+1]^(2)*imag(Ybus[ii+1,ii+1]))

                    else # Off-Diagonal Term

                        J_11[ii,jj] = -(SolutionVector_V[ii+1]*SolutionVector_V[jj+1]*abs(Ybus[ii+1,jj+1]) *sin(angle(Ybus[ii+1,jj+1]) + deg2rad(SolutionVector_Delta[jj+1]) - deg2rad(SolutionVector_Delta[ii+1])))

                    end

                end

            end

            # Computing J_22
            for ii in 1:J_22_RowNum # Through Rows

                for jj in 1:J_22_ColNum # Through Columns

                    if (ii == jj) # Diagonal Term

                        J_22[ii,jj] = (PQ_BusArray[ii+1,2]) - (SolutionVector_V[ii+1]^(2)*imag(Ybus[ii+1,ii+1]))

                    else # Off-Diagonal Term

                        J_22[ii,jj] = -(SolutionVector_V[ii+1]*SolutionVector_V[jj+1]*abs(Ybus[ii+1,jj+1]) *sin(angle(Ybus[ii+1,jj+1]) + deg2rad(SolutionVector_Delta[jj+1]) - deg2rad(SolutionVector_Delta[ii+1])))

                    end

                end

            end

            # Creating Jacobian_NR
            Jacobian_NR = [J_11, J_22]

    elseif (NR_Type == 3) # Fast Decoupled Newton-Raphson

        # Computing Sizes for Jacobian Submatrices
        J_11_RowNum = N_PQ_Bus+N_PV_Bus
        J_11_ColNum = N_PQ_Bus+N_PV_Bus

        J_22_RowNum = N_PQ_Bus
        J_22_ColNum = N_PQ_Bus

        # Initializing the Jacobian Submatrices
        J_11 = zeros(J_11_RowNum,J_11_ColNum)
        J_22 = zeros(J_22_RowNum,J_22_ColNum)

        # Computing J_11
        for ii in 1:J_11_RowNum # Through Rows

            for jj in 1:J_11_ColNum # Through Columns

                J_11 = -imag(Ybus[ii+1,jj+1])

            end

        end

        # Computing J_22
        for ii in 1:J_22_RowNum # Through Rows

            for jj in 1:J_22_ColNum # Through Columns

                J_22 = -imag(Ybus[ii+1,jj+1])

            end

        end

        # Creating Jacobian_NR
        Jacobian_NR = [J_11, J_22]

    end

    return Jacobian_NR

end

"""
    Create_Jacobian_CPF_Predict(CDF_DF_List_pu, Ybus, SolutionVector_CPF,
    NR_Type, K_Vector, Index_CPF)

Creates Jacobian Matrix for power system network for Continuation Power Flow
predict step.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
- 'Ybus': A complex array of Ybus elements ordered according to bus
type: Slack->PQ->PV.
- 'SolutionVector_CPF': Power Flow Solution Voltage and Angle at each
bus ordered according to bus type: PQ->PV->Lambda.
- 'NR_Type': 1 -> Full Newton-Raphson, 2-> Decoupled Newton-Raphson,
3 -> Fast Decoupled Newton-Raphson
- 'K_Vector': K Vector of Continuation Power Flow
- 'Index_CPF': Continuation Parameter Index
'''
'''
# Output
- 'Jacobian_CPF_Predict': Jacobian Matrix for Predict Step of Continuation Power
Flow: PQ->PV->Lambda.
'''
"""
function Create_Jacobian_CPF_Predict(CDF_DF_List_pu, Ybus, SolutionVector_CPF, NR_Type, K_Vector, Index_CPF)

    # Creating Solution vectors for V Delta
    SolutionVector_V, SolutionVector_Delta = Create_SolutionVector_VDelta_CPF(CDF_DF_List_pu, SolutionVector_CPF)

    # Computing PQ Bus Array
    PQ_BusArray = Compute_PQ_BusArray(Ybus, SolutionVector_V, SolutionVector_Delta)

    # Computing NR Jacobian
    Jacobian_NR_Predict = Create_Jacobian_NR(CDF_DF_List_pu, Ybus, SolutionVector_V, SolutionVector_Delta, PQ_BusArray, NR_Type, 1)

    # Augmenting NR Jacobian with K Vector
    Jacobian_CPF_Predict = hcat(Jacobian_NR_Predict, K_Vector)

    # Creating Last row of CPF Jacobian
    Jacobian_CPF_Predict_LastRow = zeros(1,size(Jacobian_CPF_Predict)[2])

    Jacobian_CPF_Predict_LastRow[1,Index_CPF] = 1

    # Augmenting Last Row to CPF Jacobian
    Jacobian_CPF_Predict = vcat(Jacobian_CPF_Predict, Jacobian_CPF_Predict_LastRow)

    return Jacobian_CPF_Predict

end

"""
    Create_Jacobian_CPF_Correct(CDF_DF_List_pu, Ybus, SolutionVector_CPF,
    NR_Type, K_Vector, Index_CPF)

Creates Jacobian Matrix for power system network for Continuation Power Flow
for the correct step.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
- 'Ybus': A complex array of Ybus elements ordered according to bus
type: Slack->PQ->PV.
- 'SolutionVector_CPF': Power Flow Solution Voltage and Angle at each
bus ordered according to bus type: PQ->PV->Lambda.
- 'NR_Type': 1 -> Full Newton-Raphson, 2-> Decoupled Newton-Raphson,
3 -> Fast Decoupled Newton-Raphson
- 'K_Vector': K Vector of Continuation Power Flow
- 'Index_CPF': Continuation Parameter Index
'''
'''
# Output
- 'Jacobian_CPF_Correct': Jacobian Matrix for Predict Step of Continuation Power
Flow: PQ->PV->Lambda.
'''
"""
function Create_Jacobian_CPF_Correct(CDF_DF_List_pu, Ybus, CPF_Predictor_Vector, NR_Type, K_Vector, Index_CPF)

    # Creating Solution vectors for V Delta
    SolutionVector_V, SolutionVector_Delta = Create_SolutionVector_VDelta_CPF(CDF_DF_List_pu, CPF_Predictor_Vector)

    # Computing PQ Bus Array
    PQ_BusArray = Compute_PQ_BusArray(Ybus, SolutionVector_V, SolutionVector_Delta)

    # Computing NR Jacobian
    Jacobian_NR_Correct = Create_Jacobian_NR(CDF_DF_List_pu, Ybus, SolutionVector_V, SolutionVector_Delta, PQ_BusArray, NR_Type, 1)

    # Augmenting NR Jacobian with K Vector
    Jacobian_CPF_Correct = hcat(Jacobian_NR_Predict, K_Vector)

    # Creating Last row of CPF Jacobian
    Jacobian_CPF_Correct_LastRow = zeros(1,size(Jacobian_CPF_Predict)[2])

    Jacobian_CPF_Correct_LastRow[1,Index_CPF] = 1

    # Augmenting Last Row to CPF Jacobian
    Jacobian_CPF_Correct = vcat(Jacobian_CPF_Predict, Jacobian_CPF_Predict_LastRow)

    return Jacobian_CPF_Correct

end

function constructJacobian(CDF_DF_List_pu::Vector{DataFrame},
    P::Vector{Float64}, Q::Vector{Float64}, V::Vector{Float64}, 
    delta::Vector{Float64}, ybus::Matrix{ComplexF64};
    E::Vector{Vector{Int64}},
    combinationOrder::String="hcat-then-vcat",
    constructionType::String="byYBusElements",
    verbose::Bool = false,
    save::Bool = true,
    saveLocation::String = "processedData/",
    systemName::String = "systemNOTSpecified",
    fileExtension::String = ".csv")::Matrix{Float64}

    if !isdefined(Main, :E)
        E = [[i for i in 1:N] for _ in 1:N]
    end

    powSysData = initializeVectors_pu(CDF_DF_List_pu)
    nPV = powSysData.nPV
    nPQ = powSysData.nPQ
    nSlack = powSysData.nSlack
    lPV = powSysData.listOfPVBuses
    lPQ = powSysData.listOfPQBuses
    lNonSlack = powSysData.listOfNonSlackBuses
    
    busData_pu = CDF_DF_List_pu[2]
    N = size(busData_pu, 1);
    J11 = constructJacobianSubMatrix(CDF_DF_List_pu, P, Q, V, delta, ybus, E, type="J11")
    J12 = constructJacobianSubMatrix(CDF_DF_List_pu, P, Q, V, delta, ybus, E, type="J12")
    J21 = constructJacobianSubMatrix(CDF_DF_List_pu, P, Q, V, delta, ybus, E, type="J21")
    J22 = constructJacobianSubMatrix(CDF_DF_List_pu, P, Q, V, delta, ybus, E, type="J22")

    if combinationOrder == "hcat-then-vcat"
        JTop = hcat(J11, J12)
        JBottom = hcat(J21, J22)
        J = vcat(JTop, JBottom)
    elseif combinationOrder == "vcat-then-hcat"
        JLeft = vcat(J11, J21)
        JRight = vcat(J12, J22)
        J = hcat(JLeft, JRight)
    else
        error("Unknown combination order.")
    end

    return J
    
end

function constructJacobianSubMatrix(CDF_DF_List_pu::Vector{DataFrame},
    P::Vector{Float64},
    Q::Vector{Float64},
    V::Vector{Float64},
    delta::Vector{Float64},
    ybus::Matrix{ComplexF64},
    E::Vector{Vector{Int64}};
    type::String="J11",
    constructionType::String="byYBusElements",
    verbose::Bool=false)

    busData_pu = CDF_DF_List_pu[2]
    N = size(busData_pu, 1)

    busPositions = busPositionsForJacobian(CDF_DF_List_pu)

    powSysData = initializeVectors_pu(CDF_DF_List_pu)
    nPV = powSysData.nPV
    nPQ = powSysData.nPQ
    lPV = powSysData.listOfPVBuses
    lPQ = powSysData.listOfPQBuses
    lNonSlack = powSysData.listOfNonSlackBuses

    if type == "J11"
        JSub = zeros(Float64, N-1, N-1)
    elseif type == "J12"
        JSub = zeros(Float64, N-1, nPQ)
    elseif type == "J21"
        JSub = zeros(Float64, nPQ, N-1)
    elseif type == "J22"
        JSub = zeros(Float64, nPQ, nPQ)
    else
        error("Unknown Jacobian Sub-matrix.")
    end

    busPositions = busPositionsForJacobian(CDF_DF_List_pu)

    if constructionType == "byYBusElements"
        for i = 1:N
            k = 1
            for k = E[i]
                Y_ik = ybus[i, k]
                if Y_ik != 0

                    if type == "J11"
                        row = busPositions.nonSlackBusIdx[i]
                        col = busPositions.nonSlackBusIdx[k]
                    elseif type == "J12"
                        row = busPositions.nonSlackBusIdx[i]
                        col = busPositions.PQBusIdx[k]
                    elseif type == "J21"
                        row = busPositions.PQBusIdx[i]
                        col = busPositions.nonSlackBusIdx[k]
                    elseif type == "J22"
                        row = busPositions.PQBusIdx[i]
                        col = busPositions.PQBusIdx[k]
                    else
                        error("Unknown Jacobian Submatrix type")
                    end

                    if row != -1 && col != -1

                        if i == k
                            B_ii = imag(Y_ik)
                            G_ii = real(Y_ik)
                            Vi_squared = V[i]^2
                            j11 = -Q[i] - B_ii*Vi_squared
                            j12 = P[i] + G_ii*Vi_squared
                            j21 = P[i] - G_ii*Vi_squared
                            j22 = Q[i] - B_ii*Vi_squared
                        else
                            YVkVi = abs(Y_ik)*V[k]*V[i]
                            gamma_deltak_delta_i = angle(Y_ik) + delta[k] - delta[i]
                            j11 = -YVkVi*sin(gamma_deltak_delta_i)
                            j12 = YVkVi*cos(gamma_deltak_delta_i)
                            j21 = -j12
                            j22 = j11
                        end
                        
                        if type == "J11"
                            JSub[row, col] = j11
                        elseif type == "J12"
                            JSub[row, col] = j12
                        elseif type == "J21"
                            JSub[row, col] = j21
                        elseif type == "J22"
                            JSub[row, col] = j22
                        else
                            error("Unknown Jacobian Submatrix type.")
                        end

                    end

                end
            end
        end
    end

    return JSub
end