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

                    J_12[ii,jj] = ((PQ_BusArray[ii+1,1]) + (SolutionVector_V[ii+1]^(2)*real(Ybus[ii+1,ii+1])))

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

                    J_22[ii,jj] = ((PQ_BusArray[ii+1,2]) - (SolutionVector_V[ii+1]^(2)*imag(Ybus[ii+1,ii+1])))

                else # Off-Diagonal Term

                    J_22[ii,jj] = -(SolutionVector_V[ii+1]*SolutionVector_V[jj+1]*abs(Ybus[ii+1,jj+1]) *sin(angle(Ybus[ii+1,jj+1]) + deg2rad(SolutionVector_Delta[jj+1]) - deg2rad(SolutionVector_Delta[ii+1])))

                end

            end

        end

        # Creating Jacobian_NR
        Jacobian_NR1 = hcat(J_11,J_12)
        Jacobian_NR2 = hcat(J_21,J_22)

        Jacobian_NR = vcat(Jacobian_NR1,Jacobian_NR2)

        # Addressing Machine Precision Problem
        for ii in 1:size(Jacobian_NR)[1]

            for jj in 1:size(Jacobian_NR)[2]

                if (abs(Jacobian_NR[ii,jj]) < 1e-12)

                    Jacobian_NR[ii,jj] = 0

                end

            end

        end

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

            # Addressing Machine Precision Problem
            for ii in 1:size(J_11)[1]

                for jj in 1:size(J_11)[2]

                    if (abs(J_11[ii,jj]) < 1e-12)

                        J_11[ii,jj] = 0

                    end

                end

            end

            for ii in 1:size(J_22)[1]

                for jj in 1:size(J_22)[2]

                    if (abs(J_22[ii,jj]) < 1e-12)

                        J_22[ii,jj] = 0

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

        # Addressing Machine Precision Problem
        for ii in 1:size(J_11)[1]

            for jj in 1:size(J_11)[2]

                if (abs(J_11[ii,jj]) < 1e-12)

                    J_11[ii,jj] = 0

                end

            end

        end

        for ii in 1:size(J_22)[1]

            for jj in 1:size(J_22)[2]

                if (abs(J_22[ii,jj]) < 1e-12)

                    J_22[ii,jj] = 0

                end

            end

        end

        # Creating Jacobian_NR
        Jacobian_NR = [J_11, J_22]

    end

    return Jacobian_NR

end

"""
    Create_Jacobian_NR_CPF(CDF_DF_List_pu, Ybus, SolutionVector_V,
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
function Create_Jacobian_NR_CPF(CDF_DF_List_pu, Ybus, SolutionVector_V, SolutionVector_Delta, PQ_BusArray, NR_Type, ContinuationPowerFlow_Indicator)

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

                    J_12[ii,jj] = ((PQ_BusArray[ii+1,1]) + (SolutionVector_V[ii+1]^(2)*real(Ybus[ii+1,ii+1])))/SolutionVector_V[ii+1]

                else # Off-Diagonal Term

                    J_12[ii,jj] = (SolutionVector_V[ii+1]*SolutionVector_V[jj+1]*abs(Ybus[ii+1,jj+1]) *cos(angle(Ybus[ii+1,jj+1]) + deg2rad(SolutionVector_Delta[jj+1]) - deg2rad(SolutionVector_Delta[ii+1])))/SolutionVector_V[jj+1]

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

                    J_22[ii,jj] = ((PQ_BusArray[ii+1,2]) - (SolutionVector_V[ii+1]^(2)*imag(Ybus[ii+1,ii+1])))/SolutionVector_V[ii+1]

                else # Off-Diagonal Term

                    J_22[ii,jj] = -(SolutionVector_V[ii+1]*SolutionVector_V[jj+1]*abs(Ybus[ii+1,jj+1]) *sin(angle(Ybus[ii+1,jj+1]) + deg2rad(SolutionVector_Delta[jj+1]) - deg2rad(SolutionVector_Delta[ii+1])))/SolutionVector_V[jj+1]

                end

            end

        end

        # Creating Jacobian_NR
        Jacobian_NR1 = hcat(J_11,J_12)
        Jacobian_NR2 = hcat(J_21,J_22)

        Jacobian_NR = vcat(Jacobian_NR1,Jacobian_NR2)

        # Addressing Machine Precision Problem
        for ii in 1:size(Jacobian_NR)[1]

            for jj in 1:size(Jacobian_NR)[2]

                if (abs(Jacobian_NR[ii,jj]) < 1e-12)

                    Jacobian_NR[ii,jj] = 0

                end

            end

        end

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

            # Addressing Machine Precision Problem
            for ii in 1:size(J_11)[1]

                for jj in 1:size(J_11)[2]

                    if (abs(J_11[ii,jj]) < 1e-12)

                        J_11[ii,jj] = 0

                    end

                end

            end

            for ii in 1:size(J_22)[1]

                for jj in 1:size(J_22)[2]

                    if (abs(J_22[ii,jj]) < 1e-12)

                        J_22[ii,jj] = 0

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

        # Addressing Machine Precision Problem
        for ii in 1:size(J_11)[1]

            for jj in 1:size(J_11)[2]

                if (abs(J_11[ii,jj]) < 1e-12)

                    J_11[ii,jj] = 0

                end

            end

        end

        for ii in 1:size(J_22)[1]

            for jj in 1:size(J_22)[2]

                if (abs(J_22[ii,jj]) < 1e-12)

                    J_22[ii,jj] = 0

                end

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
    Jacobian_NR_Predict = Create_Jacobian_NR_CPF(CDF_DF_List_pu, Ybus, SolutionVector_V, SolutionVector_Delta, PQ_BusArray, NR_Type, 1)

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
    Jacobian_NR_Correct = Create_Jacobian_NR_CPF(CDF_DF_List_pu, Ybus, SolutionVector_V, SolutionVector_Delta, PQ_BusArray, NR_Type, 1)

    # Augmenting NR Jacobian with K Vector (Debugger -K_Vector)
    Jacobian_CPF_Correct = hcat(Jacobian_NR_Correct, -K_Vector)

    # Creating Last row of CPF Jacobian
    Jacobian_CPF_Correct_LastRow = zeros(1,size(Jacobian_CPF_Correct)[2])

    Jacobian_CPF_Correct_LastRow[1,Index_CPF] = 1

    # Augmenting Last Row to CPF Jacobian
    Jacobian_CPF_Correct = vcat(Jacobian_CPF_Correct, Jacobian_CPF_Correct_LastRow)

    return Jacobian_CPF_Correct

end

""" 
Compute_H_Matrix_SE(CDF_DF_List_pu, SolutionVector_V, SolutionVector_Delta, Ybus, IncidenceMatrix_A, Detected_BadData_Vector)

Creates Initial Solution Vector for the Power System State Estimation.

'''     
# Arguments 

'''
'''
# Output

'''
"""
function Compute_H_Matrix_SE(CDF_DF_List_pu, SolutionVector_V, SolutionVector_Delta, Ybus, IncidenceMatrix_A, Detected_BadData_Vector)

    # Getting required data from CDF_DF_List
    BusDataCard_DF = CDF_DF_List_pu[2]
    BranchDataCard_DF = CDF_DF_List_pu[3]

    # Number of Buses and Lines
    N_Bus = nrow(BusDataCard_DF)
    N_Lines = nrow(BranchDataCard_DF)    
    
    # Getting the contents of Detected_BadData_Vector
    V_BadData_Vector = Detected_BadData_Vector[1]
    P_BadData_Vector = Detected_BadData_Vector[2]
    Q_BadData_Vector = Detected_BadData_Vector[3]
    Pij_BadData_Vector = Detected_BadData_Vector[4]
    Pji_BadData_Vector = Detected_BadData_Vector[5]
    Qij_BadData_Vector = Detected_BadData_Vector[6]
    Qji_BadData_Vector = Detected_BadData_Vector[7]

    # Computing number of Columns for the H Matrix Submatrices
    H_Matrix_Odd_ColNum = N_Bus-1
    H_Matrix_Even_ColNum = N_Bus

    # Computing number of Rows for the H Matrix Submatrices
    H_Matrix_Sub_RowNum_Vector = zeros(length(Detected_BadData_Vector),1)  # Initialization
    
    for ii in 1:length(Detected_BadData_Vector)  # For each element array in Detected_BadData_Vector

        Current_Matrix_RowNum = 0  # Initialization

        for jj in 1:length(Detected_BadData_Vector[ii])  # For each element in Detected_BadData_Vector[ii]

            # If Loop: To check for bad data
            if (Detected_BadData_Vector[ii][jj] == 0)  # Bad Data present

                # Do not Increment Size

            else  # Bad Data not present

                # Increment Current_Matrix_RowNum
                Current_Matrix_RowNum = Current_Matrix_RowNum + 1

            end

        end

        # Updating H_Matrix_Sub_RowNum_Vector
        H_Matrix_Sub_RowNum_Vector[ii,1] = Current_Matrix_RowNum

    end

    # Initializing H Matrix Submatrices
    H_1_SubMatrix = zeros(Int64(H_Matrix_Sub_RowNum_Vector[1,1]), N_Bus-1)  # Already Built
    H_2_SubMatrix = zeros(Int64(H_Matrix_Sub_RowNum_Vector[1,1]), N_Bus)
    H_3_SubMatrix = zeros(Int64(H_Matrix_Sub_RowNum_Vector[2,1]), N_Bus-1)
    H_4_SubMatrix = zeros(Int64(H_Matrix_Sub_RowNum_Vector[2,1]), N_Bus)
    H_5_SubMatrix = zeros(Int64(H_Matrix_Sub_RowNum_Vector[3,1]), N_Bus-1,)
    H_6_SubMatrix = zeros(Int64(H_Matrix_Sub_RowNum_Vector[3,1]), N_Bus)
    H_7_SubMatrix = zeros(Int64(H_Matrix_Sub_RowNum_Vector[4,1]), N_Bus-1)
    H_8_SubMatrix = zeros(Int64(H_Matrix_Sub_RowNum_Vector[4,1]), N_Bus)
    H_9_SubMatrix = zeros(Int64(H_Matrix_Sub_RowNum_Vector[5,1]), N_Bus-1)
    H_10_SubMatrix = zeros(Int64(H_Matrix_Sub_RowNum_Vector[5,1]), N_Bus)
    H_11_SubMatrix = zeros(Int64(H_Matrix_Sub_RowNum_Vector[6,1]), N_Bus-1)
    H_12_SubMatrix = zeros(Int64(H_Matrix_Sub_RowNum_Vector[6,1]), N_Bus)
    H_13_SubMatrix = zeros(Int64(H_Matrix_Sub_RowNum_Vector[7,1]), N_Bus-1)
    H_14_SubMatrix = zeros(Int64(H_Matrix_Sub_RowNum_Vector[7,1]), N_Bus)

    # Building H_2_SubMatrix
    H_2_SubMatrix_ii = 0  # Initialization
    for jj in 1:length(V_BadData_Vector)  # For each element in V_BadData_Vector[ii]

        # If Loop: To check for bad data
        if (V_BadData_Vector[jj] == 0)  # Bad Data present

            # Do not Increment Size_R_Inv_Matrix

        else  # Bad Data not present

            # Incrementing H_2_SubMatrix_ii
            H_2_SubMatrix_ii = H_2_SubMatrix_ii + 1

            # Compute H_2_SubMatrix Element
            H_2_SubMatrix[H_2_SubMatrix_ii, jj] = 1

        end

    end

    # Building H_3_SubMatrix
    H_3_SubMatrix_ii = 0  # Initialization
    for ii in 1:length(P_BadData_Vector)  # For each element in P_BadData_Vector[ii]

        # If Loop: To check for bad data
        if (P_BadData_Vector[ii] == 0)  # Bad Data present

            # Do not Increment Size_R_Inv_Matrix

        else  # Bad Data not present

            # Incrementing H_3_SubMatrix_ii
            H_3_SubMatrix_ii = H_3_SubMatrix_ii + 1

            # Compute H_3_SubMatrix Elements in the H_3_SubMatrix_ii Row
            for jj in 1:H_Matrix_Odd_ColNum  # Through each column of H3

                # Getting New_jj to correct for absence of Delta_1
                New_jj = jj + 1

                # To check for Diagonal and Off-Diagonal Terms
                if (ii == New_jj)  # Diagonal Terms

                    # Computing element H_3_SubMatrix[H_3_SubMatrix_ii, jj]
                    for kk in 1:N_Bus  # For each bus

                        if (ii == kk)

                            continue

                        else

                            H_3_SubMatrix[H_3_SubMatrix_ii, jj] = H_3_SubMatrix[H_3_SubMatrix_ii, jj] + (SolutionVector_V[ii] * SolutionVector_V[kk] * abs(Ybus[ii,kk]) * sin(angle(Ybus[ii,kk]) + deg2rad(SolutionVector_Delta[kk]) - deg2rad(SolutionVector_Delta[ii])))

                        end
                        
                    end

                else  # Off-Diagonal Terms

                    # Computing element H_3_SubMatrix[H_3_SubMatrix_ii, jj]
                    H_3_SubMatrix[H_3_SubMatrix_ii, jj] = -(SolutionVector_V[ii] * SolutionVector_V[New_jj] * abs(Ybus[ii,New_jj]) * sin(angle(Ybus[ii,New_jj]) + deg2rad(SolutionVector_Delta[New_jj]) - deg2rad(SolutionVector_Delta[ii])))

                end

            end

        end

    end

    # Building H_4_SubMatrix
    H_4_SubMatrix_ii = 0  # Initialization
    for ii in 1:length(P_BadData_Vector)  # For each element in P_BadData_Vector[ii]

        # If Loop: To check for bad data
        if (P_BadData_Vector[ii] == 0)  # Bad Data present

            # Do not Increment Size_R_Inv_Matrix

        else  # Bad Data not present

            # Incrementing H_4_SubMatrix_ii
            H_4_SubMatrix_ii = H_4_SubMatrix_ii + 1

            # Compute H_3_SubMatrix Elements in the H_4_SubMatrix_ii Row
            for jj in 1:H_Matrix_Even_ColNum  # Through each column of H4  
                
                # Getting New_jj to correct for absence of Delta_1
                New_jj = jj 

                # To check for Diagonal and Off-Diagonal Terms
                if (ii == New_jj)  # Diagonal Terms

                    # Computing element H_4_SubMatrix[H_4_SubMatrix_ii, jj]
                    H_4_SubMatrix[H_4_SubMatrix_ii, jj] = H_4_SubMatrix[H_4_SubMatrix_ii, jj] + (2 * SolutionVector_V[ii] * real(Ybus[ii,ii]))
                    for kk in 1:N_Bus  # For each bus

                        if (ii == kk)

                            continue

                        else

                            H_4_SubMatrix[H_4_SubMatrix_ii, jj] = H_4_SubMatrix[H_4_SubMatrix_ii, jj] + (SolutionVector_V[kk] * abs(Ybus[ii,kk]) * cos(angle(Ybus[ii,kk]) + deg2rad(SolutionVector_Delta[kk]) - deg2rad(SolutionVector_Delta[ii])))

                        end
                        
                    end

                else  # Off-Diagonal Terms

                    # Computing element H_4_SubMatrix[H_4_SubMatrix_ii, jj]
                    H_4_SubMatrix[H_4_SubMatrix_ii, jj] = (SolutionVector_V[ii] * abs(Ybus[ii,New_jj]) * cos(angle(Ybus[ii,New_jj]) + deg2rad(SolutionVector_Delta[New_jj]) - deg2rad(SolutionVector_Delta[ii])))

                end

            end
            
        end

    end

    # Building H_5_SubMatrix
    H_5_SubMatrix_ii = 0  # Initialization
    for ii in 1:length(Q_BadData_Vector)  # For each element in Q_BadData_Vector[ii]

        # If Loop: To check for bad data
        if (Q_BadData_Vector[ii] == 0)  # Bad Data present

            # Do not Increment Size_R_Inv_Matrix

        else  # Bad Data not present

            # Incrementing H_5_SubMatrix_ii
            H_5_SubMatrix_ii = H_5_SubMatrix_ii + 1

            # Compute H_5_SubMatrix Elements in the H_5_SubMatrix_ii Row
            for jj in 1:H_Matrix_Odd_ColNum  # Through each column of H5

                # Getting New_jj to correct for absence of Delta_1
                New_jj = jj + 1

                # To check for Diagonal and Off-Diagonal Terms
                if (ii == New_jj)  # Diagonal Terms

                    # Computing element H_5_SubMatrix[H_5_SubMatrix_ii, jj]
                    for kk in 1:N_Bus  # For each bus

                        if (ii == kk)

                            continue

                        else

                            H_5_SubMatrix[H_5_SubMatrix_ii, jj] = H_5_SubMatrix[H_5_SubMatrix_ii, jj] + (SolutionVector_V[ii] * SolutionVector_V[kk] * abs(Ybus[ii,kk]) * cos(angle(Ybus[ii,kk]) + deg2rad(SolutionVector_Delta[kk]) - deg2rad(SolutionVector_Delta[ii])))

                        end
                        
                    end

                else  # Off-Diagonal Terms

                    # Computing element H_5_SubMatrix[H_5_SubMatrix_ii, jj]
                    H_5_SubMatrix[H_5_SubMatrix_ii, jj] = -(SolutionVector_V[ii] * SolutionVector_V[New_jj] * abs(Ybus[ii,New_jj]) * cos(angle(Ybus[ii,New_jj]) + deg2rad(SolutionVector_Delta[New_jj]) - deg2rad(SolutionVector_Delta[ii])))

                end

            end

        end

    end

    # Building H_6_SubMatrix
    H_6_SubMatrix_ii = 0  # Initialization
    for ii in 1:length(Q_BadData_Vector)  # For each element in Q_BadData_Vector[ii]

        # If Loop: To check for bad data
        if (Q_BadData_Vector[ii] == 0)  # Bad Data present

            # Do not Increment Size_R_Inv_Matrix

        else  # Bad Data not present

            # Incrementing H_6_SubMatrix_ii
            H_6_SubMatrix_ii = H_6_SubMatrix_ii + 1

            # Compute H_6_SubMatrix Elements in the H_6_SubMatrix_ii Row
            for jj in 1:H_Matrix_Even_ColNum  # Through each column of H4  
                
                # Getting New_jj to correct for absence of Delta_1
                New_jj = jj 

                # To check for Diagonal and Off-Diagonal Terms
                if (ii == New_jj)  # Diagonal Terms

                    # Computing element H_6_SubMatrix[H_6_SubMatrix_ii, jj]
                    H_6_SubMatrix[H_6_SubMatrix_ii, jj] = H_6_SubMatrix[H_6_SubMatrix_ii, jj] - (2 * SolutionVector_V[ii] * imag(Ybus[ii,ii]))
                    for kk in 1:N_Bus  # For each bus

                        if (ii == kk)

                            continue

                        else

                            H_6_SubMatrix[H_6_SubMatrix_ii, jj] = H_6_SubMatrix[H_6_SubMatrix_ii, jj] - (SolutionVector_V[kk] * abs(Ybus[ii,kk]) * sin(angle(Ybus[ii,kk]) + deg2rad(SolutionVector_Delta[kk]) - deg2rad(SolutionVector_Delta[ii])))

                        end
                        
                    end

                else  # Off-Diagonal Terms

                    # Computing element H_6_SubMatrix[H_6_SubMatrix_ii, jj]
                    H_6_SubMatrix[H_6_SubMatrix_ii, jj] = -(SolutionVector_V[ii] * abs(Ybus[ii,New_jj]) * sin(angle(Ybus[ii,New_jj]) + deg2rad(SolutionVector_Delta[New_jj]) - deg2rad(SolutionVector_Delta[ii])))

                end

            end
            
        end

    end

    # Building H_7_SubMatrix
    H_7_SubMatrix_ii = 0  # Initialization
    for ii in 1:length(Pij_BadData_Vector)  # For each element in Pij_BadData_Vector[ii]

        # If Loop: To check for bad data
        if (Pij_BadData_Vector[ii] == 0)  # Bad Data present

            # Do not Increment Size_R_Inv_Matrix

        else  # Bad Data not present

            # Incrementing H_7_SubMatrix_ii
            H_7_SubMatrix_ii = H_7_SubMatrix_ii + 1

            ## Compute H_7_SubMatrix Elements in the H_7_SubMatrix_ii Row     
            
            # Getting Current IncidenceMatrix_A Row
            Current_A_Row = IncidenceMatrix_A[ii,:]

            # Computing Current measurements i and j Node Indices
            Bus_i_Index = 0
            Bus_j_Index = 0
            for kk in 1:length(Current_A_Row)  # For Bus_i_Index
                    if (Current_A_Row[kk] == 1)
                            Bus_i_Index = kk
                            break
                    end
            end
            for kk in 1:length(Current_A_Row)  # For Bus_j_Index
                    if (Current_A_Row[kk] == -1)
                            Bus_j_Index = kk
                            break
                    end
            end

            # Computing Derivative w.r.t Delta_i
            if (Bus_i_Index != 1)

                # Computing New Bus i Index (for positioning)
                New_Bus_i_Index = Bus_i_Index -1

                # Computing element H_7_SubMatrix[H_7_SubMatrix_ii, New_Bus_i_Index] 
                H_7_SubMatrix[H_7_SubMatrix_ii, New_Bus_i_Index] = (SolutionVector_V[Bus_i_Index] * SolutionVector_V[Bus_j_Index] * abs(Ybus[Bus_i_Index,Bus_j_Index]) * sin(angle(Ybus[Bus_i_Index,Bus_j_Index]) + deg2rad(SolutionVector_Delta[Bus_j_Index]) - deg2rad(SolutionVector_Delta[Bus_i_Index])))

            end

            # Computing Derivative w.r.t Delta_j
            if (Bus_j_Index != 1)

                # Computing New Bus j Index (for positioning)
                New_Bus_j_Index = Bus_j_Index -1

                # Computing element H_7_SubMatrix[H_7_SubMatrix_ii, New_Bus_j_Index] 
                H_7_SubMatrix[H_7_SubMatrix_ii, New_Bus_j_Index] = -(SolutionVector_V[Bus_i_Index] * SolutionVector_V[Bus_j_Index] * abs(Ybus[Bus_i_Index,Bus_j_Index]) * sin(angle(Ybus[Bus_i_Index,Bus_j_Index]) + deg2rad(SolutionVector_Delta[Bus_j_Index]) - deg2rad(SolutionVector_Delta[Bus_i_Index])))

            end          

        end

    end

    # Building H_8_SubMatrix
    H_8_SubMatrix_ii = 0  # Initialization
    for ii in 1:length(Pij_BadData_Vector)  # For each element in Pij_BadData_Vector[ii]

        # If Loop: To check for bad data
        if (Pij_BadData_Vector[ii] == 0)  # Bad Data present

            # Do not Increment Size_R_Inv_Matrix

        else  # Bad Data not present

            # Incrementing H_8_SubMatrix_ii
            H_8_SubMatrix_ii = H_8_SubMatrix_ii + 1

            ## Compute H_8_SubMatrix Elements in the H_8_SubMatrix_ii Row     
            
            # Getting Current IncidenceMatrix_A Row
            Current_A_Row = IncidenceMatrix_A[ii,:]

            # Computing Current measurements i and j Node Indices
            Bus_i_Index = 0
            Bus_j_Index = 0
            for kk in 1:length(Current_A_Row)  # For Bus_i_Index
                    if (Current_A_Row[kk] == 1)
                            Bus_i_Index = kk
                            break
                    end
            end
            for kk in 1:length(Current_A_Row)  # For Bus_j_Index
                    if (Current_A_Row[kk] == -1)
                            Bus_j_Index = kk
                            break
                    end
            end

            # Computing Derivative w.r.t V_i

            # Computing New Bus i Index (for positioning)
            New_Bus_i_Index = Bus_i_Index 

            # Computing element H_8_SubMatrix[H_8_SubMatrix_ii, New_Bus_i_Index] 
            H_8_SubMatrix[H_8_SubMatrix_ii, New_Bus_i_Index] = (2 * SolutionVector_V[Bus_i_Index] * real(Ybus[Bus_i_Index,Bus_j_Index])) + (SolutionVector_V[Bus_j_Index] * abs(Ybus[Bus_i_Index,Bus_j_Index]) * cos(angle(Ybus[Bus_i_Index,Bus_j_Index]) + deg2rad(SolutionVector_Delta[Bus_j_Index]) - deg2rad(SolutionVector_Delta[Bus_i_Index])))

            # Computing Derivative w.r.t V_j

            # Computing New Bus j Index (for positioning)
            New_Bus_j_Index = Bus_j_Index 

            # Computing element H_8_SubMatrix[H_8_SubMatrix_ii, New_Bus_j_Index] 
            H_8_SubMatrix[H_8_SubMatrix_ii, New_Bus_j_Index] = (SolutionVector_V[Bus_i_Index] * abs(Ybus[Bus_i_Index,Bus_j_Index]) * cos(angle(Ybus[Bus_i_Index,Bus_j_Index]) + deg2rad(SolutionVector_Delta[Bus_j_Index]) - deg2rad(SolutionVector_Delta[Bus_i_Index])))
        
        end

    end

    # Building H_9_SubMatrix
    H_9_SubMatrix_ii = 0  # Initialization
    for ii in 1:length(Pji_BadData_Vector)  # For each element in Pji_BadData_Vector[ii]

        # If Loop: To check for bad data
        if (Pji_BadData_Vector[ii] == 0)  # Bad Data present

            # Do not Increment Size_R_Inv_Matrix

        else  # Bad Data not present

            # Incrementing H_9_SubMatrix_ii
            H_9_SubMatrix_ii = H_9_SubMatrix_ii + 1

            ## Compute H_9_SubMatrix Elements in the H_9_SubMatrix_ii Row     
            
            # Getting Current IncidenceMatrix_A Row
            Current_A_Row = IncidenceMatrix_A[ii,:]

            # Computing Current measurements i and j Node Indices
            Bus_i_Index = 0
            Bus_j_Index = 0
            for kk in 1:length(Current_A_Row)  # For Bus_i_Index
                    if (Current_A_Row[kk] == 1)
                            Bus_i_Index = kk
                            break
                    end
            end
            for kk in 1:length(Current_A_Row)  # For Bus_j_Index
                    if (Current_A_Row[kk] == -1)
                            Bus_j_Index = kk
                            break
                    end
            end

            # Computing Derivative w.r.t Delta_i
            if (Bus_i_Index != 1)

                # Computing New Bus i Index (for positioning)
                New_Bus_i_Index = Bus_i_Index -1

                # Computing element H_9_SubMatrix[H_9_SubMatrix_ii, New_Bus_i_Index] 
                H_9_SubMatrix[H_9_SubMatrix_ii, New_Bus_i_Index] = -(SolutionVector_V[Bus_i_Index] * SolutionVector_V[Bus_j_Index] * abs(Ybus[Bus_i_Index,Bus_j_Index]) * sin(angle(Ybus[Bus_i_Index,Bus_j_Index]) + deg2rad(SolutionVector_Delta[Bus_i_Index]) - deg2rad(SolutionVector_Delta[Bus_j_Index])))

            end

            # Computing Derivative w.r.t Delta_j
            if (Bus_j_Index != 1)

                # Computing New Bus j Index (for positioning)
                New_Bus_j_Index = Bus_j_Index -1

                # Computing element H_9_SubMatrix[H_9_SubMatrix_ii, New_Bus_j_Index] 
                H_9_SubMatrix[H_9_SubMatrix_ii, New_Bus_j_Index] = (SolutionVector_V[Bus_i_Index] * SolutionVector_V[Bus_j_Index] * abs(Ybus[Bus_i_Index,Bus_j_Index]) * sin(angle(Ybus[Bus_i_Index,Bus_j_Index]) + deg2rad(SolutionVector_Delta[Bus_i_Index]) - deg2rad(SolutionVector_Delta[Bus_j_Index])))

            end          

        end

    end

    # Building H_10_SubMatrix
    H_10_SubMatrix_ii = 0  # Initialization
    for ii in 1:length(Pji_BadData_Vector)  # For each element in Pji_BadData_Vector[ii]

        # If Loop: To check for bad data
        if (Pji_BadData_Vector[ii] == 0)  # Bad Data present

            # Do not Increment Size_R_Inv_Matrix

        else  # Bad Data not present

            # Incrementing H_10_SubMatrix_ii
            H_10_SubMatrix_ii = H_10_SubMatrix_ii + 1

            ## Compute H_10_SubMatrix Elements in the H_10_SubMatrix_ii Row     
            
            # Getting Current IncidenceMatrix_A Row
            Current_A_Row = IncidenceMatrix_A[ii,:]

            # Computing Current measurements i and j Node Indices
            Bus_i_Index = 0
            Bus_j_Index = 0
            for kk in 1:length(Current_A_Row)  # For Bus_i_Index
                    if (Current_A_Row[kk] == 1)
                            Bus_i_Index = kk
                            break
                    end
            end
            for kk in 1:length(Current_A_Row)  # For Bus_j_Index
                    if (Current_A_Row[kk] == -1)
                            Bus_j_Index = kk
                            break
                    end
            end

            # Computing Derivative w.r.t V_i

            # Computing New Bus i Index (for positioning)
            New_Bus_i_Index = Bus_i_Index 

            # Computing element H_10_SubMatrix[H_10_SubMatrix_ii, New_Bus_i_Index] 
            H_10_SubMatrix[H_10_SubMatrix_ii, New_Bus_i_Index] =  (SolutionVector_V[Bus_j_Index] * abs(Ybus[Bus_i_Index,Bus_j_Index]) * cos(angle(Ybus[Bus_i_Index,Bus_j_Index]) + deg2rad(SolutionVector_Delta[Bus_i_Index]) - deg2rad(SolutionVector_Delta[Bus_j_Index])))

            # Computing Derivative w.r.t V_j

            # Computing New Bus j Index (for positioning)
            New_Bus_j_Index = Bus_j_Index 

            # Computing element H_10_SubMatrix[H_10_SubMatrix_ii, New_Bus_j_Index] 
            H_10_SubMatrix[H_10_SubMatrix_ii, New_Bus_j_Index] = (2 * SolutionVector_V[Bus_j_Index] * real(Ybus[Bus_i_Index,Bus_j_Index])) + (SolutionVector_V[Bus_i_Index] * abs(Ybus[Bus_i_Index,Bus_j_Index]) * cos(angle(Ybus[Bus_i_Index,Bus_j_Index]) + deg2rad(SolutionVector_Delta[Bus_i_Index]) - deg2rad(SolutionVector_Delta[Bus_j_Index])))
        
        end

    end

    # Building H_11_SubMatrix
    H_11_SubMatrix_ii = 0  # Initialization
    for ii in 1:length(Qij_BadData_Vector)  # For each element in Qij_BadData_Vector[ii]

        # If Loop: To check for bad data
        if (Qij_BadData_Vector[ii] == 0)  # Bad Data present

            # Do not Increment Size_R_Inv_Matrix

        else  # Bad Data not present

            # Incrementing H_11_SubMatrix_ii
            H_11_SubMatrix_ii = H_11_SubMatrix_ii + 1

            ## Compute H_11_SubMatrix Elements in the H_11_SubMatrix_ii Row     
            
            # Getting Current IncidenceMatrix_A Row
            Current_A_Row = IncidenceMatrix_A[ii,:]

            # Computing Current measurements i and j Node Indices
            Bus_i_Index = 0
            Bus_j_Index = 0
            for kk in 1:length(Current_A_Row)  # For Bus_i_Index
                    if (Current_A_Row[kk] == 1)
                            Bus_i_Index = kk
                            break
                    end
            end
            for kk in 1:length(Current_A_Row)  # For Bus_j_Index
                    if (Current_A_Row[kk] == -1)
                            Bus_j_Index = kk
                            break
                    end
            end

            # Computing Derivative w.r.t Delta_i
            if (Bus_i_Index != 1)

                # Computing New Bus i Index (for positioning)
                New_Bus_i_Index = Bus_i_Index -1

                # Computing element H_11_SubMatrix[H_11_SubMatrix_ii, New_Bus_i_Index] 
                H_11_SubMatrix[H_11_SubMatrix_ii, New_Bus_i_Index] = (SolutionVector_V[Bus_i_Index] * SolutionVector_V[Bus_j_Index] * abs(Ybus[Bus_i_Index,Bus_j_Index]) * cos(angle(Ybus[Bus_i_Index,Bus_j_Index]) + deg2rad(SolutionVector_Delta[Bus_j_Index]) - deg2rad(SolutionVector_Delta[Bus_i_Index])))

            end

            # Computing Derivative w.r.t Delta_j
            if (Bus_j_Index != 1)

                # Computing New Bus j Index (for positioning)
                New_Bus_j_Index = Bus_j_Index -1

                # Computing element H_11_SubMatrix[H_11_SubMatrix_ii, New_Bus_j_Index] 
                H_11_SubMatrix[H_11_SubMatrix_ii, New_Bus_j_Index] = -(SolutionVector_V[Bus_i_Index] * SolutionVector_V[Bus_j_Index] * abs(Ybus[Bus_i_Index,Bus_j_Index]) * cos(angle(Ybus[Bus_i_Index,Bus_j_Index]) + deg2rad(SolutionVector_Delta[Bus_j_Index]) - deg2rad(SolutionVector_Delta[Bus_i_Index])))

            end          

        end

    end

    # Building H_12_SubMatrix
    H_12_SubMatrix_ii = 0  # Initialization
    for ii in 1:length(Qij_BadData_Vector)  # For each element in Qij_BadData_Vector[ii]

        # If Loop: To check for bad data
        if (Qij_BadData_Vector[ii] == 0)  # Bad Data present

            # Do not Increment Size_R_Inv_Matrix

        else  # Bad Data not present

            # Incrementing H_12_SubMatrix_ii
            H_12_SubMatrix_ii = H_12_SubMatrix_ii + 1

            ## Compute H_12_SubMatrix Elements in the H_12_SubMatrix_ii Row     
            
            # Getting Current IncidenceMatrix_A Row
            Current_A_Row = IncidenceMatrix_A[ii,:]

            # Computing Current measurements i and j Node Indices
            Bus_i_Index = 0
            Bus_j_Index = 0
            for kk in 1:length(Current_A_Row)  # For Bus_i_Index
                    if (Current_A_Row[kk] == 1)
                            Bus_i_Index = kk
                            break
                    end
            end
            for kk in 1:length(Current_A_Row)  # For Bus_j_Index
                    if (Current_A_Row[kk] == -1)
                            Bus_j_Index = kk
                            break
                    end
            end

            # Computing Derivative w.r.t V_i

            # Computing New Bus i Index (for positioning)
            New_Bus_i_Index = Bus_i_Index 

            # Computing element H_12_SubMatrix[H_12_SubMatrix_ii, New_Bus_i_Index] 
            H_12_SubMatrix[H_12_SubMatrix_ii, New_Bus_i_Index] =  -(2 * SolutionVector_V[Bus_i_Index] * ((BranchDataCard_DF.B_pu[ii]/2) + imag(Ybus[Bus_i_Index,Bus_j_Index]))) - (SolutionVector_V[Bus_j_Index] * abs(Ybus[Bus_i_Index,Bus_j_Index]) * sin(angle(Ybus[Bus_i_Index,Bus_j_Index]) + deg2rad(SolutionVector_Delta[Bus_j_Index]) - deg2rad(SolutionVector_Delta[Bus_i_Index])))

            # Computing Derivative w.r.t V_j

            # Computing New Bus j Index (for positioning)
            New_Bus_j_Index = Bus_j_Index 

            # Computing element H_12_SubMatrix[H_12_SubMatrix_ii, New_Bus_j_Index] 
            H_12_SubMatrix[H_12_SubMatrix_ii, New_Bus_j_Index] =  -(SolutionVector_V[Bus_i_Index] * abs(Ybus[Bus_i_Index,Bus_j_Index]) * sin(angle(Ybus[Bus_i_Index,Bus_j_Index]) + deg2rad(SolutionVector_Delta[Bus_j_Index]) - deg2rad(SolutionVector_Delta[Bus_i_Index])))
        
        end

    end

    # Building H_13_SubMatrix
    H_13_SubMatrix_ii = 0  # Initialization
    for ii in 1:length(Qji_BadData_Vector)  # For each element in Qji_BadData_Vector[ii]

        # If Loop: To check for bad data
        if (Qji_BadData_Vector[ii] == 0)  # Bad Data present

            # Do not Increment Size_R_Inv_Matrix

        else  # Bad Data not present

            # Incrementing H_13_SubMatrix_ii
            H_13_SubMatrix_ii = H_13_SubMatrix_ii + 1

            ## Compute H_13_SubMatrix Elements in the H_13_SubMatrix_ii Row     
            
            # Getting Current IncidenceMatrix_A Row
            Current_A_Row = IncidenceMatrix_A[ii,:]

            # Computing Current measurements i and j Node Indices
            Bus_i_Index = 0
            Bus_j_Index = 0
            for kk in 1:length(Current_A_Row)  # For Bus_i_Index
                    if (Current_A_Row[kk] == 1)
                            Bus_i_Index = kk
                            break
                    end
            end
            for kk in 1:length(Current_A_Row)  # For Bus_j_Index
                    if (Current_A_Row[kk] == -1)
                            Bus_j_Index = kk
                            break
                    end
            end

            # Computing Derivative w.r.t Delta_i
            if (Bus_i_Index != 1)

                # Computing New Bus i Index (for positioning)
                New_Bus_i_Index = Bus_i_Index -1

                # Computing element H_13_SubMatrix[H_13_SubMatrix_ii, New_Bus_i_Index] 
                H_13_SubMatrix[H_13_SubMatrix_ii, New_Bus_i_Index] = -(SolutionVector_V[Bus_i_Index] * SolutionVector_V[Bus_j_Index] * abs(Ybus[Bus_i_Index,Bus_j_Index]) * cos(angle(Ybus[Bus_i_Index,Bus_j_Index]) + deg2rad(SolutionVector_Delta[Bus_i_Index]) - deg2rad(SolutionVector_Delta[Bus_j_Index])))

            end

            # Computing Derivative w.r.t Delta_j
            if (Bus_j_Index != 1)

                # Computing New Bus j Index (for positioning)
                New_Bus_j_Index = Bus_j_Index -1

                # Computing element H_13_SubMatrix[H_13_SubMatrix_ii, New_Bus_j_Index] 
                H_13_SubMatrix[H_13_SubMatrix_ii, New_Bus_j_Index] = (SolutionVector_V[Bus_i_Index] * SolutionVector_V[Bus_j_Index] * abs(Ybus[Bus_i_Index,Bus_j_Index]) * cos(angle(Ybus[Bus_i_Index,Bus_j_Index]) + deg2rad(SolutionVector_Delta[Bus_i_Index]) - deg2rad(SolutionVector_Delta[Bus_j_Index])))

            end          

        end

    end

    # Building H_14_SubMatrix
    H_14_SubMatrix_ii = 0  # Initialization
    for ii in 1:length(Qji_BadData_Vector)  # For each element in Qji_BadData_Vector[ii]

        # If Loop: To check for bad data
        if (Qji_BadData_Vector[ii] == 0)  # Bad Data present

            # Do not Increment Size_R_Inv_Matrix

        else  # Bad Data not present

            # Incrementing H_14_SubMatrix_ii
            H_14_SubMatrix_ii = H_14_SubMatrix_ii + 1

            ## Compute H_14_SubMatrix Elements in the H_14_SubMatrix_ii Row     
            
            # Getting Current IncidenceMatrix_A Row
            Current_A_Row = IncidenceMatrix_A[ii,:]

            # Computing Current measurements i and j Node Indices
            Bus_i_Index = 0
            Bus_j_Index = 0
            for kk in 1:length(Current_A_Row)  # For Bus_i_Index
                    if (Current_A_Row[kk] == 1)
                            Bus_i_Index = kk
                            break
                    end
            end
            for kk in 1:length(Current_A_Row)  # For Bus_j_Index
                    if (Current_A_Row[kk] == -1)
                            Bus_j_Index = kk
                            break
                    end
            end

            # Computing Derivative w.r.t V_i

            # Computing New Bus i Index (for positioning)
            New_Bus_i_Index = Bus_i_Index 

            # Computing element H_14_SubMatrix[H_14_SubMatrix_ii, New_Bus_i_Index] 
            H_14_SubMatrix[H_14_SubMatrix_ii, New_Bus_i_Index] =   - (SolutionVector_V[Bus_j_Index] * abs(Ybus[Bus_i_Index,Bus_j_Index]) * sin(angle(Ybus[Bus_i_Index,Bus_j_Index]) + deg2rad(SolutionVector_Delta[Bus_i_Index]) - deg2rad(SolutionVector_Delta[Bus_j_Index])))

            # Computing Derivative w.r.t V_j

            # Computing New Bus j Index (for positioning)
            New_Bus_j_Index = Bus_j_Index 

            # Computing element H_14_SubMatrix[H_14_SubMatrix_ii, New_Bus_j_Index] 
            H_14_SubMatrix[H_14_SubMatrix_ii, New_Bus_j_Index] =  -(2 * SolutionVector_V[Bus_j_Index] * ((BranchDataCard_DF.B_pu[ii]/2) + imag(Ybus[Bus_i_Index,Bus_j_Index]))) - (SolutionVector_V[Bus_i_Index] * abs(Ybus[Bus_i_Index,Bus_j_Index]) * sin(angle(Ybus[Bus_i_Index,Bus_j_Index]) + deg2rad(SolutionVector_Delta[Bus_i_Index]) - deg2rad(SolutionVector_Delta[Bus_j_Index])))
        
        end

    end

    # Developing Complete H_Matrix
    H_12_Matrix = hcat(H_1_SubMatrix, H_2_SubMatrix)
    H_34_Matrix = hcat(H_3_SubMatrix, H_4_SubMatrix)
    H_56_Matrix = hcat(H_5_SubMatrix, H_6_SubMatrix)
    H_78_Matrix = hcat(H_7_SubMatrix, H_8_SubMatrix)
    H_910_Matrix = hcat(H_9_SubMatrix, H_10_SubMatrix)
    H_1112_Matrix = hcat(H_11_SubMatrix, H_12_SubMatrix)
    H_1314_Matrix = hcat(H_13_SubMatrix, H_14_SubMatrix)

    H_Matrix = vcat(H_12_Matrix, H_34_Matrix, H_56_Matrix, H_78_Matrix, H_910_Matrix, H_1112_Matrix, H_1314_Matrix)

    # Addressing Machine Precision Problem
    for ii in 1:size(H_Matrix)[1]  # Through the Rows

        for jj in 1:size(H_Matrix)[2]  # Through the Columns

            if (abs(H_Matrix[ii, jj]) < 1e-12)

                H_Matrix[ii, jj] = 0

            end

        end

    end
    
    return H_Matrix 

end

"""
Compute_Del_g_Del_u_Matrix_OPF(CDF_DF_List_pu, u_P_Index_Vector)

Computes the Del_g/Del_u Matrix for optimal power flow of a power network.

'''
# Arguments
- '': 
'''
# Output
- '': 
'''
"""
function Compute_Del_g_Del_u_Matrix_OPF(CDF_DF_List_pu, u_P_Index_Vector)

    # Getting required data from CDF_DF_List
    BusDataCard_DF = CDF_DF_List_pu[2]

    # Number of Buses
    N_Bus = nrow(BusDataCard_DF)
    N_PQ_Bus = nrow(filter(row -> ((row.Type == 0) || (row.Type == 1)), BusDataCard_DF))
    N_PV_Bus = nrow(filter(row -> (row.Type == 2), BusDataCard_DF))
    N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

    N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))

    # Getting number of Control variables
    N_P = size(u_P_Index_Vector)[1]

    # Computing Number of load flow equations
    N_g_functions = 2*N_PQ_Bus + N_PV_Bus

    # Initializing Del_g_Del_u_Matrix
    Del_g_Del_u_Matrix = zeros(N_g_functions, N_P)

    # Updating Del_g_Del_u_Matrix
    for ii in 1:N_P  # For each element in u_P_Index_Vector

        # Computing Row Index for Del_g_Del_u_Matrix
        Del_g_Del_u_Matrix_RowIndex = Int64(u_P_Index_Vector[ii,1]) - 1

        # Updating Del_g_Del_u_Matrix element
        Del_g_Del_u_Matrix[Del_g_Del_u_Matrix_RowIndex,ii] = 1

    end

    # Addressing Machine Precision
    Del_g_Del_u_Matrix = Addressing_MachinePrecision(Del_g_Del_u_Matrix, 1e-12)

    return Del_g_Del_u_Matrix

end



"""
    Compute_Del_g_Del_x_Matrix_OPF(CDF_DF_List_pu, Ybus, SolutionVector_V,
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
function Compute_Del_g_Del_x_Matrix_OPF(CDF_DF_List_pu, Ybus, SolutionVector_V, SolutionVector_Delta, PQ_BusArray, NR_Type, ContinuationPowerFlow_Indicator)

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

                    J_12[ii,jj] = ((PQ_BusArray[ii+1,1]) + (SolutionVector_V[ii+1]^(2)*real(Ybus[ii+1,ii+1])))/SolutionVector_V[ii+1]

                else # Off-Diagonal Term

                    J_12[ii,jj] = (SolutionVector_V[ii+1]*SolutionVector_V[jj+1]*abs(Ybus[ii+1,jj+1]) *cos(angle(Ybus[ii+1,jj+1]) + deg2rad(SolutionVector_Delta[jj+1]) - deg2rad(SolutionVector_Delta[ii+1])))/SolutionVector_V[jj+1]

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

                    J_22[ii,jj] = ((PQ_BusArray[ii+1,2]) - (SolutionVector_V[ii+1]^(2)*imag(Ybus[ii+1,ii+1])))/SolutionVector_V[ii+1]

                else # Off-Diagonal Term

                    J_22[ii,jj] = -(SolutionVector_V[ii+1]*SolutionVector_V[jj+1]*abs(Ybus[ii+1,jj+1]) *sin(angle(Ybus[ii+1,jj+1]) + deg2rad(SolutionVector_Delta[jj+1]) - deg2rad(SolutionVector_Delta[ii+1])))/SolutionVector_V[jj+1]

                end

            end

        end

        # Creating Jacobian_NR
        Jacobian_NR1 = hcat(J_11,J_12)
        Jacobian_NR2 = hcat(J_21,J_22)

        Jacobian_NR11 = vcat(Jacobian_NR1,Jacobian_NR2)

        # Negating the Jacobian for Optimal Power Flow
        Jacobian_NR = -Jacobian_NR11

        # Addressing Machine Precision
        Jacobian_NR = Addressing_MachinePrecision(Jacobian_NR, 1e-12)

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

            # Addressing Machine Precision Problem
            for ii in 1:size(J_11)[1]

                for jj in 1:size(J_11)[2]

                    if (abs(J_11[ii,jj]) < 1e-12)

                        J_11[ii,jj] = 0

                    end

                end

            end

            for ii in 1:size(J_22)[1]

                for jj in 1:size(J_22)[2]

                    if (abs(J_22[ii,jj]) < 1e-12)

                        J_22[ii,jj] = 0

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

        # Addressing Machine Precision Problem
        for ii in 1:size(J_11)[1]

            for jj in 1:size(J_11)[2]

                if (abs(J_11[ii,jj]) < 1e-12)

                    J_11[ii,jj] = 0

                end

            end

        end

        for ii in 1:size(J_22)[1]

            for jj in 1:size(J_22)[2]

                if (abs(J_22[ii,jj]) < 1e-12)

                    J_22[ii,jj] = 0

                end

            end

        end

        # Creating Jacobian_NR
        Jacobian_NR = [J_11, J_22]

    end

    return Jacobian_NR

end

"""
    Compute_Del_f_Del_u_Matrix_OPF(CDF_DF_List_pu, u_P_Index_Vector)

Computes the Del_g/Del_u Matrix for optimal power flow of a power network.

'''
# Arguments
- '': 
'''
# Output
- '': 
'''
"""
function Compute_Del_f_Del_u_Matrix_OPF(SolutionVector_p_loop, Generator_CostCurve_Matrix_New)

    # Getting length of Del_f_Del_u_Matrix
    Len_Del_f_Del_u_Matrix = size(SolutionVector_p_loop)[1]

    # Getting length of Generator_CostCurve_Matrix_New
    Len_Generator_CostCurve_Matrix_New = size(Generator_CostCurve_Matrix_New)[1]

    # Initializing Del_f_Del_u_Matrix
    Del_f_Del_u_Matrix = zeros(Len_Del_f_Del_u_Matrix, 1)

    # Updating Del_f_Del_u_Matrix
    for ii in 1:Len_Del_f_Del_u_Matrix   # SolutionVector_p_loop

        # Checking if Slack-Bus is part of Cost Functions
        if (Len_Del_f_Del_u_Matrix != Len_Generator_CostCurve_Matrix_New)  # Slack-Bus is part of Cost Function

            # Updating Del_f_Del_u_Matrix element
            Del_f_Del_u_Matrix[ii,1] = (2 * (Generator_CostCurve_Matrix_New[ii+1,1]) * (SolutionVector_p_loop[ii,1])) + (Generator_CostCurve_Matrix_New[ii+1,2])

        else  # Slack-Bus is not part of Cost Function
            
            # Updating Del_f_Del_u_Matrix element
            Del_f_Del_u_Matrix[ii,1] = (2 * (Generator_CostCurve_Matrix_New[ii,1]) * (SolutionVector_p_loop[ii,1])) + (Generator_CostCurve_Matrix_New[ii,2])

        end
        
    end

    # Addressing Machine Precision
    Del_f_Del_u_Matrix = Addressing_MachinePrecision(Del_f_Del_u_Matrix, 1e-12)

    return Del_f_Del_u_Matrix

end

"""
    Compute_Del_f_Del_x_Matrix_OPF(CDF_DF_List_pu, SolutionVector_x_V_loop, SolutionVector_x_Delta_loop, SolutionVector_p_loop, SolutionVector_p_Full_loop, LineFlow_Array, Ybus, Generator_CostCurve_Matrix_New, Line_Index_Vector, Line_Bus_Index_Matrix, Line_P_Limit_Vector)

Sets the P Generation for Slack Bus in CDF DF Bus Data Card to values computed in converged Power Flow iteration.

'''
# Arguments
- '': 
'''
# Output
- '': 
'''
"""
function Compute_Del_f_Del_x_Matrix_OPF(CDF_DF_List_pu, SolutionVector_x_V_loop, SolutionVector_x_Delta_loop, SolutionVector_p_loop, SolutionVector_p_Full_loop, LineFlow_Array, Ybus, Generator_CostCurve_Matrix_New, Line_Index_Vector, Line_Bus_Index_Matrix, Line_P_Limit_Vector)

    # Getting required data from CDF_DF_List
    BusDataCard_DF = CDF_DF_List_pu[2]

    # Number of Buses
    N_Bus = nrow(BusDataCard_DF)
    N_PQ_Bus = nrow(filter(row -> ((row.Type == 0) || (row.Type == 1)), BusDataCard_DF))
    N_PV_Bus = nrow(filter(row -> (row.Type == 2), BusDataCard_DF))
    N_Slack_Bus = nrow(filter(row -> (row.Type == 3), BusDataCard_DF))
    
    # Computing Len_Delta and Len_V
    Len_Delta = N_PQ_Bus + N_PV_Bus
    Len_V = N_PQ_Bus
    
    # Initializing Del_f_Del_x_Matrix_Delta and Del_f_Del_x_Matrix_V
    Del_f_Del_x_Matrix_Delta = zeros(Len_Delta, 1)
    Del_f_Del_x_Matrix_V = zeros(Len_V, 1)

    # Getting Lengths of SolutionVector_p_loop and Generator_CostCurve_Matrix_New
    Len_SolutionVector_p_loop = size(SolutionVector_p_loop)[1]
    Len_Generator_CostCurve_Matrix_New = size(Generator_CostCurve_Matrix_New)[1]

    # Checking for Type of Cost Function  
    if ((Len_SolutionVector_p_loop == Len_Generator_CostCurve_Matrix_New) && (Line_Index_Vector == nothing) && (Line_Bus_Index_Matrix == nothing) && (Line_P_Limit_Vector == nothing))  # No Slack Bus in Cost Function and No Line Limits present

        # Constructing Del_f_Del_x_Matrix
        Del_f_Del_x_Matrix = vcat(Del_f_Del_x_Matrix_Delta, Del_f_Del_x_Matrix_V)

    elseif ((Len_SolutionVector_p_loop != Len_Generator_CostCurve_Matrix_New) && (Line_Index_Vector == nothing) && (Line_Bus_Index_Matrix == nothing) && (Line_P_Limit_Vector == nothing))  # Slack Bus in Cost Function and No Line Limits present

        ## Compute Slack Bus part of Del_f_Del_x_Matrix

        # Initializing Del_f_Del_x_Matrix_Delta_SlackBus and Del_f_Del_x_Matrix_V_SlackBus
        Del_f_Del_x_Matrix_Delta_SlackBus = zeros(Len_Delta, 1)
        Del_f_Del_x_Matrix_V_SlackBus = zeros(Len_V, 1)

        # Updating Del_f_Del_x_Matrix_Delta_SlackBus
        for ii in 1:Len_Delta  # For each element in Del_f_Del_x_Matrix_Delta_SlackBus

            Del_f_Del_x_Matrix_Delta_SlackBus[ii,1] = -((SolutionVector_x_V_loop[1,1] * SolutionVector_x_V_loop[ii+1,1] * abs(Ybus[1,ii+1])) *(sin(angle(Ybus[1,ii+1]) +  deg2rad(SolutionVector_x_Delta_loop[ii+1,1]) - deg2rad(SolutionVector_x_Delta_loop[1,1]))))

        end

        # Updating Del_f_Del_x_Matrix_V_SlackBus
        for ii in 1:Len_V  # For each element in Del_f_Del_x_Matrix_V_SlackBus

            Del_f_Del_x_Matrix_V_SlackBus[ii,1] = ((SolutionVector_x_V_loop[1,1] * abs(Ybus[1,ii+1])) *(cos(angle(Ybus[1,ii+1]) +  deg2rad(SolutionVector_x_Delta_loop[ii+1,1]) - deg2rad(SolutionVector_x_Delta_loop[1,1]))))

        end

        # Updating Del_f_Del_x_Matrix_Delta and Del_f_Del_x_Matrix_V_SlackBus
        Del_f_Del_x_Matrix_Delta = ((2*Generator_CostCurve_Matrix_New[1,1] * SolutionVector_p_Full_loop[1,1]) + (Generator_CostCurve_Matrix_New[1,2])) * Del_f_Del_x_Matrix_Delta_SlackBus
        Del_f_Del_x_Matrix_V = ((2*Generator_CostCurve_Matrix_New[1,1] * SolutionVector_p_Full_loop[1,1]) + (Generator_CostCurve_Matrix_New[1,2])) * Del_f_Del_x_Matrix_V_SlackBus

        # Computing Del_f_Del_x_Matrix
        Del_f_Del_x_Matrix = vcat(Del_f_Del_x_Matrix_Delta, Del_f_Del_x_Matrix_V) 


    elseif ((Len_SolutionVector_p_loop == Len_Generator_CostCurve_Matrix_New) && (Line_Index_Vector != nothing) && (Line_Bus_Index_Matrix != nothing) && (Line_P_Limit_Vector != nothing))  # No Slack Bus in Cost Function and Line Limits present

        # Computing Line Limit Violation indicator
        LineLimit_Violation_Indicator, Line_LimitViolated_Index_Vector, Line_LimitViolated_Bus_Index_Matrix, Line_LimitViolated_P_Limit_Vector = Compute_LineLimit_Violation_Indicator_OPF(LineFlow_Array, Line_Index_Vector, Line_Bus_Index_Matrix, Line_P_Limit_Vector)

        # Checking if Limits were violated
        if (LineLimit_Violation_Indicator == true)  # Line Limits are violated

            ## Compute Line Limits part of Del_f_Del_x_Matrix

            for ii in 1:size(Line_LimitViolated_Index_Vector)[1]  # For each element in Line_LimitViolated_Index_Vector

                # Initializing Del_f_Del_x_Matrix_Delta_Line and Del_f_Del_x_Matrix_V_Line
                Del_f_Del_x_Matrix_Delta_Line = zeros(Len_Delta, 1)
                Del_f_Del_x_Matrix_V_Line = zeros(Len_V, 1)

                # Getting desired Indices and Line Limit and Del_f_Del_Pij
                Line_Index = Int64(Line_LimitViolated_Index_Vector[ii,1])
                Bus_i_Index = Int64(Line_LimitViolated_Bus_Index_Matrix[ii,1])
                Bus_j_Index = Int64(Line_LimitViolated_Bus_Index_Matrix[ii,2])
                Line_Limit = Line_LimitViolated_P_Limit_Vector[ii,1]
                Del_f_Del_Pij = 2 * (LineFlow_Array[Line_Index,1] - Line_Limit)

                # Updating Del_f_Del_x_Matrix_Delta_Line w.r.t Delta_i
                if (Bus_i_Index != 1) # Not Slack Bus

                    Del_f_Del_x_Matrix_Delta_Line[Bus_i_Index-1,1] = ((SolutionVector_x_V_loop[Bus_i_Index,1] * SolutionVector_x_V_loop[Bus_j_Index,1] * abs(Ybus[Bus_i_Index,Bus_j_Index])) *(sin(angle(Ybus[Bus_i_Index,Bus_j_Index]) +  deg2rad(SolutionVector_x_Delta_loop[Bus_j_Index,1]) - deg2rad(SolutionVector_x_Delta_loop[Bus_i_Index,1]))))

                end

                # Updating Del_f_Del_x_Matrix_Delta_Line w.r.t Delta_j
                if (Bus_j_Index != 1) # Not Slack Bus

                    Del_f_Del_x_Matrix_Delta_Line[Bus_j_Index-1,1] = -((SolutionVector_x_V_loop[Bus_i_Index,1] * SolutionVector_x_V_loop[Bus_j_Index,1] * abs(Ybus[Bus_i_Index,Bus_j_Index])) *(sin(angle(Ybus[Bus_i_Index,Bus_j_Index]) +  deg2rad(SolutionVector_x_Delta_loop[Bus_j_Index,1]) - deg2rad(SolutionVector_x_Delta_loop[Bus_i_Index,1]))))

                end

                # Updating Del_f_Del_x_Matrix_Delta_Line with Del_f_Del_Pij
                Del_f_Del_x_Matrix_Delta_Line = Del_f_Del_Pij * Del_f_Del_x_Matrix_Delta_Line

                # Updating Del_f_Del_x_Matrix_V_Line w.r.t V_i
                if ((Bus_i_Index > 1) && (Bus_i_Index <= (N_PQ_Bus + 1)))  # Not a Slack or PV Bus

                    Del_f_Del_x_Matrix_V_Line[Bus_i_Index-1,1] = (2 *  SolutionVector_x_V_loop[Bus_i_Index,1] * real(Ybus[Bus_i_Index,Bus_j_Index])) + ((SolutionVector_x_V_loop[Bus_j_Index,1] * abs(Ybus[Bus_i_Index,Bus_j_Index])) *(cos(angle(Ybus[Bus_i_Index,Bus_j_Index]) +  deg2rad(SolutionVector_x_Delta_loop[Bus_j_Index,1]) - deg2rad(SolutionVector_x_Delta_loop[Bus_i_Index,1]))))

                end

                # Updating Del_f_Del_x_Matrix_V_Line w.r.t V_j
                if ((Bus_j_Index > 1) && (Bus_j_Index <= (N_PQ_Bus + 1)))  # Not a Slack or PV Bus

                    Del_f_Del_x_Matrix_V_Line[Bus_j_Index-1,1] = ((SolutionVector_x_V_loop[Bus_i_Index,1] * abs(Ybus[Bus_i_Index,Bus_j_Index])) *(cos(angle(Ybus[Bus_i_Index,Bus_j_Index]) +  deg2rad(SolutionVector_x_Delta_loop[Bus_j_Index,1]) - deg2rad(SolutionVector_x_Delta_loop[Bus_i_Index,1]))))

                end

                # Updating Del_f_Del_x_Matrix_V_Line with Del_f_Del_Pij
                Del_f_Del_x_Matrix_V_Line = Del_f_Del_Pij * Del_f_Del_x_Matrix_Delta_Line

                # Updating Del_f_Del_x_Matrix_Delta and Del_f_Del_x_Matrix_V_SlackBus
                Del_f_Del_x_Matrix_Delta = Del_f_Del_x_Matrix_Delta + Del_f_Del_x_Matrix_Delta_Line
                Del_f_Del_x_Matrix_V = Del_f_Del_x_Matrix_V + Del_f_Del_x_Matrix_V_Line

            end

            # Computing Del_f_Del_x_Matrix
            Del_f_Del_x_Matrix = vcat(Del_f_Del_x_Matrix_Delta, Del_f_Del_x_Matrix_V) 

        else  # Line Limits are not violated

            # Constructing Del_f_Del_x_Matrix
            Del_f_Del_x_Matrix = vcat(Del_f_Del_x_Matrix_Delta, Del_f_Del_x_Matrix_V)

        end

    elseif ((Len_SolutionVector_p_loop != Len_Generator_CostCurve_Matrix_New) && (Line_Index_Vector != nothing) && (Line_Bus_Index_Matrix != nothing) && (Line_P_Limit_Vector != nothing))  # Slack Bus in Cost Function and Line Limits present

        # Computing Line Limit Violation indicator
        LineLimit_Violation_Indicator, Line_LimitViolated_Index_Vector, Line_LimitViolated_Bus_Index_Matrix, Line_LimitViolated_P_Limit_Vector = Compute_LineLimit_Violation_Indicator_OPF(LineFlow_Array, Line_Index_Vector, Line_Bus_Index_Matrix, Line_P_Limit_Vector)

        # Checking if Limits were violated
        if (LineLimit_Violation_Indicator == true)  # Line Limits are violated

            ## Compute Slack Bus part of Del_f_Del_x_Matrix

            # Initializing Del_f_Del_x_Matrix_Delta_SlackBus and Del_f_Del_x_Matrix_V_SlackBus
            Del_f_Del_x_Matrix_Delta_SlackBus = zeros(Len_Delta, 1)
            Del_f_Del_x_Matrix_V_SlackBus = zeros(Len_V, 1)

            # Updating Del_f_Del_x_Matrix_Delta_SlackBus
            for ii in 1:Len_Delta  # For each element in Del_f_Del_x_Matrix_Delta_SlackBus

                Del_f_Del_x_Matrix_Delta_SlackBus[ii,1] = -((SolutionVector_x_V_loop[1,1] * SolutionVector_x_V_loop[ii+1,1] * abs(Ybus[1,ii+1])) *(sin(angle(Ybus[1,ii+1]) +  deg2rad(SolutionVector_x_Delta_loop[ii+1,1]) - deg2rad(SolutionVector_x_Delta_loop[1,1]))))

            end

            # Updating Del_f_Del_x_Matrix_V_SlackBus
            for ii in 1:Len_V  # For each element in Del_f_Del_x_Matrix_V_SlackBus

                Del_f_Del_x_Matrix_V_SlackBus[ii,1] = ((SolutionVector_x_V_loop[1,1] * abs(Ybus[1,ii+1])) *(cos(angle(Ybus[1,ii+1]) +  deg2rad(SolutionVector_x_Delta_loop[ii+1,1]) - deg2rad(SolutionVector_x_Delta_loop[1,1]))))

            end

            # Updating Del_f_Del_x_Matrix_Delta and Del_f_Del_x_Matrix_V_SlackBus
            Del_f_Del_x_Matrix_Delta = ((2*Generator_CostCurve_Matrix_New[1,1] * SolutionVector_p_Full_loop[1,1]) + (Generator_CostCurve_Matrix_New[1,2])) * Del_f_Del_x_Matrix_Delta_SlackBus
            Del_f_Del_x_Matrix_V = ((2*Generator_CostCurve_Matrix_New[1,1] * SolutionVector_p_Full_loop[1,1]) + (Generator_CostCurve_Matrix_New[1,2])) * Del_f_Del_x_Matrix_V_SlackBus

            ## Compute Line Limits part of Del_f_Del_x_Matrix
            for ii in 1:size(Line_LimitViolated_Index_Vector)[1]  # For each element in Line_LimitViolated_Index_Vector

                # Initializing Del_f_Del_x_Matrix_Delta_Line and Del_f_Del_x_Matrix_V_Line
                Del_f_Del_x_Matrix_Delta_Line = zeros(Len_Delta, 1)
                Del_f_Del_x_Matrix_V_Line = zeros(N_PQ_Bus, 1)

                # Getting desired Indices and Line Limit and Del_f_Del_Pij
                Line_Index = Int64(Line_LimitViolated_Index_Vector[ii,1])
                Bus_i_Index = Int64(Line_LimitViolated_Bus_Index_Matrix[ii,1])
                Bus_j_Index = Int64(Line_LimitViolated_Bus_Index_Matrix[ii,2])
                Line_Limit = Line_LimitViolated_P_Limit_Vector[ii,1]
                Del_f_Del_Pij = 2 * (LineFlow_Array[Line_Index,1] - Line_Limit)

                # Updating Del_f_Del_x_Matrix_Delta_Line w.r.t Delta_i
                if (Bus_i_Index != 1) # Not Slack Bus

                    Del_f_Del_x_Matrix_Delta_Line[Bus_i_Index-1,1] = ((SolutionVector_x_V_loop[Bus_i_Index,1] * SolutionVector_x_V_loop[Bus_j_Index,1] * abs(Ybus[Bus_i_Index,Bus_j_Index])) *(sin(angle(Ybus[Bus_i_Index,Bus_j_Index]) +  deg2rad(SolutionVector_x_Delta_loop[Bus_j_Index,1]) - deg2rad(SolutionVector_x_Delta_loop[Bus_i_Index,1]))))

                end

                # Updating Del_f_Del_x_Matrix_Delta_Line w.r.t Delta_j
                if (Bus_j_Index != 1) # Not Slack Bus

                    Del_f_Del_x_Matrix_Delta_Line[Bus_j_Index-1,1] = -((SolutionVector_x_V_loop[Bus_i_Index,1] * SolutionVector_x_V_loop[Bus_j_Index,1] * abs(Ybus[Bus_i_Index,Bus_j_Index])) *(sin(angle(Ybus[Bus_i_Index,Bus_j_Index]) +  deg2rad(SolutionVector_x_Delta_loop[Bus_j_Index,1]) - deg2rad(SolutionVector_x_Delta_loop[Bus_i_Index,1]))))

                end

                # Updating Del_f_Del_x_Matrix_Delta_Line with Del_f_Del_Pij
                Del_f_Del_x_Matrix_Delta_Line = Del_f_Del_Pij * Del_f_Del_x_Matrix_Delta_Line

                # Updating Del_f_Del_x_Matrix_V_Line w.r.t V_i
                if ((Bus_i_Index > 1) && (Bus_i_Index <= (N_PQ_Bus + 1)))  # Not a Slack or PV Bus

                    Del_f_Del_x_Matrix_V_Line[Bus_i_Index-1,1] = (2 *  SolutionVector_x_V_loop[Bus_i_Index,1] * real(Ybus[Bus_i_Index,Bus_j_Index])) + ((SolutionVector_x_V_loop[Bus_j_Index,1] * abs(Ybus[Bus_i_Index,Bus_j_Index])) *(cos(angle(Ybus[Bus_i_Index,Bus_j_Index]) +  deg2rad(SolutionVector_x_Delta_loop[Bus_j_Index,1]) - deg2rad(SolutionVector_x_Delta_loop[Bus_i_Index,1]))))

                end

                # Updating Del_f_Del_x_Matrix_V_Line w.r.t V_j
                if ((Bus_j_Index > 1) && (Bus_j_Index <= (N_PQ_Bus + 1)))  # Not a Slack or PV Bus

                    Del_f_Del_x_Matrix_V_Line[Bus_j_Index-1,1] = ((SolutionVector_x_V_loop[Bus_i_Index,1] * abs(Ybus[Bus_i_Index,Bus_j_Index])) *(cos(angle(Ybus[Bus_i_Index,Bus_j_Index]) +  deg2rad(SolutionVector_x_Delta_loop[Bus_j_Index,1]) - deg2rad(SolutionVector_x_Delta_loop[Bus_i_Index,1]))))

                end

                # Updating Del_f_Del_x_Matrix_V_Line with Del_f_Del_Pij
                Del_f_Del_x_Matrix_V_Line = Del_f_Del_Pij * Del_f_Del_x_Matrix_V_Line

                # Updating Del_f_Del_x_Matrix_Delta and Del_f_Del_x_Matrix_V_SlackBus
                Del_f_Del_x_Matrix_Delta = Del_f_Del_x_Matrix_Delta + Del_f_Del_x_Matrix_Delta_Line
                Del_f_Del_x_Matrix_V = Del_f_Del_x_Matrix_V + Del_f_Del_x_Matrix_V_Line

            end

            # Computing Del_f_Del_x_Matrix
            Del_f_Del_x_Matrix = vcat(Del_f_Del_x_Matrix_Delta, Del_f_Del_x_Matrix_V) 

        else  # Line Limits are not violated

            ## Compute Slack Bus part of Del_f_Del_x_Matrix

            # Initializing Del_f_Del_x_Matrix_Delta_SlackBus and Del_f_Del_x_Matrix_V_SlackBus
            Del_f_Del_x_Matrix_Delta_SlackBus = zeros(Len_Delta, 1)
            Del_f_Del_x_Matrix_V_SlackBus = zeros(N_PQ_Bus, 1)

            # Updating Del_f_Del_x_Matrix_Delta_SlackBus
            for ii in 1:Len_Delta  # For each element in Del_f_Del_x_Matrix_Delta_SlackBus

                Del_f_Del_x_Matrix_Delta_SlackBus[ii,1] = -((SolutionVector_x_V_loop[1,1] * SolutionVector_x_V_loop[ii+1,1] * abs(Ybus[1,ii+1])) *(sin(angle(Ybus[1,ii+1]) +  deg2rad(SolutionVector_x_Delta_loop[ii+1,1]) - deg2rad(SolutionVector_x_Delta_loop[1,1]))))

            end

            # Updating Del_f_Del_x_Matrix_V_SlackBus
            for ii in 1:Len_V  # For each element in Del_f_Del_x_Matrix_V_SlackBus

                Del_f_Del_x_Matrix_v_SlackBus[ii,1] = ((SolutionVector_x_V_loop[1,1] * abs(Ybus[1,ii+1])) *(cos(angle(Ybus[1,ii+1]) +  deg2rad(SolutionVector_x_Delta_loop[ii+1,1]) - deg2rad(SolutionVector_x_Delta_loop[1,1]))))

            end

            # Updating Del_f_Del_x_Matrix_Delta and Del_f_Del_x_Matrix_V_SlackBus
            Del_f_Del_x_Matrix_Delta = ((2*Generator_CostCurve_Matrix_New[1,1] * SolutionVector_p_Full_loop[1,1]) + (Generator_CostCurve_Matrix_New[1,2])) * Del_f_Del_x_Matrix_Delta_SlackBus
            Del_f_Del_x_Matrix_V = ((2*Generator_CostCurve_Matrix_New[1,1] * SolutionVector_p_Full_loop[1,1]) + (Generator_CostCurve_Matrix_New[1,2])) * Del_f_Del_x_Matrix_V_SlackBus

            # Computing Del_f_Del_x_Matrix
            Del_f_Del_x_Matrix = vcat(Del_f_Del_x_Matrix_Delta, Del_f_Del_x_Matrix_V) 

        end

    end

    # Addressing Machine Precision
    Del_f_Del_x_Matrix = Addressing_MachinePrecision(Del_f_Del_x_Matrix, 1e-12)

    return Del_f_Del_x_Matrix 

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

