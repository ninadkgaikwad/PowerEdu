# Ybus_Builder.jl

"""
    Create_Ybus_WithoutTaps(CDF_DF_List)

Creates Ybus without taps for a power system network.

'''
# Arguments
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
'''
'''
# Output
- 'Ybus_WithoutTaps': A complex array of Ybus elements ordered according to bus
type: Slack->PQ->PV.
'''
"""
function Create_Ybus_WithoutTaps(CDF_DF_List_pu)

    # Getting required data from CDF_DF_List
    BusDataCard_DF = CDF_DF_List_pu[2]
    BranchDataCard_DF = CDF_DF_List_pu[3]

    # Getting Size of Ybus
    Size_Ybus = length(BusDataCard_DF.Bus_Num)

    # Initializing Ybus Complex Array
    Ybus_WithoutTaps = Array{Complex{Float64}}(undef, Size_Ybus,Size_Ybus)

    # Computing Ybus Off-Diagonal elements
    for ii in 1:Size_Ybus # Through Rows

        for jj in 1:1:Size_Ybus # Through Columns

            if (ii == jj) # Diagonal Element

                continue

            elseif (ii < jj) # Off-Diagonal elements upper triangle

                # Getting currentBus Numbers from BusDataCard_DF
                Bus1_Num = BusDataCard_DF.Bus_Num[ii]

                Bus2_Num = BusDataCard_DF.Bus_Num[jj]

                # Finding Row in BranchDataCard_DF based on current Bus Numbers
                BranchDataCard_FilterRow = filter(row -> ((row.Tap_Bus_Num == Bus1_Num) && (row.Z_Bus_Num == Bus2_Num)) || ((row.Tap_Bus_Num == Bus2_Num) && (row.Z_Bus_Num == Bus1_Num)), BranchDataCard_Row)

                BranchDataCard_FilterRow_Num = nrow(BranchDataCard_FilterRow)

                if (BranchDataCard_FilterRow_Num == 0) # There is no connection between buses

                    # Filling up Branch Admittance in Ybus_WithoutTaps
                    Ybus_WithoutTaps[ii,jj] = complex(0,0)

                    # Ybus is Symmetrical
                    Ybus_WithoutTaps[jj,ii] = Ybus_WithoutTaps[ii,jj]

                elseif (BranchDataCard_FilterRow_Num > 0) # There is connection between buses

                    # Creating Line Series Admittance
                    Line_SeriesAdmittance = 1/complex(BranchDataCard_FilterRow.R_pu[1],BranchDataCard_FilterRow.X_pu[1])

                    # Filling up Branch Admittance in Ybus_WithoutTaps
                    Ybus_WithoutTaps[ii,jj] = -Line_SeriesAdmittance

                    # Ybus is Symmetrical
                    Ybus_WithoutTaps[jj,ii] = Ybus_WithoutTaps[ii,jj]

                end


            else # Off-Diagonal elements lower triangle

                continue

            end

        end

    end

    # Computing Ybus Diagonal elements
    for ii in 1:Size_Ybus # Through Diagonal Elements Row-wise

        # Getting effect of Off-Diagonal Terms
        OffDiagonal_Terms = 0

        for jj = 1:Size_Ybus # Thorough Columns

            if (ii == jj) # Diagonal Term

                continue

            else

                OffDiagonal_Terms = OffDiagonal_Terms + (-Ybus_WithoutTaps[ii,jj])

            end

        end

        # Getting Effect of Bus Shunt Admittance
        BusAdmittance_Shunt = complex(BusDataCard_DF.G_pu[ii],-BusDataCard_DF.B_pu[ii])

        # Getting Effect of Line Shunt Admittance connected to the Bus
        BusLineAdmittance_Shunt = 0

        Bus_Num = BusDataCard_DF.Bus_Num[ii]

        BranchDataCard_Filter = filter(row -> (row.Tap_Bus_Num == Bus_Num) || (row.Z_Bus_Num == Bus_Num), BranchDataCard_Row)

        BranchDataCard_Filter_Num = nrow(BranchDataCard_Filter)

        if (BranchDataCard_Filter_Num == 0) # Bus not connected to any other bus through a transmission line

            BusLineAdmittance_Shunt = 0

        elseif (BranchDataCard_Filter_Num > 0) # Bus connected to other buses through a transmission lines

            for kk in 1:length(BranchDataCard_Filter.Tap_Bus_Num)

                BusLineAdmittance_Shunt = BusLineAdmittance_Shunt + complex(0,-(BranchDataCard_Filter.B_pu[kk]/2))

            end

        end

        # Total effect oin Ybus Diagonal term
        Ybus_WithoutTaps[ii,ii] = OffDiagonal_Terms + BusAdmittance_Shunt + BusLineAdmittance_Shunt

    end

    # Rearranging Create_Ybus_WithoutTaps in the order Slack->PQ->PV
    Ybus_WithoutTaps_PQ_PV = Ybus_WithoutTaps[1:end-1,1:end-1]

    Ybus_WithoutTaps_Slack1 = Ybus_WithoutTaps[1:end-1,end]

    Ybus_WithoutTaps_Slack2 = Ybus_WithoutTaps[end,1:end-1]

    Ybus_WithoutTaps_Slack3 = Ybus_WithoutTaps[end,end]

    Ybus_WithoutTaps_1 = vcat(Ybus_WithoutTaps_Slack2,Ybus_WithoutTaps_PQ_PV)

    Ybus_WithoutTaps_2 = vcat(Ybus_WithoutTaps_Slack3,Ybus_WithoutTaps_Slack1)

    Ybus_WithoutTaps = hcat(Ybus_WithoutTaps_1,Ybus_WithoutTaps_2)

    return Ybus_WithoutTaps

end

"""
    Create_Ybus_WithoutTaps(CDF_DF_List)

Creates Ybus with taps for a power system network.

'''
# Arguments
- 'Ybus_WithoutTaps': A complex array of Ybus elements ordered according to bus
type: Slack->PQ->PV.
- 'CDF_DF_List_pu': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [BusDataCard_DF, BranchDataCard_DF,
LossZonesCard_DF, InterchangeDataCard_DF, TieLinesDataCard_DF].
'''
'''
# Output
- 'Ybus_WithTaps': A complex array of Ybus elements ordered according to bus
type: Slack->PV->PQ.
'''
"""
function Create_Ybus_WithTaps(Ybus_WithoutTaps,CDF_DF_List_pu)

    # Getting required data from CDF_DF_List
    BusDataCard_DF = CDF_DF_List_pu[2]

    BranchDataCard_DF = CDF_DF_List_pu[3]

    # Getting Last Row number of BusDataCard_DF to locate Slack Bus Row number
    SlackBus_RowNumber = length(BusDataCard_DF.Bus_Num)

    # Initializing Ybus_WithTaps
    Ybus_WithTaps = Ybus_WithoutTaps

    # Getting Subset of BranchDataCard_DFfor lines with Tap Changing Transformers
    BranchDataCard_Filter = filter(row -> ((row.Transformer_t != 0) || (row.Transformer_ps != 0), BranchDataCard_DF))

    BranchDataCard_Filter_Num = nrow(BranchDataCard_Filter)

    if (BranchDataCard_Filter_Num == 0) # No Tap Changing Transformers present

        Ybus_WithTaps = Ybus_WithoutTaps

    elseif (BranchDataCard_Filter_Num > 0) # Tap Changing Transformers present

        for ii in 1:BranchDataCard_Filter_Num

            # Creating Tap value 'a'
            a = BranchDataCard_Filter.Transformer_t[ii] * cis(deg2rad(BranchDataCard_Filter.Transformer_ps[ii]))

            # Getting Bus Numbers 'i': Z_Bus_Num (Impedance Side) , 'j': Tap_Bus_Num (non-unity tap Side)
            Bus_Num_i = BranchDataCard_Filter.Z_Bus_Num[ii]

            Bus_Num_j = BranchDataCard_Filter.Tap_Bus_Num[ii]

            # Getting associated 'Bus_i_Index' and 'Bus_j_Index' from BusDataCard_DF to access correct location within Ybus_WithoutTaps
            for jj in 1:SlackBus_RowNumber

                if (Bus_Num_i == BusDataCard_DF.Bus_Num[jj])

                    if (jj == SlackBus_RowNumber)

                        Bus_i_Index = 1

                    else

                        Bus_i_Index = jj+1

                    end

                elseif (Bus_Num_j == BusDataCard_DF.Bus_Num[jj])

                    if (jj == SlackBus_RowNumber)

                        Bus_j_Index = 1

                    else

                        Bus_j_Index = jj+1

                    end

                else

                    continue

                end

            end

            # Changing the [Bus_i_Index, Bus_j_Index] in Ybus_WithTaps based on 'a'

            # Changing [Bus_i_Index, Bus_i_Index]
            Ybus_WithTaps[Bus_i_Index, Bus_i_Index] = Ybus_WithoutTaps[Bus_i_Index, Bus_i_Index]

            # Changing [Bus_i_Index, Bus_j_Index]
            Ybus_WithTaps[Bus_i_Index, Bus_j_Index] = Ybus_WithoutTaps[Bus_i_Index, Bus_j_Index]/a

            # Changing [Bus_j_Index, Bus_i_Index]
            Ybus_WithTaps[Bus_j_Index, Bus_i_Index] = Ybus_WithoutTaps[Bus_j_Index, Bus_i_Index]/conj(a)

            # Changing [Bus_j_Index, Bus_j_Index]
            Ybus_WithTaps[Bus_j_Index, Bus_j_Index] = (Ybus_WithoutTaps[Bus_j_Index, Bus_j_Index]) - (-Ybus_WithoutTaps[Bus_j_Index, Bus_i_Index]) + (-Ybus_WithoutTaps[Bus_j_Index, Bus_i_Index]/abs2(a))

        end

    end

    return Ybus_WithTaps

end
