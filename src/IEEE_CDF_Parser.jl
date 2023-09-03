# IEEE_CDF_Parser.jl

#Uncomment next lines if you want to test the parser here itself.
# using CSV
# using DataFrames
# include("Ybus_Builder.jl")

"""
    CDF_Parser(CDF_FilePath)

Creates Julia Dataframe from the IEEE Common Data Format (CDF) text file.

'''
# Arguments
- 'CDF_FilePath': File path to the IEEE CDF text file.
- 'SortValue': 1 - Sort CDF File according to Bus Type PQ->PV->Slack, any other value - do not sort
'''
'''
# Output
- 'CDF_DF_List': IEEE CDF file in List of Dataframe format according to
Data Card types in IEEE CDF file : [TitleCard_DF, BusDataCard_DF,
BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF,
TieLinesDataCard_DF].
'''
"""
function CDF_Parser(CDF_FilePath,
        SortValue;
        saveTables::Bool=false, 
        saveLocation = "processedData/")
        
        # Read lines of IEEE CDF File in an Array
        CDF_Text_File = open(CDF_FilePath)

        CDF_Text_Array = readlines(CDF_Text_File)

        # Gathering Data of Different Cards

        # Gathering Data From Title Card into a DataFrame
        TitleCard_DF = DataFrame(DateTime = CDF_Text_Array[1][2:9],
                                 Origin = CDF_Text_Array[1][11:30],
                                 MVA_Base = parse(Float64, CDF_Text_Array[1][32:37]),
                                 Year = parse(Int64, CDF_Text_Array[1][39:42]),
                                 Season = CDF_Text_Array[1][44],
                                 Case_ID = CDF_Text_Array[1][46:end])

        # Initializing Data Card Arrays
        BusDataCard_Array=[]
        BranchDataCard_Array=[]
        LossZonesCard_Array=[]
        InterchangeDataCard_Array=[]
        TieLinesCard_Array=[]

        # Gathering Data from Data Cards into a Array
        for ii in 1:length(CDF_Text_Array)

                # Creating Regular Expressions for matching CDF File Card Titles
                BusCard_Title_RegEx = match(r"BUS DATA FOLLOWS",CDF_Text_Array[ii])
                if (isnothing(BusCard_Title_RegEx))
                        BusCard_Title_RegEx = match(r"",CDF_Text_Array[ii])
                end
                BranchCard_Title_RegEx = match(r"BRANCH DATA FOLLOWS",CDF_Text_Array[ii])
                if (isnothing(BranchCard_Title_RegEx))
                        BranchCard_Title_RegEx = match(r"",CDF_Text_Array[ii])
                end
                ZoneCard_Title_RegEx = match(r"LOSS ZONES FOLLOWS",CDF_Text_Array[ii])
                if (isnothing(ZoneCard_Title_RegEx))
                        ZoneCard_Title_RegEx = match(r"",CDF_Text_Array[ii])
                end
                InterchangeCard_Title_RegEx = match(r"INTERCHANGE DATA FOLLOWS",CDF_Text_Array[ii])
                if (isnothing(InterchangeCard_Title_RegEx))
                        InterchangeCard_Title_RegEx = match(r"",CDF_Text_Array[ii])
                end
                TieLinesCard_Title_RegEx = match(r"TIE LINES FOLLOWS",CDF_Text_Array[ii])
                if (isnothing(TieLinesCard_Title_RegEx))
                        TieLinesCard_Title_RegEx = match(r"",CDF_Text_Array[ii])
                end

                # Gathering Data from Bus Data Card into a Array
                if (!Bool(cmp("BUS DATA FOLLOWS",BusCard_Title_RegEx.match)))


                        for jj in ii+1:length(CDF_Text_Array)

                                End_Card_RegEx = match(r"-999",CDF_Text_Array[jj])
                                if (isnothing(End_Card_RegEx))
                                        End_Card_RegEx = match(r"",CDF_Text_Array[jj])
                                end 

                                if (!Bool(cmp("-999",End_Card_RegEx.match)))

                                        break

                                else

                                        append!(BusDataCard_Array,[CDF_Text_Array[jj]])

                                end

                        end

                # Gathering Data from Branch Data Card into a Array
                elseif ((!Bool(cmp("BRANCH DATA FOLLOWS",BranchCard_Title_RegEx.match))))

                        for jj in ii+1:length(CDF_Text_Array)

                                End_Card_RegEx = match(r"-999",CDF_Text_Array[jj])
                                if (isnothing(End_Card_RegEx))
                                        End_Card_RegEx = match(r"",CDF_Text_Array[jj])
                                end 

                                if (!Bool(cmp("-999",End_Card_RegEx.match)))

                                        break

                                else

                                        append!(BranchDataCard_Array,[CDF_Text_Array[jj]])

                                end

                        end

                # Gathering Data from Line Zones Card into a Array
                elseif ((!Bool(cmp("LOSS ZONES FOLLOWS",ZoneCard_Title_RegEx.match))))

                        for jj in ii+1:length(CDF_Text_Array)

                                End_Card_RegEx = match(r"-99",CDF_Text_Array[jj])
                                if (isnothing(End_Card_RegEx))
                                        End_Card_RegEx = match(r"",CDF_Text_Array[jj])
                                end 

                                if (!Bool(cmp("-99",End_Card_RegEx.match)))

                                        break

                                else

                                        append!(LossZonesCard_Array,[CDF_Text_Array[jj]])

                                end

                        end

                # Gathering Data from Interchange Data Card into a Array
                elseif ((!Bool(cmp("INTERCHANGE DATA FOLLOWS",InterchangeCard_Title_RegEx.match))))

                        for jj in ii+1:length(CDF_Text_Array)

                                End_Card_RegEx = match(r"-9",CDF_Text_Array[jj])
                                if (isnothing(End_Card_RegEx))
                                        End_Card_RegEx = match(r"",CDF_Text_Array[jj])
                                end 

                                if (!Bool(cmp("-9",End_Card_RegEx.match)))

                                        break

                                else

                                        append!(InterchangeDataCard_Array,[CDF_Text_Array[jj]])

                                end

                        end

                # Gathering Data from Tie Lines Data Card into a Array
                elseif ((!Bool(cmp("TIE LINES FOLLOWS",TieLinesCard_Title_RegEx.match))))

                        for jj in ii+1:length(CDF_Text_Array)

                                End_Card_RegEx = match(r"-999",CDF_Text_Array[jj])
                                if (isnothing(End_Card_RegEx))
                                        End_Card_RegEx = match(r"",CDF_Text_Array[jj])
                                end 

                                if (!Bool(cmp("-999",End_Card_RegEx.match)))

                                        break

                                else

                                        append!(TieLinesCard_Array,[CDF_Text_Array[jj]])

                                end

                        end

                else

                        continue

                end

        end

        # Initializing Data DataFrames
        BusDataCard_DF=DataFrame(Bus_Num = Int64[],
                                 Name = String[],
                                 LF_Area_Num = Int64[],
                                 LZ_Num = Int64[],
                                 Type = Int64[],
                                 Type_Original = Int64[],
                                 Final_V_pu = Float64[],
                                 Final_A_deg = Float64[],
                                 Final_V_pu_Original = Float64[],
                                 Final_A_deg_Original = Float64[],
                                 Load_MW = Float64[],
                                 Load_MVAR = Float64[],
                                 Gen_MW = Float64[],
                                 Gen_MVAR = Float64[],
                                 Base_KV = Float64[],
                                 Desired_V_pu = Float64[],
                                 Max_MVAR_V = Float64[],
                                 Min_MVAR_V = Float64[],
                                 MVAR_V_Limit = Float64[],
                                 G_pu = Float64[],
                                 B_pu = Float64[],
                                 Remote_Con_Bus_Num = Int64[],
                                 V_Measured = Float64[],
                                 P_Measured = Float64[],
                                 Q_Measured = Float64[],
                                 ErrorVariance_V = Float64[],
                                 ErrorVariance_P = Float64[],
                                 ErrorVariance_Q = Float64[])

        BranchDataCard_DF=DataFrame(Tap_Bus_Num = Int64[],
                                    Z_Bus_Num = Int64[],
                                    LF_Area = Int64[],
                                    LZ = Int64[],
                                    Circuit = Int64[],
                                    Type = Int64[],
                                    R_pu = Float64[],
                                    X_pu =Float64[],
                                    B_pu = Float64[],
                                    MVA_Rating_1 = Int64[],
                                    MVA_Rating_2 = Int64[],
                                    MVA_Rating_3 = Int64[],
                                    Con_Bus_Num = Int64[],
                                    Side = Int64[],
                                    Transformer_t = Float64[],
                                    Transformer_ps = Float64[],
                                    Min_t_ps = Float64[],
                                    Max_t_ps = Float64[],
                                    StepSize = Float64[],
                                    Min_MVAR_MW_V = Float64[],
                                    Max_MVAR_MW_V = Float64[],
                                    Line_Flow_P = Float64[],
                                    Line_Flow_Q = Float64[],
                                    Line_Flow_Pos_P = Float64[],
                                    Line_Flow_Neg_P = Float64[],
                                    Line_Flow_Pos_Q = Float64[],
                                    Line_Flow_Neg_Q = Float64[],
                                    ErrorVariance_Pos_P = Float64[],
                                    ErrorVariance_Neg_P = Float64[],
                                    ErrorVariance_Pos_Q = Float64[],
                                    ErrorVariance_Neg_Q = Float64[])

        LossZonesCard_DF=DataFrame(LZ_Num = Int64[],
                                   LZ_Name = String[])

        InterchangeDataCard_DF=DataFrame(Area_Num = Int64[],
                                         Inter_SB_Num = Int64[],
                                         Alt_SB_Name = String[],
                                         Area_Inter_export = Float64[],
                                         Area_Inter_tol = Float64[],
                                         Area_Code = String[],
                                         Area_Name = String[])

        TieLinesDataCard_DF=DataFrame(Meter_Bus_Num = Int64[],
                                  Meter_Area_Num = Int64[],
                                  NonMeter_Bus_Num = Int64[],
                                  NonMeter_Area_Num = Int64[],
                                  Circuit_Num = Int64[])

        # Filling up BusDataCard_DF
        for ii in 1:length(BusDataCard_Array)

                push!(BusDataCard_DF,(parse(Int64,BusDataCard_Array[ii][1:4]),
                                      BusDataCard_Array[ii][6:17],
                                      parse(Int64,BusDataCard_Array[ii][19:20]),
                                      parse(Int64,BusDataCard_Array[ii][21:23]),
                                      parse(Int64,BusDataCard_Array[ii][25:26]),
                                      parse(Int64,BusDataCard_Array[ii][25:26]),
                                      parse(Float64,BusDataCard_Array[ii][28:33]),
                                      parse(Float64,BusDataCard_Array[ii][34:40]),
                                      parse(Float64,BusDataCard_Array[ii][28:33]),
                                      parse(Float64,BusDataCard_Array[ii][34:40]),
                                      parse(Float64,BusDataCard_Array[ii][41:49]),
                                      parse(Float64,BusDataCard_Array[ii][50:59]),
                                      parse(Float64,BusDataCard_Array[ii][60:67]),
                                      parse(Float64,BusDataCard_Array[ii][68:75]),
                                      parse(Float64,BusDataCard_Array[ii][77:83]),
                                      parse(Float64,BusDataCard_Array[ii][85:90]),
                                      parse(Float64,BusDataCard_Array[ii][91:98]),
                                      parse(Float64,BusDataCard_Array[ii][99:106]),
                                      -9999,
                                      parse(Float64,BusDataCard_Array[ii][107:114]),
                                      parse(Float64,BusDataCard_Array[ii][115:122]),
                                      parse(Int64,BusDataCard_Array[ii][124:end]),
                                      -9999,
                                      -9999,
                                      -9999,
                                      -9999,
                                      -9999,
                                      -9999))

        end

        # Ordering BusDataCard_DF: PQ->PV->Slack
        if (SortValue == 1)
        
                sort!(BusDataCard_DF, [order(:Type)])

        end


        # Filling up BranchDataCard_DF
        for ii in 1:length(BranchDataCard_Array)

                push!(BranchDataCard_DF,(parse(Int64,BranchDataCard_Array[ii][1:4]),
                                      parse(Int64,BranchDataCard_Array[ii][6:9]),
                                      parse(Int64,BranchDataCard_Array[ii][11:12]),
                                      parse(Int64,BranchDataCard_Array[ii][13:15]),
                                      parse(Int64,BranchDataCard_Array[ii][17]),
                                      parse(Int64,BranchDataCard_Array[ii][19]),
                                      parse(Float64,BranchDataCard_Array[ii][20:29]),
                                      parse(Float64,BranchDataCard_Array[ii][30:40]),
                                      parse(Float64,BranchDataCard_Array[ii][41:50]),
                                      parse(Int64,BranchDataCard_Array[ii][51:55]),
                                      parse(Int64,BranchDataCard_Array[ii][57:61]),
                                      parse(Int64,BranchDataCard_Array[ii][63:67]),
                                      parse(Int64,BranchDataCard_Array[ii][69:72]),
                                      parse(Int64,BranchDataCard_Array[ii][74]),
                                      parse(Float64,BranchDataCard_Array[ii][77:82]),
                                      parse(Float64,BranchDataCard_Array[ii][84:90]),
                                      parse(Float64,BranchDataCard_Array[ii][91:97]),
                                      parse(Float64,BranchDataCard_Array[ii][98:104]),
                                      parse(Float64,BranchDataCard_Array[ii][106:111]),
                                      parse(Float64,BranchDataCard_Array[ii][113:118]),
                                      parse(Float64,BranchDataCard_Array[ii][120:end]),
                                      0.0,
                                      0.0,
                                      -9999,
                                      -9999,
                                      -9999,
                                      -9999,
                                      -9999,
                                      -9999,
                                      -9999,
                                      -9999))

        end

        # Filling up LossZonesCard_DF
        for ii in 1:length(LossZonesCard_Array)

                push!(LossZonesCard_DF,(parse(Int64,LossZonesCard_Array[ii][1:3]),
                                      LossZonesCard_Array[ii][5:end]))

        end


        # Filling up InterchangeDataCard_DF
        for ii in 1:length(InterchangeDataCard_Array)

                push!(InterchangeDataCard_DF,(parse(Int64,InterchangeDataCard_Array[ii][1:2]),
                                      parse(Int64,InterchangeDataCard_Array[ii][4:7]),
                                      InterchangeDataCard_Array[ii][9:20],
                                      parse(Float64,InterchangeDataCard_Array[ii][21:28]),
                                      parse(Float64,InterchangeDataCard_Array[ii][30:35]),
                                      InterchangeDataCard_Array[ii][38:43],
                                      InterchangeDataCard_Array[ii][46:end]))

        end

        # Filling up TieLinesDataCard_DF
        for ii in 1:length(TieLinesCard_Array)

                push!(TieLinesDataCard_DF,(parse(Int64,TieLinesCard_Array[ii][1:4]),
                                      parse(Int64,TieLinesCard_Array[ii][7:8]),
                                      parse(Int64,TieLinesCard_Array[ii][11:14]),
                                      parse(Int64,TieLinesCard_Array[ii][17:18]),
                                      parse(Int64,TieLinesCard_Array[ii][21])))

        end

        CDF_DF_List = [TitleCard_DF, BusDataCard_DF, BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF, TieLinesDataCard_DF]
        systemName = extractSystemName(CDF_DF_List)

        filenames = ["TitleCard.csv", "BusDataCard_MVA.csv", "BranchDataCard_MVA.csv", "LossZonesCard_MVA.csv", "InterchangeDataCard_MVA.csv", "TieLinesDataCard_MVA.csv"]

        if saveTables
                for (df, filename) in zip(CDF_DF_List, filenames)
                        CSV.write(saveLocation*systemName*"/"*filename, df)
                end
        end

        return CDF_DF_List

end

"""
    CDF_pu_Converter(CDF_DF_List; saveTables::Bool=false, saveLocation="processedData/")

Converts the actual value columns of a Common Data Format (CDF) DataFrame list to per unit (pu) values.

## Arguments
- `CDF_DF_List`: A list of DataFrames representing the CDF data. The list should contain the following DataFrames in the specified order:
    - `TitleCard_DF`: DataFrame representing the title card data.
    - `BusDataCard_DF`: DataFrame representing the bus data card.
    - `BranchDataCard_DF`: DataFrame representing the branch data card.
    - `LossZonesCard_DF`: DataFrame representing the loss zones card.
    - `InterchangeDataCard_DF`: DataFrame representing the interchange data card.
    - `TieLinesDataCard_DF`: DataFrame representing the tie lines data card.
- `saveTables::Bool` (optional, default=false): A flag indicating whether to save the converted DataFrames as CSV files.
- `saveLocation` (optional, default="processedData/"): The directory path where the converted CSV files will be saved.

## Returns
- `CDF_DF_List_pu`: A list of DataFrames with the actual value columns converted to per unit (pu) values.

## Note
This function assumes that the input DataFrames have specific column names and structures. Make sure the input DataFrames match the expected format.

## Example
```julia
# Assuming you have loaded the CDF data into the CDF_DF_List

# Convert the actual value columns to per unit (pu) values
CDF_DF_List_pu = CDF_pu_Converter(CDF_DF_List, saveTables=true, saveLocation="processedData/")

# CDF_DF_List_pu contains the converted DataFrames
# The converted DataFrames are also saved as CSV files in the "processedData" directory
"""
function CDF_pu_Converter(CDF_DF_List;
        saveTables::Bool=false, 
        saveLocation = "processedData/")

        # Getting required data from CDF_DF_List
        # WARNING: Be very very careful when copying arrays in Julia
        # Here without deepcopy, the original arguments were modified too,
        # leading to values like Gen_MW = 0.02324 in CDF_DF_List[2]
        TitleCard_DF = deepcopy(CDF_DF_List[1])
        BusDataCard_DF = deepcopy(CDF_DF_List[2])
        BranchDataCard_DF = deepcopy(CDF_DF_List[3])
        LossZonesCard_DF = deepcopy(CDF_DF_List[4])
        InterchangeDataCard_DF = deepcopy(CDF_DF_List[5])
        TieLinesDataCard_DF = deepcopy(CDF_DF_List[6])

        systemName = extractSystemName(CDF_DF_List)

        # Getting Base MVA
        Base_MVA = TitleCard_DF.MVA_Base[1]

        # Converting Actual Value Columns in BusDataCard_DF to pu
        for ii in 1:length(BusDataCard_DF.Bus_Num)

                # Converting Load MW/MVAR and Gen MW/MVAR to pu
                BusDataCard_DF.Load_MW[ii] /= Base_MVA

                BusDataCard_DF.Load_MVAR[ii] /= Base_MVA

                BusDataCard_DF.Gen_MW[ii] /= Base_MVA

                BusDataCard_DF.Gen_MVAR[ii] /= Base_MVA

                # Converting Max_MVAR_V and Min_MVAR_V to pu
                if (BusDataCard_DF.Type[ii] == 0) # Unregulated PQ

                        # Do Nothing

                elseif (BusDataCard_DF.Type[ii] == 1) # Voltage Limit PQ

                        BusDataCard_DF.Max_MVAR_V[ii] = BusDataCard_DF.Max_MVAR_V[ii]/BusDataCard_DF.Base_KV[ii]

                        BusDataCard_DF.Min_MVAR_V[ii] = BusDataCard_DF.Min_MVAR_V[ii]/BusDataCard_DF.Base_KV[ii]

                elseif (BusDataCard_DF.Type[ii] == 2) # VAR Limit PV

                        BusDataCard_DF.Max_MVAR_V[ii] = BusDataCard_DF.Max_MVAR_V[ii]/Base_MVA

                        BusDataCard_DF.Min_MVAR_V[ii] = BusDataCard_DF.Min_MVAR_V[ii]/Base_MVA

                elseif (BusDataCard_DF.Type[ii] == 3) # Slack

                        # Do Nothing

                end

        end

        # Converting Actual Value Columns in BranchDataCard_DF to pu
        for ii in 1:length(BranchDataCard_DF.Tap_Bus_Num)

                # Converting Max_MVAR_MW_V and Min_MVAR_MW_V to pu
                if (BranchDataCard_DF.Type[ii] == 0) # Transmission Line

                        # Do Nothing

                elseif (BranchDataCard_DF.Type[ii] == 1) # Fixed Tap

                        # Do Nothing

                elseif (BranchDataCard_DF.Type[ii] == 2) # Variable Tap V Control

                        # Getting Base KV from BusDataCard_DF
                        Base_KV = filter(row -> row.Bus_Num == BranchDataCard_DF.Tap_Bus_Num[ii], BusDataCard_DF).Base_KV

                        # Converting to pu
                        BranchDataCard_DF.Max_MVAR_MW_V[ii] = BranchDataCard_DF.Max_MVAR_MW_V[ii]/Base_KV

                        BranchDataCard_DF.Min_MVAR_MW_V[ii] = BranchDataCard_DF.Min_MVAR_MW_V[ii]/Base_KV

                elseif (BranchDataCard_DF.Type[ii] == 3) # Variable Tap MVAR Control

                        BranchDataCard_DF.Max_MVAR_MW_V[ii] = BranchDataCard_DF.Max_MVAR_MW_V[ii]/Base_MVA

                        BranchDataCard_DF.Min_MVAR_MW_V[ii] = BranchDataCard_DF.Min_MVAR_MW_V[ii]/Base_MVA

                elseif (BranchDataCard_DF.Type[ii] == 4) # Variable Phase Angle MW Control

                        BranchDataCard_DF.Max_MVAR_MW_V[ii] = BranchDataCard_DF.Max_MVAR_MW_V[ii]/Base_MVA

                        BranchDataCard_DF.Min_MVAR_MW_V[ii] = BranchDataCard_DF.Min_MVAR_MW_V[ii]/Base_MVA

                end

        end

        # Converting Actual Value Columns in InterchangeDataCard_DF to pu
        for ii in 1:length(InterchangeDataCard_DF.Area_Num)

                # Converting Area_Inter_export and Area_Inter_tol to pu
                InterchangeDataCard_DF.Area_Inter_export[ii] = InterchangeDataCard_DF.Area_Inter_export[ii]/Base_MVA

                InterchangeDataCard_DF.Area_Inter_tol[ii] = InterchangeDataCard_DF.Area_Inter_tol[ii]/Base_MVA

        end

        CDF_DF_List_pu = [TitleCard_DF, BusDataCard_DF, BranchDataCard_DF, LossZonesCard_DF, InterchangeDataCard_DF, TieLinesDataCard_DF]

        filenames = ["TitleCard.csv", "BusDataCard_pu.csv", "BranchDataCard_pu.csv", "LossZonesCard_pu.csv", "InterchangeDataCard_pu.csv", "TieLinesDataCard_pu.csv"]

        if saveTables
                for (df, filename) in zip(CDF_DF_List_pu, filenames)
                        CSV.write(saveLocation*systemName*"/"*filename, df)
                end
        end

        return CDF_DF_List_pu
end

# # test the parser
# CDF_DF_List = CDF_Parser("data/IEEE_14/IEEE_14_Data.txt")
# busData = CDF_DF_List[2]
# busData.Gen_MW

# CDF_DF_List_pu = CDF_pu_Converter(CDF_DF_List)
# busData_pu = CDF_DF_List_pu[2]
# busData_pu.Gen_MW