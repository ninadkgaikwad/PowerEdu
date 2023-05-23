using CSV
using DelimitedFiles, DataFrames

folderInput = "rawData/"
folder_processedData = "processedData/"
systemName = "ieee14"

fileType_busData = ".csv"
filename_busData = folder_processedData*systemName*"/busData"*fileType_busData
println(filename_busData)

fileType_branchData = ".csv"
filename_branchData = folder_processedData*systemName*"/branchData"*fileType_branchData
println(filename_branchData)

# Load busData into a DataFrame
busData = CSV.read(filename_busData, DataFrame)

# Load branchData into a DataFrame
branchData = CSV.read(filename_branchData, DataFrame)

# vscodedisplay(busData)
# vscodedisplay(branchData)

ybus, BMatrix, b, A, branchNames, E = ybusGenerator(busData, branchData, saveTables=true, systemName=systemName)
ybusByTypes, BMatrixByTypes = ybusGenerator(busData, branchData, sortBy="busTypes", verbose=true, saveTables=true, saveLocation=folder_processedData, systemName=systemName)   