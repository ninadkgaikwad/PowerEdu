function sortMatrixByBusTypes(busData, ybus)

    @show outputs  = initializeVectors(busData)
    @show listOfSlackBuses, listOfPVBuses, listOfPQBuses = outputs[5], outputs[6], outputs[7]
    # @show typeof(listOfPVBuses)
    @show newOrder = vcat(listOfSlackBuses, listOfPVBuses, listOfPQBuses)
    # @show newOrder = [listOfSlackBuses listOfPVBuses listOfPQBuses]
    ybusByTypes = ybus[newOrder, newOrder]

    @show typeof(newOrder)
    @show rowNamesByTypes = [string(i) for i in newOrder]

    return ybusByTypes, rowNamesByTypes
end
