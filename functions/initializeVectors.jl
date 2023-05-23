function initializeVectors(busData; MVAb=100)

    N = size(busData, 1)
    PSpecified = zeros(N)
    QSpecified = zeros(N)
    V = zeros(N)
    delta = zeros(N)
    listOfPQBuses = zeros(Int64, N)
    listOfPVBuses = zeros(Int64, N)
    nPQ = 0
    nPV = 0
    n = 0
    nSlack = 0
    listOfNonSlackBuses = zeros(Int64, N)
    listOfSlackBuses = zeros(Int64, N)

    for i = 1:N
        bus = busData.bus[i]
        delta[bus] = 0.0000
        if busData.busType[i] == 0
            nPQ += 1
            listOfPQBuses[nPQ] = bus
            n += 1
            listOfNonSlackBuses[n] = bus
            V[bus] = 1.0000
        elseif busData.busType[i] == 2
            nPV += 1
            listOfPVBuses[nPV] = bus
            n += 1
            listOfNonSlackBuses[n] = bus
            V[bus] = busData.Vset[i]
        elseif busData.busType[i] == 3
            nSlack += 1
            listOfSlackBuses[nSlack] = bus
            V[bus] = busData.Vset[i]
        end
        PSpecified[bus] = busData.PG[i] / MVAb - busData.PL[i] / MVAb
        QSpecified[bus] = busData.QG[i] / MVAb - busData.QL[i] / MVAb
    end

    listOfSlackBuses = reshape(listOfSlackBuses[1:nSlack], nSlack)
    listOfSlackBuses = reshape(listOfSlackBuses[1:nSlack], nSlack)
    listOfPVBuses = reshape(listOfPVBuses[1:nPV], nPV)
    listOfPQBuses = reshape(listOfPQBuses[1:nPQ], nPQ)
    # @show typeof(listOfPQBuses)

    return [PSpecified, QSpecified, V, delta, listOfSlackBuses, listOfPVBuses, listOfPQBuses, listOfNonSlackBuses, nSlack, nPV, nPQ]
end
