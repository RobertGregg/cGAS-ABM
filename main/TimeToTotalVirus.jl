function cellStates(t)
    #Number of healthy cells at time t
    totaHealthy = sum(θ.cellsInfected .> t)
    #Number of dead cells at time t
    totalDead = sum(θ.cellsDead .< t)
    #Number of infected cells at time t
    totalInfected = nCells - totaHealthy - totalDead
    return [totaHealthy,totalInfected,totalDead]
end

allstates = zeros(Int64,length(sol.t),3)

for i=1:length(sol.t)
    allstates[i,:] = cellStates(sol.t[i])
end

plot(sol.t,allstates,framestyle = :box,linewidth=2,
    labels=[:Healthy,:Infected,:Dead],size=(500,300))
xticks!(0:12:48)
xlabel!("Time (hr)")
ylabel!("Number of Cells")

savefig("../Figures/CellStates.pdf")
