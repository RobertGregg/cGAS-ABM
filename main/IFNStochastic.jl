include("Driver.jl")


#Give a unique parameter set for each cells (randomly choosen), ignore virus parameters?

stepSize = 0.1
percentIFN = range(0.0,1.0,step=0.1)
simRepeat = 10

saveTimePoints = range(tspan[1],tspan[2],step=stepSize)
allStates = [zeros(Int64,length(saveTimePoints),3) for i=1:length(percentIFN),j=1:simRepeat]


for (i,percentIFN) in enumerate(percentIFN)
    for j=1:simRepeat
        @show (percentIFN,j)

        #Reset the parameters
        θ.cellsInfected = fill(Inf,N,N) #Make constant when not testing
        θ.cellsInfected[findall(u0[:,:,2] .> 0.0), 1] .= 0.0
        θ.deathParameter =  ones(N,N)
        θ.cellsDead = fill(Inf,N,N)
        θ.infectFirstAttempt = trues(N,N)
        #Change the parameters
        if percentIFN == 0.0
            θ.par[13] .= 0.0
        else
            #kcat7 produces IFN, make it nonzero ~percentIFN of the time
            kcat7Vals = zeros(N,N)
            θ.par[13][rand(N,N) .< percentIFN] .= 13.006
        end

        probNew = remake(prob; p=θ)
        solNew = solve(probNew,CVODE_BDF(linear_solver=:GMRES),callback=cb)

        for k=1:length(saveTimePoints)
            @inbounds allStates[i,j][k,:] = cellStates(saveTimePoints[k],θ)
        end
    end
end


#=
plot(saveTimePoints,mean(allStates),ribbon=std(allStates),framestyle = :box,linewidth=2,
    labels=[:Healthy,:Infected,:Dead])
xticks!(0:12:48)
xlabel!("Time (hr)")
ylabel!("Number of Cells")
title!("MOI: $(moi)")

savefig("../Figures/CellStates.pdf")


deadMeans = zeros(length(vars))
deadStds = zeros(length(vars))

for i=1:length(vars)
    deadEnd = [allStates[i,j][end,3] for j=1:simRepeat]
    deadMeans[i] = mean(deadEnd)
    deadStds[i] = std(deadEnd)
end
=#



dotArray = [allStates[i,j][end,3] for j=1:simRepeat,i=1:length(percentIFN)]
boxplot(100 .*dotArray./nCells,framestyle=:box,labels=collect(percentIFN),legend=false)
ylabel!("Percentage of Dead Cells")
xlabel!("Percentage of Cells Producing IFNβ")
savefig("../Figures/DeadCells.pdf")

violin(dotArray,framestyle=:box,legend=false)

using DataFrames
CSV.write("BoxData.csv",  DataFrame(dotArray), writeheader=false)

scatter(100 .*percentIFN,100 .*vec(dotArray)./nCells,legend=false,framestyle=:box,xaxis=0:10:100)
ylabel!("Percentage of Dead Cells")
xlabel!("Percentage of Cells Producing IFNb")
