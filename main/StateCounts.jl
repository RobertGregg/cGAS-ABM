using ProgressMeter
include("Driver.jl")


function cellStates(t,θ)
    #Number of healthy cells at time t
    totaHealthy = sum(θ.cellsInfected .> t)
    #Number of dead cells at time t
    totalDead = sum(θ.cellsDead .< t)
    #Number of infected cells at time t
    totalInfected = nCells - totaHealthy - totalDead
    return [totaHealthy,totalInfected,totalDead]
end

#Give a unique parameter set for each cells (randomly choosen), ignore virus parameters?

stepSize = 0.1
vars = range(0.0,0.5,step=stepSize)
simRepeat = 30

saveTimePoints = range(tspan[1],tspan[2],step=stepSize)
allStates = [zeros(Int64,length(saveTimePoints),3) for i=1:length(vars),j=1:simRepeat]


for (i,newVar) in enumerate(vars)
    for j=1:simRepeat
        @show (newVar,j)
        #Change the initial condition
        #probDistInfected = Poisson(moi)
        #u0[:,:,2] = @. m2c(1e3*rand(probDistInfected,N,N))

        #Reset the parameters
        θ.cellsInfected = fill(Inf,N,N) #Make constant when not testing
        θ.cellsInfected[findall(u0[:,:,2] .> 0.0), 1] .= 0.0
        θ.deathParameter =  ones(N,N)
        θ.cellsDead = fill(Inf,N,N)
        θ.infectFirstAttempt = trues(N,N)
        #Change the parameters
        if newVar == 0.0
            θ.par = fill.(θVals,N,N)
        else
            sampleDist = @. Uniform((1.0-newVar)*θVals,(1.0+newVar)*θVals)
            θ.par = reshape.(rand.(sampleDist,nCells),N,N)
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
=#

deadMeans = zeros(length(vars))
deadStds = zeros(length(vars))

for i=1:length(vars)
    deadEnd = [allStates[i,j][end,3] for j=1:simRepeat]
    deadMeans[i] = mean(deadEnd)
    deadStds[i] = std(deadEnd)
end




dotArray = [allStates[i,j][end,3] for j=1:simRepeat,i=1:length(vars)]
baxplot(dotArray,framestyle=:box,labels=collect(vars))
ylabel!("Number of Dead Cells")
xlabel!("Parameter Standard Deviation")
savefig("../Figures/DeadCells.pdf")

#violin(dotArray,framestyle=:box,legend=false)

using DataFrames
CSV.write("BoxData.csv",  DataFrame(dotArray), writeheader=false)
