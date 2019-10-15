include("Driver.jl")

moiRange = exp10.(range(-3.0,2.0,length=10))
repetitions = 10
vars = range(0.0,0.5,step=0.1)

cellInfectedTot = [zeros(length(moiRange),repetitions) for i=1:length(vars)]


for (i,newVar) in enumerate(vars)
    for (j,newMOI) in enumerate(moiRange)
        for k=1:repetitions
            #Change the initial condition
            probDistInfected = Poisson(newMOI)
            u0[:,:,2] = @. m2c(1e3*rand(probDistInfected,N,N))

            #reset the matrix keeping track of the infected cells
            cellsInfected = fill(Inf,N,N) #Make constant when not testing
            cellsInfected[findall(u0[:,:,2] .> 0.0), 1] .= 0.0

            #Change the parameters
            if newVar == 0.0
                θ.par = fill.(θVals,N,N)
            else
                sampleDist = @. Uniform((1-newVar)*θVals,(1+newVar)*θVals)
                θ.par = reshape.(rand.(sampleDist,nCells),N,N)
            end

            probNew = remake(prob; u0=u0, p=θ)
            solNew = @time solve(probNew,CVODE_BDF(linear_solver=:GMRES),saveat=0.1,callback=cbVec)

            cellInfectedTot[i][j,k] = 1.0 - sum(solNew(tspan[2])[:,:,2] .== 0.0)/N^2
        end
    end
end


cellInfectedTot = newTot
plot(log10.(moiRange),100 .*mean.(cellInfectedTot,dims=2),
            marker=:auto,frame=:box,label=string.(collect(vars)),legendtitle="Variability")
xlabel!("log10(MOI)")
ylabel!("Total Cells Infected (%)")

savefig("../Figures/MOI.pdf")
