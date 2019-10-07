if "main" ∈ readdir(pwd()) #Are we in the base folder?
    cd(pwd()*"\\main") #If so go into the main folder
end

using Distributed
addprocs(10)

@everywhere include("Driver.jl")
#@everywhere using JLD2, FileIO, UnicodePlots


#=
#Not working because of callback?
function VirusRandomPar(prob,i,repeat)
    prob.p.par = reshape.(rand.(sampleDist,nCells),N,N)
    return prob
end

ensembleProb = EnsembleProblem(prob,prob_func = VirusRandomPar)

sim = solve(ensembleProb,CVODE_BDF(linear_solver=:GMRES),callback=cb_infect,
            tstops=tstop,trajectories=10)
=#

#How many simulations to run
const simNum = 5

@everywhere VirusSolver(prob) = solve(prob,CVODE_BDF(linear_solver=:GMRES),callback=cb_infect,tstops=tstop)

@everywhere function VirusBatch(σ,θ)
	problems = Vector(undef,simNum)

	sampleDist = @. Uniform((1-σ)*θVals,(1+σ)*θVals)
	for i=1:simNum
	    θnew = deepcopy(θ)
	    θnew.par = reshape.(rand.(sampleDist,nCells),N,N)
	    problems[i] = ODEProblem(Model!,u0,tspan,θnew)
	end

	sim = pmap(VirusSolver,problems)

	timePoints = [sim[i].t for i=1:simNum]
	infectPercent = [zeros(length(sim[i].t)) for i=1:simNum]

	for (i,currentSol) in enumerate(sim)
		for t = 1:length(currentSol.t)
			infectPercent[i][t] = 1.0 - sum(currentSol[:,:,2,t] .== 0.0)/N^2
		end
	end

	return (timePoints,infectPercent)
end

percents = 0.05:0.2
simPercent = Vector(undef,length(percents))
plt = plot()

for (i,percent) in enumerate(percents)
	sim[i] = VirusBatch(percent,θ)
	plot!(sim[i][1],mean(sim[i][2]),ribbon=std(sim[i][2]),label="$(percent)")
end


#@save "infectPercent.jld2" sim

plt = plot(timePoints,infectPercent,legend=false)
savefig("../Figures/InfectionVary.pdf")
