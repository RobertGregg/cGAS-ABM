if "main" ∈ readdir(pwd()) #Are we in the base folder?
    cd(pwd()*"\\main") #If so go into the main folder
end

using Distributed
addprocs(10)

@everywhere include("Driver.jl")
@everywhere using JLD2, FileIO, UnicodePlots


#=
function VirusRandomPar(prob,i,repeat)
    prob.p.par = reshape.(rand.(sampleDist,nCells),N,N)
    return prob
end

ensembleProb = EnsembleProblem(prob,prob_func = VirusRandomPar)

sim = solve(ensembleProb,CVODE_BDF(linear_solver=:GMRES),callback=cb_infect,
            tstops=tstop,trajectories=10)
=#

#How many simulations to run
simNum = 5
problems = Vector(undef,simNum)

for i=1:simNum
    θnew = deepcopy(θ)
    θnew.par = reshape.(rand.(sampleDist,nCells),N,N)
    problems[i] = ODEProblem(Model!,u0,tspan,θnew)
end

@everywhere VirusSolver(prob) = solve(prob,CVODE_BDF(linear_solver=:GMRES),callback=cb_infect,tstops=tstop)

sim = pmap(VirusSolver,problems)


timePoints = [sim[i].t for i=1:simNum]
infectPercent = [zeros(length(sim[i].t)) for i=1:simNum]

for (i,currentSol) in enumerate(sim)
	for t = 1:length(currentSol.t)
		infectPercent[i][t] = 1.0 - sum(currentSol[:,:,2,t] .== 0.0)/N^2
	end
end

@save "infectPercent.jld2" infectPercent, timePoints

plt = plot(timePoints,infectPercent,legend=false)
savefig("../Figures/InfectionVary.pdf")


