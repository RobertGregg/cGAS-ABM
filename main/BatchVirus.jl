if "main" ∈ readdir(pwd()) #Are we in the base folder?
    cd(pwd()*"\\main") #If so go into the main folder
end

using Distributed
addprocs(10)

@everywhere include("Driver.jl")

#How many simulations to run

#=
function VirusRandomPar(prob,i,repeat)
    prob.p.par = reshape.(rand.(sampleDist,nCells),N,N)
    return prob
end

ensembleProb = EnsembleProblem(prob,prob_func = VirusRandomPar)

sim = solve(ensembleProb,CVODE_BDF(linear_solver=:GMRES),callback=cb_infect,
            tstops=tstop,trajectories=10)
=#

simNum = 10
problems = Vector(undef,simNum)

for i=1:simNum
    θnew = deepcopy(θ)
    θnew.par = reshape.(rand.(sampleDist,nCells),N,N)
    problems[i] = ODEProblem(Model!,u0,tspan,θnew)
end

@everywhere VirusSolver(prob) = solve(prob,CVODE_BDF(linear_solver=:GMRES),callback=cb_infect,tstops=tstop)

sim = pmap(VirusSolver,problems)
