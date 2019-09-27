using TransformVariables, LogDensityProblems, DynamicHMC,
    DynamicHMC.Diagnostics, Parameters, Statistics, Random, ModelingToolkit

using DifferentialEquations, DiffEqParamEstim, Distributions,Plots


## Define the ODE model
@parameters t a b c
@variables u1(t) u2(t)
@derivatives D'~t

eqs = [
D(u1) ~ a*u1 - b*u1*u2
D(u2) ~ -c*u2 + u1*u2
]

de = ODESystem(eqs)
f1 = ODEFunction(de,[u1,u2],[a,b,c])




p = [1.5,1.0,3.0]
u0 = ones(length(eqs))
tspan = (0.0,10.0)
prob = ODEProblem(f1,u0,tspan,p)
alg = Tsit5()
sol = solve(prob)

tData = range(0.0,10.0,length=20)
solData = [sol(t) .+ rand(Normal(0.0,0.1),2) for t in tData]

#=
plot(sol)
scatter!(tData,[s[1] for s in solData])
scatter!(tData,[s[2] for s in solData])
=#

struct DataODE{T,D,P}
  times::T #Time points for data
  data::D # LFC of data
  prob::P
end


function (problem::DataODE)(θ)
    @unpack times,data,prob = problem       # extract the data

    #Solve the ODE
    currentprob = remake(prob;u0=convert.(eltype(θ),prob.u0),p=θ)
    currentsol = solve(currentprob)

    if sol.retcode != :Success
        return Inf
    else
        # SSE
        return mapreduce(x -> sum(x.^2), +, currentsol(times).u .- data)
    end
end




lossFunc = DataODE(tData,solData,prob)

SumSqErr = build_loss_objective(prob,alg,L2Loss(tData,hcat(solData...)))






trans = as((a = as(Real, 0, 10), b = as(Real, 0, 10),c = as(Real, 0, 10)))
ℓ = TransformedLogDensity(trans, SumSqErr)
∇ℓ = ADgradient(:ForwardDiff, ℓ)

results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇ℓ, 10_000)


posterior = transform.(trans, results.chain)
posterior_α = first.(posterior)
mean(posterior_α)


summarize_tree_statistics(results.tree_statistics)

EBFMI(results.tree_statistics)
