using TransformVariables, LogDensityProblems, DynamicHMC,
    DynamicHMC.Diagnostics, Parameters, Statistics, Random

using DifferentialEquations, Distributions,Plots

f1 = function (du,u,p,t)
  du[1] = p[1] * u[1] - p[2] * u[1]*u[2]
  du[2] = -3.0 * u[2] + u[1]*u[2]
end


p = [1.5,1.0]
u0 = [1.0;1.0]
tspan = (0.0,11.0)
prob = ODEProblem(f1,u0,tspan,p)
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
    currentsol = solve(currentprob,Rodas4())

    if sol.retcode != :Success
        return Inf
    else
        # SSE
        err = 0.0
            for (i,t) in enumerate(times)
                err += sum((currentsol(t) .- data[i]).^2)
            end
        return err
    end
end

p = DataODE(tData,solData,prob)
α = (α1=1.5,α2=0.5)

trans = as((α1 = as(Real, 0, 10), α2 = as(Real, 0, 10)))
ℓ = TransformedLogDensity(trans, p)
∇ℓ = ADgradient(:ForwardDiff, ℓ)

results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇ℓ, 10_000)


posterior = transform.(trans, results.chain)
posterior_α = first.(posterior)
mean(posterior_α)


summarize_tree_statistics(results.tree_statistics)

EBFMI(results.tree_statistics)
