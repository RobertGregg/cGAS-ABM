#Importing neccessary Packages
using  DifferentialEquations, ModelingToolkit, RecursiveArrayTools, DiffEqBayes
using Turing, Distributions, StatsPlots

#Using ModelingToolkit to define the ODE equations
@parameters t α β
@variables x(t) y(t)
@derivatives D'~t

#Define the Differential Equations
eqs = [
 D(x) ~ α*x - x*y
 D(y) ~ -β*y + x*y
]

de = ODESystem(eqs)
f = ODEFunction(de,[x,y],[α,β])

p = [1.5,3.0]
u0 = [1.0,1.0]
tspan = (0.0,10.0)
prob = ODEProblem(f,u0,tspan,p)

σ = 0.3                     # noise, fixed for now
t = collect(range(1,10,length=20))   # observation times
sol = solve(prob,Tsit5())

randomized = VectorOfArray([(sol(t[i]) + σ * randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)

scatter(t,data')
plot!(sol)

priors = [Uniform(0.0,10.0),Uniform(0.0,10.0)]
bayesian_result_turing = turing_inference(prob,Tsit5(),t,data,priors;num_samples=1_000)
