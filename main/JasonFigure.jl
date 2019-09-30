#Import Packages needed for simulation
using DifferentialEquations, Sundials #For implementing Differential equations
using LinearAlgebra, SparseArrays, Distributions, Statistics, CSV #Linear Algebra and Statistics
using StatsPlots #For graphing

#Include callback functions
include("ModelSetup.jl")
include("Callbacks.jl")

#sol = @time solve(prob,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)
sol = @time solve(prob,CVODE_BDF(linear_solver=:GMRES),callback=cb_infect,tstops=tstop)



Maxval = maximum(sol[:,:,7,:])
i=25.0
heatmap(sol(i)[:,:,7],clims=(0.0,Maxval))


rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
plot!(rectangle(20,20,45,45),legend=false,color = false,linecolor =:green,axis=nothing,
    linewidth = 3)
savefig("grantfig.pdf")

heatmap(sol(i)[45:65,45:65,7],clims=(0.0,Maxval))
savefig("grantfig1.pdf")


A = [i for i=1:100, j=1:100]

heatmap(A, c=ColorGradient([:blue, :lightgreen,:yellow]))
