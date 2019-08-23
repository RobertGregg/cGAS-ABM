#Import Packages needed for simulation
using DifferentialEquations, DiffEqOperators, Sundials #For implementing Differential equations
using LinearAlgebra, SparseArrays, Distributions, Statistics #Linear Algebra and Statistics
using StatsPlots #For graphing

#Include callback functions
include("ModelSetup.jl")
include("Callbacks.jl")

#sol = @time solve(prob,CVODE_BDF(linear_solver=:GMRES),saveat=0.1,callback=cb_infect,tstops=tstop)
#sol = @time solve(prob,CVODE_BDF(linear_solver=:GMRES),callback=cb_infect,tstops=tstop)

#sol = @time solve(prob,ESERK5(),saveat=0.1,callback=cb_all)
#sol = @time solve(prob,CVODE_BDF(linear_solver=:GMRES),saveat=0.1,callback=cb_12,tstops=tstop,progress=true)
