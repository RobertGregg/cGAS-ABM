#Import Packages needed for simulation
using DifferentialEquations, DiffEqOperators #For implementing Differential equations
using LinearAlgebra, SparseArrays, Distributions, Statistics #Linear Algebra and Statistics
using StatsPlots #For graphing

#Include callback functions
include("Callbacks.jl")
include("ModelSetup.jl")

sol = @time solve(prob,ESERK5(),saveat=0.1,callback=cb_infect,tstops=tstop)
