#Import Packages needed for simulation
using DifferentialEquations, DiffEqOperators, Sundials #For implementing Differential equations
using LinearAlgebra, SparseArrays, Distributions, Statistics #Linear Algebra and Statistics
using StatsPlots #For graphing

#Include callback functions
include("Callbacks.jl")
include("ModelSetup.jl")

#sol = @time solve(prob,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)
