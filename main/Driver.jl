#Check if we are in the main folder
if "main" ∈ readdir(pwd()) #Are we in the base folder?
    cd(pwd()*"\\main") #If so go into the main folder
end

#Import Packages needed for simulation
using DifferentialEquations, Sundials #For implementing Differential equations
using LinearAlgebra, SparseArrays, Distributions, Statistics, CSV, DataFrames #Linear Algebra and Statistics
using StatsPlots #For graphing
using JLD2, FileIO

#Include Model and Callback functions
include("ModelSetup.jl")
include("NewVirusCB.jl")
#include("Callbacks.jl")

#What type of infection are we performing?
infectionType = :Virus #ISD, Virus, None

#HEY YOU, did you check the DNA ODE???

if infectionType == :ISD
    sol = @time solve(prob,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)

elseif infectionType == :Virus
    sol = @time solve(prob,CVODE_BDF(linear_solver=:GMRES),callback=cb)
end


saveIFN = vec(sol[:,:,7,:])


@save "saveIFN.jld2" saveIFN
