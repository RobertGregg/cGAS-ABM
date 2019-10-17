#Check if we are in the main folder
if "main" ∈ readdir(pwd()) #Are we in the base folder?
    cd(pwd()*"\\main") #If so go into the main folder
end

#Import Packages needed for simulation
using DifferentialEquations, Sundials #For implementing Differential equations
using LinearAlgebra, SparseArrays, Distributions, Statistics, CSV #Linear Algebra and Statistics
using StatsPlots #For graphing

#Include Model and Callback functions
include("ModelSetup.jl")
include("NewVirusCB.jl")
#include("Callbacks.jl")

#What type of infection are we performing?
infectionType = :Virus #ISD, Virus, None

if infectionType == :ISD
    sol = @time solve(prob,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)

elseif infectionType == :Virus
    sol = @time solve(prob,CVODE_BDF(linear_solver=:GMRES),callback=cb)
end


#TODO Where should this go?

function cellStates(t)
    #Number of healthy cells at time t
    totaHealthy = sum(θ.cellsInfected .>= t)
    #Number of dead cells at time t
    totalDead = sum(θ.cellsDead .<= t)
    #Number of infected cells at time t
    totalInfected = nCells - totaHealthy - totalDead
    return [totaHealthy,totalInfected,totalDead]
end
