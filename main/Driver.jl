#Check if we are in the main folder
if "main" ∈ readdir(pwd()) #Are we in the base folder?
    cd(pwd()*"\\main") #If so go into the main folder
end

#Import Packages needed for simulation
using DifferentialEquations, Sundials #For implementing Differential equations
using LinearAlgebra, SparseArrays, Distributions, Statistics, CSV, DataFrames #Linear Algebra and Statistics
using StatsPlots #For graphing

#Include Model and Callback functions
include("ModelSetup.jl")
include("NewVirusCB.jl")
#include("Callbacks.jl")

#What type of infection are we performing?
infectionType = :ISD #ISD, Virus, None

#HEY YOU, did you check the DNA ODE???

if infectionType == :ISD
    sol = @time solve(prob,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)

elseif infectionType == :Virus
    sol = @time solve(prob,CVODE_BDF(linear_solver=:GMRES),callback=cb)
end
<<<<<<< HEAD
=======

#Save the simulation for heatmap
heat10hrVirus = DataFrame()

heat10hrVirus.x = repeat(1:N,N)
heat10hrVirus.y = repeat(1:N,inner=N)
heat10hrVirus.IFN = vec(sol(10.0)[:,:,7])
nfksjd
CSV.write("heat10hrVirus.csv", heat10hrVirus)
>>>>>>> 8a9f9a9969dfc7889aa02c519e6d10ed50132366
