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
    sol = @time solve(prob,CVODE_BDF(linear_solver=:GMRES),saveat=0.1,callback=cb)
end


#TODO Remove later, just to test larger simulation on server

#Get the max value for IFNβ concentration for whole simulation
Maxval = maximum(sol[:,:,7,:])
#Loop through time points to make an animation
anim = @animate for i = tspan[1]:0.1:tspan[2]
    heatmap(sol(i)[:,:,7],clims=(0.0,Maxval),
            title="Time = " * string(i) * " hrs")
end

gif(anim,"../Figures/InterferonAnimation.gif")

function cellStates(t)
    #Number of healthy cells at time t
    totaHealthy = sum(θ.cellsInfected .> t)
    #Number of dead cells at time t
    totalDead = sum(θ.cellsDead .< t)
    #Number of infected cells at time t
    totalInfected = nCells - totaHealthy - totalDead
    return [totaHealthy,totalInfected,totalDead]
end

allstates = zeros(Int64,length(sol.t),3)

for i=1:length(sol.t)
    allstates[i,:] = cellStates(sol.t[i])
end

plot(sol.t,allstates,framestyle = :box,linewidth=2,
    labels=[:Healthy,:Infected,:Dead],size=(500,300))
xticks!(0:12:48)
xlabel!("Time (hr)")
ylabel!("Number of Cells")

savefig("../Figures/CellStates.pdf")
