#Import Packages needed for simulation
using DifferentialEquations, Sundials #For implementing Differential equations
using LinearAlgebra, SparseArrays, Distributions, Statistics, CSV, DataFrames #Linear Algebra and Statistics
using StatsPlots #For graphing

#Constants for cell
const cellVol = 3e-12 #Cell Volume (liters)
const Na = 6.02e23 #Avagadro's number
const species = 13 #Number of states within each cell (including virus)

#Functtion that converts molecules to nM
m2c(molecule) = @. 1e9*molecule/(cellVol*Na)

#These species are not starting at zero
nonZeroSpeciesIdx = [1,3,5] #cGAS, Sting, IRF3
nonZeroSpeciesValues = m2c([1e3, 1e3, 1e4]) #convert to concentration

const cGAStot, Stingtot, IRF3tot = nonZeroSpeciesValues

#Paramter values for the ODEs
θNames = [:k1f, :k1r, :k3f, :k3r, :k4f, :kcat5, :Km5, :k5r, :kcat6, :Km6, :kcat7,
:Km7, :kcat8, :Km8, :k8f, :k9f, :k10f1, :k10f2, :k11f, :k12f, :k13f, :k6f, :kcat2,
:Km2, :τ4, :τ6, :τ7, :τ8, :τ9, :τ10, :τ11, :τ12, :τ13, :k14f,:τ14]
θ = [2.6899, 4.8505, 0.0356, 7.487, 517.4056, 22328.3852, 11226.3682,0.9341,
         206.9446, 10305.461, 47639.70295,3.8474, 13.006, 78.2048, 0.0209,
         0.0059, 0.001, 0.0112, 0.001, 99.9466, 15.1436,0.0276, 237539.3249,
         61688.259, 0.96, 0.347, 12.2428736,1.2399, 1.5101, 0.347, 0.165, 6.9295,
         0.0178]


# Define the discretized PDE as an ODE function
function Model!(du,u,p,t)
  #Species
  cGAS, DNA, Sting, cGAMP, IRF3, IFNβm, IFNβ, STAT, SOCSm, IRF7m, TREX1m, IRF7, TREX1 = u

  #Parameters
  k1f, k1r, k3f, k3r, k4f, kcat5, Km5, k5r, kcat6, Km6, kcat7, Km7, kcat8, Km8, k8f, k9f, k10f1, k10f2, k11f, k12f, k13f, k6f, kcat2, Km2, τ4, τ6, τ7, τ8, τ9, τ10, τ11, τ12, τ13 = p

  #Update derivatives for each species according to model
  du[1] = -k1f*cGAS*DNA + k1r*(cGAStot - cGAS)
  du[2] = -k1f*cGAS*DNA + k1r*(cGAStot - cGAS) - kcat2*TREX1*DNA / (Km2 + DNA)
  du[3] = -k3f*cGAMP*Sting + k3r*(Stingtot - Sting)
  du[4] = k4f*(cGAStot - cGAS) - k3f*cGAMP*Sting + k3f*(Stingtot - Sting) - τ4*cGAMP
  du[5] = -kcat5*IRF3*(Stingtot - Sting) / (Km5 +IRF3) + k5r*(IRF3tot - IRF3)
  du[6] = kcat6*(IRF3tot - IRF3) / (Km6 + (IRF3tot - IRF3)) + k6f*IRF7 - τ6*IFNβm
  du[7] = kcat7*IFNβm / (Km7 + IFNβm) - τ7*IFNβ
  du[8] = kcat8*IFNβ / (Km8 + IFNβ) * 1.0/(1.0+k8f*SOCSm) - τ8*STAT
  du[9] = k9f*STAT - τ9*SOCSm
  du[10] = k10f1*STAT + k10f2*IRF7 - τ10*IRF7m
  du[11] = k11f*STAT - τ11*TREX1m
  du[12] = k12f*IRF7m - τ12*IRF7
  du[13] = k13f*TREX1m - τ13*TREX1
end


#Define the initial conditions
u0 = zeros(species)
u0[nonZeroSpeciesIdx] .= nonZeroSpeciesValues
#Add some DNA into the mix
u0[2] = m2c(1000)

#Time span for simulation
const tspan = (0.0,48.0)

#Contruct the ODE problem
prob = ODEProblem(Model!,u0,tspan,θ)



function KD(parChange)

    KnockDownVals=1:-0.25:0
    percentLabels = [string(i)*"%" for i=0:25:100]
    kdSamples = length(KnockDownVals)

    solKD =  Vector(undef,kdSamples)

    for (i,percent) in enumerate(KnockDownVals)
        θcurrent = deepcopy(θ)
        θcurrent[parChange] = θcurrent[parChange] * percent

        println(θcurrent[[20,23]])

        probKD = remake(prob; p=θcurrent)
        sol = solve(probKD,saveat=0.1)
        solKD[i] = vec(sol[7,:])
    end

    #return p
    return solKD
end


KDSols = map(KD,[20,23])


#Save the simulation for KD
KDDataODE = DataFrame()
KDDataODE.IRF7KD = vcat(KDSols[1]...)
KDDataODE.TREXKD = vcat(KDSols[2]...)

timeLength = 481
KnockDownVals=1:-0.25:0
percentLabels = [string(i)*"%" for i=0:25:100]
kdSamples = length(KnockDownVals)

KDDataODE.Time = repeat(0:0.1:48,outer=kdSamples)
KDDataODE.Percent = repeat(KnockDownVals,inner=timeLength)

CSV.write("KDDataODE.csv",KDDataODE)
