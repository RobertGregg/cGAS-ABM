using DifferentialEquations, Plots

const cellVol = 3e-12 #Cell Volume (liters)
const Na = 6.02e23 #Avagadro's number
const species = 13

#Functtion that converts molecules to nM
m2c(molecule) = @. 1e9*molecule/(cellVol*Na)

const tspan = (0.0,48.0) #Time span for simulation
const statesNames = ["cGAS","DNA","Sting","cGAMP","IRF3","IFNbm","IFNb","STAT",
                     "SOCSm","IRF7m","TREX1m","IRF7","TREX1","Virus"] #for plotting


θvals = [2.6899, 4.8505, 0.0356, 7.487, 517.4056, 22328.3852, 11226.3682,0.9341,
     206.9446, 10305.461, 47639.95,3.8474, 13.006, 78.2048, 0.0209,
     0.0059, 0.001, 0.0112, 0.001, 99.9466, 15.1436,0.0276, 237539.3249,
     61688.259, 0.96, 0.347, 12.8736,1.2399, 1.5101, 0.347, 0.165, 6.9295,
     0.0178]

#Defining the model
function Model(dx,x,p,t)

    #Parameters
    k1f, k1r, k3f, k3r, k4f, kcat5, Km5, k5r, kcat6, Km6, kcat7, Km7, kcat8, Km8, k8f, k9f, k10f1, k10f2, k11f, k12f, k13f, k6f, kcat2, Km2, τ4, τ6, τ7, τ8, τ9, τ10, τ11, τ12, τ13 = p.par
    #Constants from the mass balances
    cGAStot, Stingtot, IRF3tot = p.mass
    #Loop through all grid points with cells

        #Unpack all of the States
            cGAS = x[1]
            DNA = x[2]
            Sting = x[3]
            cGAMP = x[4]
            IRF3 = x[5]
            IFNbm = x[6]
            IFNb = x[7]
            STATP2n = x[8]
            SOCSm = x[9]
            IRF7m = x[10]
            TREX1m = x[11]
            IRF7Pn = x[12]
            TREX1 = x[13]

        # cGAS
            dx[1]=-k1f*cGAS*DNA + k1r*(cGAStot - cGAS)
        # DNA
            dx[2]=-k1f*cGAS*DNA + k1r*(cGAStot - cGAS) - kcat2*TREX1*DNA / (Km2 + DNA)
        # Sting
            dx[3]=-k3f*cGAMP*Sting + k3r*(Stingtot - Sting)
        # cGAMP
            dx[4]=k4f*(cGAStot - cGAS) - k3f*cGAMP*Sting + k3f*(Stingtot - Sting) - τ4*cGAMP
        # IRF3
            dx[5]=-kcat5*IRF3*(Stingtot - Sting) / (Km5 +IRF3) + k5r*(IRF3tot - IRF3)
        # IFNbm
            dx[6]=kcat6*(IRF3tot - IRF3) / (Km6 + (IRF3tot - IRF3)) + k6f*IRF7Pn - τ6*IFNbm
        # IFNb
            dx[7]=kcat7*IFNbm / (Km7 + IFNbm) - τ7*IFNb #Note Cellflux
        # STATP2n
            dx[8]=kcat8*IFNb / (Km8 + IFNb) * 1.0/(1.0+k8f*SOCSm) - τ8*STATP2n
        # SOCSm
            dx[9]=k9f*STATP2n - τ9*SOCSm
        # IRF7m
            dx[10]=k10f1*STATP2n + k10f2*IRF7Pn - τ10*IRF7m
        # TREX1m
            dx[11]=k11f*STATP2n - τ11*TREX1m
        # IRF7Pn
            dx[12]=k12f*IRF7m - τ12*IRF7Pn
        # TREX1
            dx[13]=k13f*TREX1m - τ13*TREX1
end #Model

#These species are not starting at zero
nonZeroSpeciesIdx = [1,3,5] #cGAS, Sting, IRF3
nonZeroSpeciesValues = m2c([1e3, 1e3, 1e4]) #convert to concentration
cGAStot, Stingtot, IRF3tot = nonZeroSpeciesValues #unpack parameters

#Define the initial conditions
u0 = zeros(species)

#Loop through non-zero species and update their concentrations
for (idx,val) in zip(nonZeroSpeciesIdx,nonZeroSpeciesValues)
    u0[idx] =  val
end

#Add some of that dna boi
u0[2] = m2c(1000.0)

#Define a structure to hold all the parameters for the ODE solver
mutable struct ParContainer{T}
  par::T #Rate Constants
  mass::Vector #Mass balances
end

#Create an instance of the structure
θ = ParContainer(θvals,nonZeroSpeciesValues)

const timeBreak = [8.0]

function condition(u,t,integrator)
    t in timeBreak
end

function affect(integrator)
    integrator.p.par[9] = 0.0
end

deadcb = DiscreteCallback(condition,affect)


#Contruct the ODE problem
prob = ODEProblem(Model,u0,tspan,θ)


sol = solve(prob,tstops = timeBreak)

plot!(sol,layout=species,legend=false,title=[t for j=1:1, t in statesNames],size=(800,600))
