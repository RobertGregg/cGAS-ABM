#Some options to choose in the setup
infectionMethod = :wash #wash or drop
parameterVary = true #All cells have different parameters?

#Constants for all cell
const N=5 #number of grid points along one dimensions
const nCells = N^2 #number of cells in the simulation
const cellVol = 3e-12 #Cell Volume (liters)
const Na = 6.02e23 #Avagadro's number
const D=1.0 #Diffusion coefficient (μm^2/s)
const species = 14 #Number of states within each cell (including virus)
const moi = 1.0 #Multicity of infection

#Converts molecules to nM
m2c(molecule) = @. 1e9*molecule/(cellVol*Na)

#Paramter values for the ODEs
θVals = [2.6899, 4.8505, 0.0356, 7.487, 517.4056, 22328.3852, 11226.3682,0.9341,
         206.9446, 10305.461, 47639702.95,3.8474, 50.006, 78.2048, 0.0209,
         0.0059, 0.001, 0.0112, 0.001, 99.9466, 15.1436,0.0276, 237539.3249,
         61688.259, 0.96, 1.347, 12242.8736,1.2399, 1.5101, 0.347, 0.165, 6.9295,
         0.0178]
θVirus = [2.0, 1.0] # k14f τ14 (Virus Parameters)


if parameterVary
  #Append the virus parameters to the orginal parameters
  append!(θVals,θVirus)
  #Give a unique parameter set for each cells (randomly choosen)
  percent = 0.05
  sampleDist = @. Uniform((1-percent)*θVals,(1+percent)*θVals)
  θ = reshape.(rand.(sampleDist,nCells),N,N)
else
  θ = append!(θVals,θVirus)
end

const virusThresold = m2c(1e3)

const tspan = (0.0,48.0) #Time span for simulation
const tstop = 1:tspan[2] #Times where the simulation stops and check for virus movement
const statesNames = ["cGAS","DNA","Sting","cGAMP","IRF3","IFNbm","IFNb","STAT",
                     "SOCSm","IRF7m","TREX1m","IRF7","TREX1","Virus"] #for plotting


 function Laplacian1D(n)
     #Derivative order
         DerOrder = 2
     #Approximation order
         ApproxOrder = 2
     #Grid spacing (diameter of cell in μm)
         h = 1.0
     #Second order approximation of the second derivative
         L = DerivativeOperator{Float64}(DerOrder,ApproxOrder,h,n,:Neumann0,:Neumann0)
     return sparse(L)
 end

const Lx = Laplacian1D(N) #x direction of the 1-D Laplacian
const Ly = Lx' #y direction of the 1-D Laplacian

#Create containers to hold the matrix multiplications
Ly●IFNβ = zeros(N,N)
IFNβ●Lx = zeros(N,N)

# Define the discretized PDE as an ODE function
function Model!(du,u,p,t)

  #Species
  cGAS = @view u[:,:,1]
  DNA = @view u[:,:,2]
  Sting = @view u[:,:,3]
  cGAMP = @view u[:,:,4]
  IRF3 = @view u[:,:,5]
  IFNβm = @view u[:,:,6]
  IFNβ = @view u[:,:,7]
  STAT = @view u[:,:,8]
  SOCSm = @view u[:,:,9]
  IRF7m = @view u[:,:,10]
  TREX1m = @view u[:,:,11]
  IRF7 = @view u[:,:,12]
  TREX1 = @view u[:,:,13]
  Virus = @view u[:,:,14]

  #Derivatives
  d_cGAS = @view du[:,:,1]
  d_DNA = @view du[:,:,2]
  d_Sting = @view du[:,:,3]
  d_cGAMP = @view du[:,:,4]
  d_IRF3 = @view du[:,:,5]
  d_IFNβm = @view du[:,:,6]
  d_IFNβ = @view du[:,:,7]
  d_STAT = @view du[:,:,8]
  d_SOCSm = @view du[:,:,9]
  d_IRF7m = @view du[:,:,10]
  d_TREX1m = @view du[:,:,11]
  d_IRF7 = @view du[:,:,12]
  d_TREX1 = @view du[:,:,13]
  d_Virus = @view du[:,:,14]

  #Parameters
  k1f, k1r, k3f, k3r, k4f, kcat5, Km5, k5r, kcat6, Km6, kcat7, Km7, kcat8, Km8, k8f, k9f, k10f1, k10f2, k11f, k12f, k13f, k6f, kcat2, Km2, τ4, τ6, τ7, τ8, τ9, τ10, τ11, τ12, τ13, k14f, τ14  = p

  #Calculate diffusion of interferon (Ly*IFN + IFN*Lx)
  mul!(Ly●IFNβ, Ly, IFNβ)
  mul!(IFNβ●Lx,IFNβ,Lx)

  #Sum the x and y componentsd
   L●IFNβ = @. D*(Ly●IFNβ + IFNβ●Lx)

  #Update derivatives for each species according to model
  @. d_cGAS = -k1f*cGAS*DNA + k1r*(cGAStot - cGAS)
  @. d_DNA = -k1f*cGAS*DNA + k1r*(cGAStot - cGAS) - kcat2*TREX1*DNA / (Km2 + DNA)
  @. d_Sting = -k3f*cGAMP*Sting + k3r*(Stingtot - Sting)
  @. d_cGAMP = k4f*(cGAStot - cGAS) - k3f*cGAMP*Sting + k3f*(Stingtot - Sting) - τ4*cGAMP
  @. d_IRF3 = -kcat5*IRF3*(Stingtot - Sting) / (Km5 +IRF3) + k5r*(IRF3tot - IRF3)
  @. d_IFNβm = kcat6*(IRF3tot - IRF3) / (Km6 + (IRF3tot - IRF3)) + k6f*IRF7 - τ6*IFNβm
  @. d_IFNβ = kcat7*IFNβm / (Km7 + IFNβm) - τ7*IFNβ + L●IFNβ #Add the diffusion in here
  @. d_STAT = kcat8*IFNβ / (Km8 + IFNβ) * 1.0/(1.0+k8f*SOCSm) - τ8*STAT
  @. d_SOCSm = k9f*STAT - τ9*SOCSm
  @. d_IRF7m = k10f1*STAT + k10f2*IRF7 - τ10*IRF7m
  @. d_TREX1m = k11f*STAT - τ11*TREX1m
  @. d_IRF7 = k12f*IRF7m - τ12*IRF7
  @. d_TREX1 = k13f*TREX1m - τ13*TREX1
  @. d_Virus = k14f*DNA - τ14*Virus
  #@. d_Virus = DNA*(1.0-DNA/k14f)
end

#Loop though all of the species and determine initial conditions

u0 = zeros(N,N,species) #Initial Condition

#These species are not starting at zero
nonZeroSpeciesIdx = [1,3,5] #cGAS, Sting, IRF3
nonZeroSpeciesValues = m2c([1e3, 1e3, 1e4])
const cGAStot, Stingtot, IRF3tot = nonZeroSpeciesValues

#Loop through non-zero species and update their concentrations
for (idx,val) in zip(nonZeroSpeciesIdx,nonZeroSpeciesValues)
    u0[:,:,idx] .=  val
end

#Define a set of indices for looping through every cell
const cellIndicies = CartesianIndices(u0[:,:,1])


if infectionMethod == :wash

  probDistInfected = Poisson(moi)
  u0[:,:,2] = @. m2c(1e3*rand(probDistInfected,N,N))

elseif infectionMethod == :drop

  circleOrigin = [10,10] #Where is the drop
  circleRadiusSquared = 15^2 #How big is the drop
  #Calculate squared distances
  sqDist(x,c) = reduce(+, @. (x-c)^2)
  #Loop though cells and check if they are infected
  for (i,currentCell) in enumerate(cellIndicies)
      #Are the cells inside the infected region?
      if sqDist([currentCell[1],currentCell[2]],circleOrigin) <= circleRadiusSquared
          u0[currentCell,2] = m2c(1e3)
      end
  end

end

#Keep track of infected cells (save time when infected, Inf means not infected)
const cellsInfected = fill(Inf,N,N) #Make constant when not testing
cellsInfected[findall(u0[:,:,2] .> 0.0), 1] .= 0.0


#Contruct the ODEs
prob = ODEProblem(Model!,u0,tspan,θ)
#Solve the problem
#sol = @time solve(prob,ROCK4(),saveat=0.1)
