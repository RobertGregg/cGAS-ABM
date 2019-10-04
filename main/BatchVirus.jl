using Distributed

addproc(20)

percent = 0.05 #current value ± percent
sampleDist = @. Uniform((1-percent)*θVals,(1+percent)*θVals)
θ = reshape.(rand.(sampleDist,nCells),N,N)
  
  
ODEProblem(Model!,u0,tspan,θ)