#Callback 1: healthy → Infect

function InfectCondition(u,t,integrator)
    t ∈ tstop
end

function InfectAffect(integrator)

    #The sundials algorithm flattens out the solution
    if isa(integrator.alg,CVODE_BDF)
        u =  reshape(integrator.u,N,N,species)
    else
        u = integrator.u
    end

    #Get the index for the current time
    tIdx = searchsortedfirst(tstop,integrator.t)

    #calculate prob of getting infected
    Probs = InfectProbability(u,tIdx)

    #Cells from the previous time point are still infected
    cellsInfected[:,:,tIdx+1] .= cellsInfected[:,:,tIdx] #faster in loop?

    #Detemine if prob is high enough to infect cell
    for coord in CartesianIndices(Probs)
        #Probs[coord]
        if rand() < Probs[coord] #infect the cell
            #Add Viral DNA to the cell
            u[coord,2] = m2c(1e3)
            #Mark that the cell was infected
            cellsInfected[coord,tIdx+1] = 1
        end
    end

    #update the species concentrations (ode solver dependent)
    if isa(integrator.alg,CVODE_BDF)
        integrator.u = u[:]
    else
        integrator.u = u
    end


end

function InfectProbability(u,tIdx)

    #Takes in a viral concentration and outputs a probability for infection
    d(vConc) = ccdf(Normal(0.8,0.1),vConc)

    virusConc = u[:,:,end] #current viral concentrations in each cell

    #Initialize an array of probabilities
    probInfected = zeros(N,N)

    #Retrieve the indices to iterate over
    Index = CartesianIndices(probInfected)
    #Get the first and last index
    Ifirst, Ilast = first(Index), last(Index)

    #Loop through all of the grid points
    for I in Index
            probNotInfected = 1.0

            if cellsInfected[I,tIdx]==0 #If it is a healthy cell...
            #Loop through all of the neighbors (and current grid point)
                for J in max(Ifirst, I-Ifirst):min(Ilast, I+Ifirst)
                    if cellsInfected[J,tIdx]==1 & (I != J) #count only infected neighbors, skip self
                        probNotInfected *= d(virusConc[J])
                    end
                end
            end
            #Probability of getting infected
            probInfected[I] = 1.0 - probNotInfected
    end

    return probInfected
end

cb_infect = DiscreteCallback(InfectCondition,InfectAffect)
