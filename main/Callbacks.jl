#The sundials algorithms convert the solution into a one dimensional vector.
#To avoid annoying indexing convert back to 3D shape (xdir,ydir,species)
function SunDials_Reshape(integrator)
    if isa(integrator.alg,CVODE_BDF)
        u =  reshape(integrator.u,N,N,species)
    else
        u = integrator.u
    end
    return u
end
#This function just reverses the process to put back into the solver
function SunDials_Flatten(integrator,u)
    #update the species concentrations (ode solver dependent)
    if isa(integrator.alg,CVODE_BDF)
        u = u[:]
    end
    return u
end



#############################
#Callback 1: healthy → Infect
#############################


function InfectCondition(out,u,t,integrator)
    #println(size(out))
    u = SunDials_Reshape(integrator)

    infectiousState = virusThresold .- u[:,:,14]
    for (idx,virusConc) in enumerate(infectiousState)
        out[idx] = virusConc
    end
    println(out)
    return out
end

function InfectAffect!(integrator,idx)

    #Check if the sundials solver is used, extract current states
    u = SunDials_Reshape(integrator)
    println(idx)
    #Get the index for cells above virusThresold
    #vIdx = findall(u[:,:,14] .> virusThresold)
    vIdx = cellIndicies[idx]

    #calculate prob of getting infected
    Probs = InfectProbability(u,vIdx)

    #Detemine if prob is high enough to infect cell
    for coord in cellIndicies
        #Probs[coord]
        if rand() < Probs[coord] #infect the cell
            #Add Viral DNA to the cell
            u[coord,2] = m2c(1e3)
            #Mark that the cell was infected
            cellsInfected[coord] = integrator.t
        end
    end

    #Update the current solver solution
    integrator.u = SunDials_Flatten(integrator,u)
end

function InfectProbability(u,vIdx)

    #Takes in a viral concentration and outputs a probability for infection
    d(vConc) = ccdf(Normal(m2c(1e3),0.1),vConc)

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

            if cellsInfected[I]==Inf #If it is a healthy cell...
            #Loop through all of the neighbors (and current grid point)
                for J in max(Ifirst, I-Ifirst):min(Ilast, I+Ifirst)
                    if (J in vIdx) & (I != J) #count only infectous neighbors, skip self
                        probNotInfected *= d(virusConc[J])
                    end
                end
            end
            #Probability of getting infected
            probInfected[I] = 1.0 - probNotInfected
    end

    return probInfected
end

cb_infect =VectorContinuousCallback(InfectCondition,InfectAffect!,N^2,rootfind=false)

#############################
#Callback 2: Re-challenge
#############################

function ChellangeCondition(u,t,integrator)
    t ∈ 24.0
end


function ChellangeAffect!(integrator)

    #The sundials algorithm flattens out the solution
    u = SunDials_Reshape(integrator)

    for (i,currentCell) in enumerate(cellIndicies)
        #Are the cells inside the infected region?
        if sqDist([currentCell[1],currentCell[2]],circleOrigin) <= circleRadiusSquared
            u[currentCell,2] = m2c(1e3)
        end
    end

    integrator.u = SunDials_Flatten(integrator,u)
end

cb_challenge = DiscreteCallback(ChellangeCondition,ChellangeAffect!)


#############################
#Callback 3: Infect → dead
#############################

const InfectTimeframe = rand(TruncatedNormal(3.0,0.5,0.0,Inf),N,N)

function DeadCondition(out,u,t,integrator)
    #When does the cell go beyond the alloted tim to be infected?
    out = @. InfectTimeframe - (integrator.t - cellsInfected)
    out = out[:]
end

function DeadAffect!(integrator,idx)
    u = SunDials_Reshape(integrator)
    u[cellIndicies[idx],:] .= 0.0
    integrator.u = SunDials_Flatten(integrator,u)
end

cb_dead = VectorContinuousCallback(DeadCondition,DeadAffect!,N^2)




#Collect all of the callbacks
cb_13 = CallbackSet(cb_infect,cb_dead)
