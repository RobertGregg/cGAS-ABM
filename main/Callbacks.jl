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

#Trading readbility for speed here
function InfectCondition(u,t,integrator)
    #println(integrator.t)

    t ∈ tstop

    # SunDials_Reshape(integrator)
    #cutoff = m2c(1000.)
    #@inbounds for i = 1:size(u)[1], j=1:size(u)[2]
    #    if u[i,j,end] > cutoff
    #        return true
    #    end
    #end
    #return false
end

function InfectAffect(integrator)

    #Check if the sundials solver is used, extract current states
    u = SunDials_Reshape(integrator)

    #Get the index for the current time
    #tIdx = searchsortedfirst(tstop,integrator.t)

    #calculate prob of getting infected
    Probs = InfectProbability(u)

    #Cells from the previous time point are still infected
    #cellsInfected[:,:,tIdx+1] .= cellsInfected[:,:,tIdx] #faster in loop?

    #Detemine if prob is high enough to infect cell
    @inbounds for coord in cellIndicies
        #Probs[coord]
        if rand() < Probs[coord] #infect the cell
            #Add Viral DNA to the cell
            u[coord,2] = rand(Uniform(0.0,m2c(1e3)))
            #Mark that the cell was infected
            cellsInfected[coord] = integrator.t
        end
    end

    #Update the current solver solution
    integrator.u = SunDials_Flatten(integrator,u)


end

function InfectProbability(u)

    #Takes in a viral concentration and outputs a probability for infection
    d(vConc) = ccdf(Normal(0.5*m2c(1e3),0.1),vConc)

    #Initialize an array of probabilities
    probInfected = zeros(N,N)

    #Retrieve the indices to iterate over
    Index = CartesianIndices(probInfected)
    #Get the first and last index
    Ifirst, Ilast = first(Index), last(Index)

    #Loop through all of the grid points
    @inbounds for I in Index
            probNotInfected = 1.0

            if cellsInfected[I]==Inf #If it is a healthy cell...
            #Loop through all of the neighbors (and current grid point)
                for J in max(Ifirst, I-Ifirst):min(Ilast, I+Ifirst)
                    if !isinf(cellsInfected[J]) & (I != J) #count only infected neighbors, skip self
                        probNotInfected *= d(u[J,14])
                    end
                end
            end
            #Probability of getting infected
            probInfected[I] = 1.0 - probNotInfected
    end

    return probInfected
end

cb_infect = DiscreteCallback(InfectCondition,InfectAffect)

#############################
#Callback 2: Re-challenge
#############################

function ChellangeCondition(u,t,integrator)
    t ∈ 24.0
end


function ChellangeAffect(integrator)

    #The sundials algorithm flattens out the solution
    u = SunDials_Reshape(integrator)

    for (i,currentCell) in enumerate(cellIndicies)
        #Are the cells inside the infected region?
        if sqDist([currentCell[1],currentCell[2]],[50,50]) <= circleRadiusSquared
            u[currentCell,2] = m2c(1e3)
        end
    end

    integrator.u = SunDials_Flatten(integrator,u)
end

cb_challenge = DiscreteCallback(ChellangeCondition,ChellangeAffect)


#############################
#Callback 3: Infect → dead
#############################

const InfectTimeframe = rand(TruncatedNormal(3.0,0.5,0.0,Inf),N,N)

function DeadCondition(u,t,integrator)
    #When does the cell go beyond the alloted tim to be infected?
    any( InfectTimeframe .< (t .- cellsInfected) )
end

function DeadAffect!(integrator)
    u = SunDials_Reshape(integrator)

    deadIdx = findall(InfectTimeframe .< (integrator.t .- cellsInfected))

    for idx in deadIdx
        u[cellIndicies[idx],:] .= 0.0
    end

    integrator.u = SunDials_Flatten(integrator,u)
end

cb_dead = DiscreteCallback(DeadCondition,DeadAffect!)



#Collect all of the callbacks
cb_12 = CallbackSet(cb_infect,cb_challenge)
cb_13 = CallbackSet(cb_infect,cb_dead)
