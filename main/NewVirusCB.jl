#Continuous callback

#Create an array that keeps track of whether or not a cell has tried to infect neighbors
infectFirstAttempt = trues(N,N)

#############################
#Callback 1: healthy → Infect
#############################
function VecCondition(u,t,integrator)
    #Reshape the array into 3D
    uReshaped =  reshape(u,N,N,species)

    #Check the whole virus grid for cells with viral conc. greater than 0.4 nM
    return any( @. (uReshaped[:,:,end] > 0.4) & infectFirstAttempt)
end


#These don't need to be constructed every loop
    const Index = CartesianIndices(zeros(N,N))
    #Get the first and last index
    const Ifirst, Ilast = first(Index), last(Index)

function VecAffect(integrator)
    #Probability that an infected cell will spread infection
    chanceOfInfection = Bernoulli(0.5)
    @show integrator.t

    #Sundials does not presereve problem shape...
    uReshaped = reshape(integrator.u,N,N,species)


    #Loop through all of the infectious cells
    for I in findall(uReshaped[:,:,end] .> 0.4)
        #Loop through all neighbors to see if they get infected
        for J in max(Ifirst, I-Ifirst):min(Ilast, I+Ifirst)
            #If the cell is healthy and not the current cell try to infect
            if isinf(cellsInfected[J]) && (I != J) && rand(chanceOfInfection)
                #Add Viral DNA to the cell
                uReshaped[J,2] = rand(Uniform(0.0,m2c(1e3)))
                #Mark that the cell was infected
                cellsInfected[J] = integrator.t
            end
        end
        #This cell has used its one attempt to infect neighboring cells
        infectFirstAttempt[I] = false
    end

    integrator.u = uReshaped[:]
end

cbVec = DiscreteCallback(VecCondition,VecAffect)
