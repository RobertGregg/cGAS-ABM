########################################################################
#Create an animation of Interferon production across the cell population
########################################################################
    #Get the max value for IFNβ concentration for whole simulation
    Maxval = maximum(sol[:,:,2,:])
    #Loop through time points to make an animation
    anim = @animate for i = tspan[1]:0.1:tspan[2]
        heatmap(sol(i)[:,:,2],clims=(0.0,Maxval),
                title="Time = " * string(i) * " hrs")
    end

    gif(anim,"../Figures/DNAAnimation.gif")



########################################################################
#Cross section of IFNb and Virus propagation
########################################################################

Maxval = maximum(sol[:,:,[7,14],:]) #largest value of IFN or Virus
anim = @animate for i = tspan[1]:0.1:tspan[2]
    ifn = sol(i)[:,:,7]
    vir = sol(i)[:,:,14]
    plot(diag(ifn),ylims=(0,Maxval),title="Time = " * string(i) * " hrs",legend=false)
    plot!(diag(vir))
end

gif(anim,"../Figures/CrossSectionAnimation.gif")



########################################################################
#Plot dynamics for a particular state
########################################################################

stateToPlot = 13
plotState=[sol[coord,stateToPlot,:] for coord in cellIndicies]
plot(sol.t,plotState[:],leg=false,framestyle=:box,xtickfontsize=16,ytickfontsize=16)
title!(statesNames[stateToPlot])
xlabel!("Time (hrs)")
ylabel!("(nM)")
savefig("../Figures/" * statesNames[stateToPlot] * "Dynamics.png")






########################################################################
#Plot dynamics for all states
########################################################################

allStates = Vector(undef,species)
for (i,name) in enumerate(statesNames)

    plotStateInfected = Vector(undef,sum(u0[:,:,2] .> 0.0))
    plotStateHealthy = Vector(undef,sum(u0[:,:,2] .== 0.0))
    j=1
    k=1
    for coord in cellIndicies
        if u0[coord,2] > 0.0
            plotStateInfected[j] = sol[coord,i,:]
            j += 1
        else
            plotStateHealthy[k] = sol[coord,i,:]
            k += 1
        end
    end
    allStates[i]=plot(sol.t,mean(plotStateInfected[:]),ribbon=std(plotStateInfected[:]),leg=false,framestyle=:box)
    plot!(sol.t,mean(plotStateHealthy[:]),ribbon=std(plotStateHealthy[:]),leg=false,framestyle=:box)
    title!(statesNames[i])
    xlabel!("Time (hrs)")
    ylabel!("(nM)")
end

plot(allStates...,size=(1200,800))
savefig("../Figures/AllStatesSecondary.pdf")


########################################################################
#Plot cell States (needs cell states function from StateCount.jl)
########################################################################

saveTimePoints = range(tspan[1],tspan[2],step=0.1)
allCellStates = zeros(Int64,length(saveTimePoints),3)

for k=1:length(saveTimePoints)
    allCellStates[k,:] = cellStates(saveTimePoints[k],θ)
end

plot(saveTimePoints,allCellStates,framestyle = :box,linewidth=2,
    labels=[:Healthy,:Infected,:Dead])
