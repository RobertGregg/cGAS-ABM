#Get the max value for IFNb concentration for whole simulation
Maxval = maximum(sol[:,:,7,:])

#Loop through all time points to make an animation
anim = @animate for i = tspan[1]:0.1:tspan[2]
    G = sol(i)[:,:,7]
    heatmap(G,clims=(0,Maxval),title="Time = " * string(i) * " hrs")
end

gif(anim,"CellGridAnimation.gif")



stateToPlot = 14
plotState=[sol[coord,stateToPlot,:] for coord in cellIndicies]
plot(sol.t,plotState[:],leg=false)
title!(statesNames[stateToPlot])
xlabel!("Time (hrs)")
ylabel!("(nM)")




#Plotting just the infected cells

p = Vector(undef,species)
for (i,name) in enumerate(statesNames)

    plotState = []
    for coord in cellIndicies
        if cellsInfected[coord,1] == 1
            push!(plotState,sol[coord,i,:])
        end
    end
    p[i]=plot(sol.t,plotState[:],leg=false,framestyle=:box)
    title!(statesNames[i])
    xlabel!("Time (hrs)")
    ylabel!("(nM)")
end

plot(p...,size=(1200,800))
savefig("AllStates.pdf")
