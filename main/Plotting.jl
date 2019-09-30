#Get the max value for IFNb concentration for whole simulation
Maxval = maximum(sol[:,:,7,:])

#Loop through all time points to make an animation
anim = @animate for i = tspan[1]:0.1:tspan[2]
    G = sol(i)[:,:,7]
    H = 10.0 .* sol(i)[:,:,14]
    plot(diag(G),ylims=(0,Maxval),title="Time = " * string(i) * " hrs",legend=false)
    plot!(diag(H))
end

gif(anim,"CellGridAnimation.gif")



Maxval = maximum(sol[:,:,7,:])
#Loop through all time points to make an animation
anim = @animate for i = tspan[1]:0.1:tspan[2]
    G = sol(i)[:,:,7]
    heatmap(G,clims=(0.0,Maxval),title="Time = " * string(i) * " hrs")
end

gif(anim,"CellGridAnimationlow.gif")




stateToPlot = 14
plotState=[sol[coord,stateToPlot,:] for coord in cellIndicies]
plot(sol.t,plotState[:],leg=false,framestyle=:box)
title!(statesNames[stateToPlot])
xlabel!("Time (hrs)")
ylabel!("(nM)")






#Plotting just the infected cells

plt = Vector(undef,species)
for (i,name) in enumerate(statesNames)

    plotState = Vector(undef,sum(u0[:,:,2] .> 0.0))
    j=1
    for coord in cellIndicies
        if u0[coord,2] > 0.0
            plotState[j] = sol[coord,i,:]
            j += 1
        end
    end
    plt[i]=plot(sol.t,plotState[:],leg=false,framestyle=:box)
    title!(statesNames[i])
    xlabel!("Time (hrs)")
    ylabel!("(nM)")
end

plot(plt...,size=(1200,800))
savefig("AllStates.pdf")
