#Get the max value for IFNb concentration for whole simulation
Maxval = maximum(sol[:,:,7,:])

#Loop through all time points to make an animation
anim = @animate for i = tspan[1]:0.1:tspan[2]
    G = sol(i)[:,:,7]
    heatmap(G,clims=(0,Maxval),title="Time = " * string(i) * " hrs")
end

gif(anim,"CellGridAnimation.gif")
