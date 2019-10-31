using RCall
using RDatasets
using DataFrames
@rlibrary ggplot2
@rlibrary ggpubr


d = DataFrame(v = [3,4,5], w = [5,6,7], x = [1,2,3], y = [4,5,6], z = [1,1,2])

ggplot(d, aes(x=:x,y=:y)) + geom_line()
ggsave("MakePlot.pdf")


#Save for plotting a state over time
allStateData = DataFrame()

for (i,st) in enumerate(statesNames)
    allStateData[!,Symbol(st)] = vec(sol[:,:,i,:])
end

allStateData.Cell = repeat(1:nCells,length(sol.t))
allStateData.Time = repeat(sol.t,inner=nCells)
allStateData.Infected = repeat(vec(u0[:,:,2] .> 0.0),length(sol.t))

CSV.write("allStateData.csv", allStateData)

#Save the simulation for heatmap
heat10hr = DataFrame()

heat10hr.x = repeat(1:N,N)
heat10hr.y = repeat(1:N,inner=N)
heat10hr.IFN = vec(sol(10.0)[:,:,7])

CSV.write("heat10hr.csv", heat10hr)
