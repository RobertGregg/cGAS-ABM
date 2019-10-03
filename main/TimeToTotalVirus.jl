include("Driver.jl")

sol = @time solve(prob,CVODE_BDF(linear_solver=:GMRES),callback=cb_14,tstops=tstop)

f(x) =  1 - sum(sol(x)[:,:,2] .== 0.0)/N^2
plot(tstop,f.(tstop),framestyle = :box,legend=false)
xticks!(0:12:48)
ylims!((0,1))
xlabel!("Time (hrs)")
ylabel!("Infected Cell %")


savefig("../Figures/HeteroVirus.pdf")

cellsInfected_plot = copy(cellsInfected)
cellsInfected_plot[isinf.(cellsInfected_plot)] .= -1.0
heatmap(cellsInfected_plot,color=:viridis, clims=(0,10))
title!("Time of Infection")

savefig("../Figures/InfectionTime.pdf")
