#Run with:
# tspan[2] = 168
# Wash
# no par vary
using ProgressMeter

include("Driver.jl")

parRange = 20
solutionArray = Vector(undef,parRange^2)

#Vary STAT production and IFN Degradation
kcat8True,τ7True = θ.par[[13,27]]
kcat8Vals = range(0.5*kcat8True,2.0*kcat8True,length=parRange)
τ7Vals = range(0.5*τ7True,2.0*τ7True,length=parRange)

#Create an iterator to loop through the parameters (iterator loops through columns, need to switch τ7 and kcat8)
iter = zip(Iterators.product(kcat8Vals,τ7Vals),1:length(solutionArray))
#This could be parallelized if its too slow
for ((kcat8,τ7),i) in iter
    @show kcat8,τ7
    θ.par[13] = kcat8
    θ.par[27] = τ7
    probNew = remake(prob; p=θ)

    solutionArray[i] = solve(probNew,CVODE_BDF(linear_solver=:GMRES),saveat=1.0)
    meanIFNβ = mean(solutionArray[i][:,:,7,end])
    @show meanIFNβ
end

#Reduce the solution
finalIFNb = [mean(sol[:,:,7,end]) for sol in solutionArray]

plt = heatmap(kcat8Vals,τ7Vals,reverse(finalIFNb),color=:viridis)

title!("Chronic Inflammation")
xlabel!("Interferon Degradation tau7")
ylabel!("JAK/STAT Activity kcat8")

savefig(plt,"../Figures/IFNbChronic1.pdf")


@save "savefinalIFNb.jld2" finalIFNb
#=
θ.par[[13,27]] = [20.,20.]

probNew = remake(prob; p=θ)

sol = solve(probNew,CVODE_BDF(linear_solver=:GMRES))
=#
