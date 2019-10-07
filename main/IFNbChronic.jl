using ProgressMeter

include("Driver.jl")

parRange = 20
solutionArray = Vector(undef,parRange^2)

kcat8Vals = range(1.0,50.0,length=parRange)
τ7Vals = exp10.(range(0.0,2.0,length=parRange))

#Create an iterator to loop through the parameters (iterator loops through columns, need to switch τ7 and kcat8)
iter = zip(Iterators.product(τ7Vals, kcat8Vals),1:length(solutionArray))
#This could be parallelized if its too slow
@showprogress 1 "Computing..." for ((kcat8,τ7),i) in iter
    θ.par[13] = kcat8
    θ.par[27] = τ7
    probNew = remake(prob; p=θ)

    solutionArray[i] = solve(probNew,CVODE_BDF(linear_solver=:GMRES),save_everystep=false)
end

#Reduce the solution
finalIFNβ = [mean(sol[:,:,7,end]) for sol in solutionArray]

plt = heatmap(log10.(τ7Vals),kcat8Vals,reverse(finalIFNβ),color=:viridis,framestyle=:box)

title!("Chronic Inflammation")
xlabel!("Interferon Degradation tau7")
ylabel!("JAK/STAT Activity kcat8")

savefig(plt,"../Figures/IFNbChronic.pdf")
