using ProgressMeter

include("Driver.jl")

solutionArray = []
parRange = 10

kcat8Vals = range(1.0,100.0,length=parRange)
τ7Vals = exp10.(range(3.0,6.0,length=parRange))

@showprogress 1 "Computing..." for kcat8 in kcat8Vals, τ7 in τ7Vals
    θVals[13] = kcat8
    θVals[27] = τ7
    remake(prob; p=θVals)

    push!(solutionArray,solve(prob,CVODE_BDF(linear_solver=:GMRES),save_everystep=false))
end

#Reduce the solution
finalIFNβ = [mean(sol[:,:,7,end]) for sol in solutionArray]

p = heatmap(log10.(τ7Vals),kcat8Vals,finalIFNβ,color=:viridis,framestyle=:box)

title!("Chronic Inflammation")
xlabel!("Interferon Degradation log(tau7)")
ylabel!("JAK/STAT Activity kcat8")

savefig(p,"./Figures/IFNb_Chronic.pdf")
