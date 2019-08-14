include("Driver.jl")
using ProgressMeter


#Give a unique parameter set for each cells (randomly choosen), ignore virus parameters?
parcentRange = 10
percents = range(0.01,1.0,length=parcentRange)

virusPeaks = Vector(undef,length(percents))

@showprogress 1 "Computing..." for (i,percent) in enumerate(percents)
    sampleDist = @. Uniform((1-percent)*θVals,(1+percent)*θVals)
    θcurrent = reshape.(rand.(sampleDist,nCells),N,N)

    #Contruct the ODEs
    probNew = remake(prob; p=θcurrent)
    #Solve the problem
    #sol = @time solve(prob,ROCK4(),saveat=0.1)
    sol = solve(probNew,CVODE_BDF(linear_solver=:GMRES),saveat=0.1,callback=cb_infect,tstops=tstop)

    virusPeaks[i] = [maximum(sol[coord,14,:]) for coord in cellIndicies]

end

p=boxplot(vec.(virusPeaks),legend=false,outliers=false,frame=:box,color=:orange)
xticks!(1:parcentRange,string.(percents))
xlabel!("Parameter Variability")
ylabel!("Maximum Viral Production")

savefig(p,"./Figures/HeteroVirusProd.pdf")
