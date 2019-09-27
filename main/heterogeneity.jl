include("Driver.jl")
using ProgressMeter, DataFrames


#Give a unique parameter set for each cells (randomly choosen), ignore virus parameters?
percentRange =11
percents = range(0.01,1.0,length=percentRange)

virusPeaks = Vector(undef,length(percents))

@showprogress 1 "Computing..." for (i,percent) in enumerate(percents)
    global cellsInfected

    sampleDist = @. Uniform(θVals.*(1-percent),θVals.*(1+percent))
    θcurrent = reshape.(rand.(sampleDist,nCells),N,N)

    #Reset the infection matrix
    cellsInfected = fill(Inf,N,N) #Make constant when not testing
    cellsInfected[findall(u0[:,:,2] .> 0.0), 1] .= 0.0

    #Keep the virus parameters the same
    for par in [34,35]
    θcurrent[par] = fill(θVals[par],N,N)
    end

    #Contruct the ODEs
    probNew = remake(prob; p=θcurrent)
    #Solve the problem
    #sol = @time solve(prob,ROCK4(),saveat=0.1)
    sol = solve(probNew,CVODE_BDF(linear_solver=:GMRES),saveat=0.1,callback=cb_infect,tstops=tstop)
    virusPeaks[i] = [maximum(sol[coord,14,:]) for coord in cellIndicies]

end

#p=boxplot(vec.(virusPeaks),legend=false,outliers=false,frame=:box,color=:orange)
#p=violin(vec.(virusPeaks),legend=false,frame=:box,color=:orange)
#xticks!(1:percentRange,string.(percents))
#xlabel!("Parameter Variability")
#ylabel!("Maximum Viral Production")

#savefig(p,"./Figures/HeteroVirusProd.pdf")

#I have found no way to make a asymmetric violin plot besides converting the data
#into dataframes...so here I am...converting the data into dataframes...so ugly


#Flatten out the arrays
infectidx = u0[:,:,2] .> 0.0
#Seperate out the initial and secondary infections
primary = [ virusPeaksCurrent[infectidx] for virusPeaksCurrent in virusPeaks]
secondary = [ virusPeaksCurrent[.!infectidx] for virusPeaksCurrent in virusPeaks]

#Convert into dataframes
primary = DataFrame(percent=repeat(1:percentRange,inner = sum(infectidx)),vConc=vcat(primary...))
secondary = DataFrame(percent=repeat(1:percentRange,inner = sum(.!infectidx)),vConc=vcat(secondary...))

#Plot the data
@df primary violin(:percent, :vConc,side=:left,frame=:box,label=:primary,legend=:topleft,show_median=true,trim=false)
@df secondary violin!(:percent, :vConc,side=:right,label=:secondary,show_median=true,trim=false)
xticks!((1:percentRange),string.(percents))
xlabel!("Parameter Variability")
ylabel!("Maximum Viral Production")

savefig("./Figures/HeteroVirusProd.pdf")
