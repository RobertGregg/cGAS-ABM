#Additional Files and Packages
include("Driver.jl")
using ProgressMeter, DataFrames

#Give a unique parameter set for each cells (randomly choosen), ignore virus parameters?
percentRange = 11
percents = range(0.0,1.0,length=percentRange)

sols = Vector(undef,length(percents))
u0New = copy(u0)

@showprogress 1 "Computing..." for (i,percent) in enumerate(percents)

    for (idx,val,θidx) in zip(nonZeroSpeciesIdx,nonZeroSpeciesValues,1:3)
        #Update the non-zero species
        u0New[:,:,idx] = rand(TruncatedNormal(val,percent,0.0,Inf),N,N)
        #Change the total amount of non-zero species inside the cell
        θ.mass[θidx] = u0New[:,:,idx]
    end

    #Contruct the ODEs
    probNew = remake(prob; u0=u0New,p=θ)
    #Solve the problem
    sol = solve(probNew,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)
    sols[i] = [maximum(sol[coord,7,:]) for coord in cellIndicies]

end

#Determine what cells were initiall infected
infectidx = u0[:,:,2] .> 0.0
#Seperate out the initial and secondary infections
primary = [ solsCurrent[infectidx] for solsCurrent in sols]
secondary = [ solsCurrent[.!infectidx] for solsCurrent in sols]

#Convert into dataframes
primary = DataFrame(percent=repeat(1:percentRange,inner = sum(infectidx)),IFN=vcat(primary...),CellState = fill("Primary",length(vcat(primary...))))
secondary = DataFrame(percent=repeat(1:percentRange,inner = sum(.!infectidx)),IFN=vcat(secondary...),CellState = fill("Secondary",length(vcat(secondary...))))

#Plot the data
@df primary violin(:percent, :vConc,side=:left,frame=:box,label=:primary,legend=:topleft,show_median=true)
@df secondary violin!(:percent, :vConc,side=:right,label=:secondary,show_median=true)
xticks!((1:percentRange),string.(round.(percents,digits=2)))
xlabel!("Initial Condition Variance")
ylabel!("Peak Interferon Production")

savefig("../Figures/InitialCondVary.pdf")


percent =0.05
for (idx,val,θidx) in zip(nonZeroSpeciesIdx,nonZeroSpeciesValues,1:3)
    #Update the non-zero species
    u0New[:,:,idx] = rand(TruncatedNormal(val,percent,0.0,Inf),N,N)
    #Change the total amount of non-zero species inside the cell
    θ.mass[θidx] = u0New[:,:,idx]
end

#Contruct the ODEs
probNew = remake(prob; u0=u0New,p=θ)
#Solve the problem
sol = solve(probNew,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)
