include("Driver.jl")

function KD(parChange)

    KnockDownVals=1:-0.25:0
    percentLabels = [string(i)*"%" for i=0:25:100]
    kdSamples = length(KnockDownVals)

    IFNAve = Vector(undef,kdSamples)
    IFNStd = Vector(undef,kdSamples)
    solKD = []

    for (i,percent) in enumerate(KnockDownVals)
        θcurrent = deepcopy(θ)
        θcurrent.par[parChange] .= θcurrent.par[parChange] .* percent
        println(mean(θcurrent.par[parChange]))

        probKD = remake(prob; p=θcurrent)
        solKD = solve(probKD,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)

        IFNAve[i] = [mean(solKD[:,:,7,t]) for t in 1:length(solKD.t)]
        IFNStd[i] = [std(solKD[findall(u0[:,:,2] .> 0),7,t]) for t in 1:length(solKD.t)]
        IFNall[i] = [solKD[coord,7,:] for coord in cellIndicies]
    end

    p = plot()
    linCol = [:red,:gold,:green,:blue,:purple]

    for (IFNncurr,i) in zip(IFNall,1:kdSamples)
        plot!(p,tspan[1]:0.1:tspan[2],IFNncurr[:],label=percentLabels[i],
        legend=:false,framestyle = :box,
        linewidth = 2,linecolor=linCol[i],linealpha = 0.1)
    end

    title!(statesNames[7])
    xlabel!("Time (hrs)")
    ylabel!("$(statesNames[7]) (nM)")
    xticks!(0:12:48)

    #return p
    return (solKD,p)
end


plotKD = map(KD,[23])
plotKD[1][2]
