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
        θcurrent.par[parChange] = θcurrent.par[parChange] * percent

        println(θcurrent.par[[20,23]])

        probKD = remake(prob; p=θcurrent)
        solKD = solve(probKD,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)

        IFNAve[i] = [mean(solKD[:,:,7,t]) for t in 1:length(solKD.t)]
        IFNStd[i] = [std(solKD[findall(u0[:,:,2] .> 0),7,t]) for t in 1:length(solKD.t)]
    end

    plt = plot()

    for (μ,σ,i) in zip(IFNAve,IFNStd,1:kdSamples)
        #ribbon =σ
        plot!(plt,tspan[1]:0.1:tspan[2],μ,ribbon =σ,label=percentLabels[i], fillalpha=0.3,
        legend=:right,framestyle = :box,legendtitle="Percent Knockdown",linewidth = 2,linealpha=1.0)
    end

    title!(statesNames[7])
    xlabel!("Time (hrs)")
    ylabel!("Average $(statesNames[7]) (nM)")
    xticks!(0:12:48)

    #return p
    return (solKD,plt)
end


plotKD = map(KD,[20,23])
savefig(plotKD[1][2],"../Figures/IRF7_KD.pdf")
savefig(plotKD[2][2],"../Figures/TREX1_KD.pdf")
