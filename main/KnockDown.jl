include("Driver.jl")

function KD(parChange)

    KD=1:-0.25:0
    percentLabels = [string(i)*"%" for i=0:25:100]
    kdSamples = length(KD)

    IFNAve = Vector(undef,kdSamples)
    IFNStd = Vector(undef,kdSamples)
    #solKD = []
    for (i,percent) in enumerate(KD)
        θcurrent = copy(prob.p)
        θcurrent[parChange] *= percent
        probKD = remake(prob; p=θcurrent)

        solKD = solve(probKD,ESERK5(),saveat=0.1)

        IFNAve[i] = [mean(solKD[:,:,7,t]) for t in 1:length(solKD.t)]
        IFNStd[i] = [std(solKD[findall(u0[:,:,2] .> 0),7,t]) for t in 1:length(solKD.t)]
    end

    p = plot()

    for (μ,σ,i) in zip(IFNAve,IFNStd,1:kdSamples)
        plot!(p,tspan[1]:0.1:tspan[2],μ,label=percentLabels[i],legend=:right,framestyle = :box,legendtitle="Percent Knockdown",linewidth = 2)
    end

    title!(statesNames[7])
    xlabel!("Time (hrs)")
    ylabel!("Average $(statesNames[7]) (nM)")
    xticks!(0:12:48)

    return p
    #return (solKD,p)
end


plotKD = map(KD,[20,23])
savefig(plotKD[1],"./Figures/IRF7_KD.pdf")
savefig(plotKD[2],"./Figures/TREX1_KD.pdf")
