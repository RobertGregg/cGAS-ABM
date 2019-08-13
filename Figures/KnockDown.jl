
TrexKD=1:-0.25:0
percentLabels = [string(i)*"%" for i=0:25:100]

p = plot()

for percent in TrexKD
    θcurrent = prob.p
    θcurrent[23] *= percent
    probKD = remake(prob; p=θcurrent)

    solKD = solve(prob,ESERK5(),saveat=0.1)

    IFNAve = [mean(sol[:,:,statePlot,i]) for i in 1:length(sol.t)]
    plot!(p,)
end


plot(sol.t,IFNAve,legend=false)
title!(statesNames[statePlot])
xlabel!("Time (hrs)")
ylabel!("Average $(statesNames[statePlot]) (nM)")
