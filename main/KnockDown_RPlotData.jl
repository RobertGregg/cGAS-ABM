include("Driver.jl")
function KD(parChange)

    KnockDownVals=1:-0.25:0
    percentLabels = [string(i)*"%" for i=0:25:100]
    kdSamples = length(KnockDownVals)

    solKD =  Vector(undef,kdSamples)

    for (i,percent) in enumerate(KnockDownVals)
        θcurrent = deepcopy(θ)
        θcurrent.par[parChange] = θcurrent.par[parChange] * percent

        println(θcurrent.par[[20,23]])

        probKD = remake(prob; p=θcurrent)
        sol = solve(probKD,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)
        solKD[i] = vec(sol[:,:,7,:])
    end

    #return p
    return solKD
end


KDSols = map(KD,[20,23])

#Save the simulation for KD
KDData = DataFrame()

KDData.IRF7KD = vcat(KDSols[1]...)
KDData.TREXKD = vcat(KDSols[2]...)

timeLength = 481
KnockDownVals=1:-0.25:0
percentLabels = [string(i)*"%" for i=0:25:100]
kdSamples = length(KnockDownVals)

KDData.Cell = repeat(1:nCells,timeLength*kdSamples)
KDData.Time = repeat(0:0.1:48,inner=nCells,outer=kdSamples)
KDData.Percent = repeat(KnockDownVals,inner=nCells*timeLength)
