using JLD2, FileIO

@load "saveIFN.jld2" saveIFN

heatAll = DataFrame()
heatAll.x = repeat(repeat(1:N,N),timeLength)
heatAll.y = repeat(repeat(1:N,inner=N),timeLength)
heatAll.t = repeat(t,inner=N^2)
heatAll.IFN = saveIFN
