using BenchmarkTools


function doubler(vec)
    for i=1:length(vec)
        vec[i] *= 2.0
    end
    return vec
end

function doublerbounds(vec)
    @inbounds for i=1:length(vec)
        vec[i] *= 2.0
    end
    return vec
end

a = rand(100_000)

@benchmark doubler(a)




@benchmark doublerbounds(a)
