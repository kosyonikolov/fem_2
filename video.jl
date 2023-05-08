using Plots

function dumpFrames(dir::AbstractString, ts::AbstractArray, qs::AbstractArray)
    if !isdir(dir)
        mkdir(dir)
    end

    nPts = length(qs[1])
    nAxis = round(Integer, sqrt(nPts))
    @assert nAxis^2 == nPts

    nQ = length(ts)
    for i = 1:nQ
        t = ts[i]
        p = surface(rotl90(reshape(qs[i], (nAxis,nAxis))), zlims=(-5,125), title = "$t", size = (1280, 720))
        savefig(p, "$dir/$i.png")
    end
end