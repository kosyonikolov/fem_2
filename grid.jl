# Make a regular grid on [0,1]^2 with square cells
function makeGrid(ptsPerAxis::Integer)
    xBase = range(0, 1, ptsPerAxis)
    nPts = ptsPerAxis^2
    elemPerAxis = ptsPerAxis - 1
    nElem = elemPerAxis^2

    pts = Array{Float64, 2}(undef, (nPts, 2))
    elems = Array{Integer, 2}(undef, (nElem, 4))
    for iy = 0:ptsPerAxis - 1, ix = 0:ptsPerAxis - 1
        ptIdx = iy * ptsPerAxis + ix + 1
        pts[ptIdx, :] = [xBase[ix + 1], xBase[iy + 1]]
        if iy < elemPerAxis && ix < elemPerAxis
            elIdx = iy * elemPerAxis + ix + 1
            elems[elIdx, :] = [0, 1, ptsPerAxis, ptsPerAxis + 1] .+ ptIdx
        end
    end

    return pts, elems
end