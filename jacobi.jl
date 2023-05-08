function jacobi(f::Function, x::AbstractArray; ε::Number = 1e-3)
    y0 = f(x)
    nX = lastindex(x)
    nY = lastindex(y0)
    j = zeros(nY, nX)
    for i in 1:nX
        old = x[i]
        x[i] += ε
        j[:, i] = (f(x) .- y0) ./ ε
        x[i] = old
    end
    return j
end

function fTest(p::AbstractArray)
    x = p[1]
    y = p[2]
    return [x^2 + 3y + 2 * x * y, y^3 - 3x^2 + x]
end