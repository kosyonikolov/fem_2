using SparseArrays
using LinearAlgebra

include("bilinear.jl")
include("grid.jl")
include("jacobi.jl")
include("kahan.jl")

# Constants for equation
struct QuasiFem
    h::Number
    pts::AbstractArray
    elements::AbstractArray
    m0::AbstractArray
    m1::AbstractArray
end

function setupEquation(ptsPerAxis::Integer) :: QuasiFem
    pts, elems = makeGrid(ptsPerAxis)
    h = pts[2,1] - pts[1,1]
    lm0 = localMassMatrix(h)
    lm1 = localStiffnessMatrix()
    # Assemble both matrices
    nPts = size(pts)[1]
    m0 = spzeros(nPts, nPts)
    m1 = spzeros(nPts, nPts)
    nElem = size(elems)[1]
    for i in 1:nElem
        ids = elems[i,:]
        m0[ids, ids] .+= lm0
        m1[ids, ids] .+= lm1
    end

    return QuasiFem(h, pts, elems, m0, m1)
end

function initialState(x)
    rad = x[1]^2 + x[2]^2
    return rad <= 1/16 ? 100 : 0
end

function initialQ(pts::AbstractArray)
    n = size(pts)[1]
    return [initialState(pts[i,:]) for i = 1:n]
end

# Kirchoff change
function g(x::Number)
    return x^2 / 2
end

# Solve the FEM problem using the explicit Euler method
function solveEuler(fem::QuasiFem, Tmax::Number, τ::Number)
    nT = ceil(Int, Tmax / τ)
    ts = range(0, Tmax, nT)
    τ = ts[2] - ts[1]

    nPts = size(fem.pts)[1]
    qs = Array{Float64, 2}(undef, (nPts, nT))
    qs[:,1] = initialQ(fem.pts)

    for i = 2:nT
        # qp = qs[[i-1]];
        # qs[[i]] = LinearSolve[m0, m0 . qp - tau * m1 . Map[g, qp]];
        prev = qs[:,i-1]
        rhs = fem.m0 * prev - τ .* (fem.m1 * map(g, prev))
        qs[:,i] = fem.m0 \ rhs
    end

    return ts, qs
end

# LHS of Crank-Nicolson ODE method for the equation
function crankNicolsonLhs(fem::QuasiFem, qPrev::AbstractArray, qCurr::AbstractArray, τ::Number)
    a = fem.m0 * ((qCurr - qPrev) ./ τ)
    b = fem.m1 * map(g, 0.5 .* (qCurr + qPrev))
    return a .+ b 
end

# Solve CN system using simple iteration
function simpleIteration(fem::QuasiFem, qPrev::AbstractArray, qCurr::AbstractArray, τ::Number, z;
                         maxIters::Integer = 10, minErr::Number = 1e-3, dbgPrint::Bool = false)
    err = x -> norm(crankNicolsonLhs(fem, qPrev, x, τ))
    lastErr = err(qCurr)
    if dbgPrint
        println("Start = $lastErr")
    end

    for i in 1:maxIters
        qNew = qPrev - τ .* (z * map(g, (qPrev + qCurr) ./ 2))
        newErr = err(qNew)
        if dbgPrint
            println("Iter $i = $newErr")
        end

        if newErr >= lastErr
            return qCurr
        elseif newErr < minErr
            return qNew
        end
        qCurr = qNew
        lastErr = newErr
    end
    return qCurr
end

function solveCnSimpleIteration(fem::QuasiFem, Tmax::Number, τ::Number)
    nT = ceil(Int, Tmax / τ)
    ts = range(0, Tmax, nT)
    τ = ts[2] - ts[1]

    nPts = size(fem.pts)[1]
    qs = Array{Float64, 2}(undef, (nPts, nT))
    qs[:,1] = initialQ(fem.pts)

    # M0^-1 * M1
    z = collect(fem.m0) \ collect(fem.m1)

    for i = 2:nT
        prev = qs[:,i-1]
        qs[:,i] = simpleIteration(fem, prev, prev, τ, z)
    end

    return ts, qs
end

# Solve CN system using Newton's method
# Calculate the Jacobi matrix only once
function newton1(fem::QuasiFem, qPrev::AbstractArray, qCurr::AbstractArray, τ::Number;
                 maxIters::Integer = 10, minErr::Number = 1e-3, 
                 maxItersIsFail::Bool = false, dbgPrint::Bool = false)
    f = x -> crankNicolsonLhs(fem, qPrev, x, τ)

    lastErr = norm(f(qCurr))
    if dbgPrint
        println("Start = $lastErr")
    end

    j = jacobi(f, qCurr)
    invJ = inv(j)

    for i in 1:maxIters
        qNew = qCurr - invJ * f(qCurr) 
        newErr = norm(f(qNew))
        if dbgPrint
            println("Iter $i = $newErr")
        end

        if newErr >= lastErr
            return qCurr, false
        elseif newErr < minErr
            return qNew, true
        end
        qCurr = qNew
        lastErr = newErr
    end
    return qCurr, !maxItersIsFail
end

# Solve CN system using Newton's method
# Calculate the Jacobi matrix at every iteration
function newton2(fem::QuasiFem, qPrev::AbstractArray, qCurr::AbstractArray, τ::Number;
                 maxIters::Integer = 10, minErr::Number = 1e-3, 
                 maxItersIsFail::Bool = false, dbgPrint::Bool = false)
    f = x -> crankNicolsonLhs(fem, qPrev, x, τ)

    lastErr = norm(f(qCurr))
    if dbgPrint
        println("Start = $lastErr")
    end

    for i in 1:maxIters
        j = jacobi(f, qCurr)
        qNew = j \ (j * qCurr - f(qCurr)) 
        newErr = norm(f(qNew))
        if dbgPrint
            println("Iter $i = $newErr")
        end

        if newErr >= lastErr
            return qCurr, false
        elseif newErr < minErr
            return qNew, true
        end
        qCurr = qNew
        lastErr = newErr
    end
    return qCurr, !maxItersIsFail
end

function solveCnNewton(fem::QuasiFem, Tmax::Number, τ::Number)
    nT = ceil(Int, Tmax / τ)
    ts = range(0, Tmax, nT)
    τ = ts[2] - ts[1]

    nPts = size(fem.pts)[1]
    qs = Array{Float64, 2}(undef, (nPts, nT))
    qs[:,1] = initialQ(fem.pts)

    for i = 2:nT
        prev = qs[:,i-1]
        qs[:,i], _ = newton1(fem, prev, copy(prev), τ)
    end

    return ts, qs
end

function solveCnNewtonAdaptive(fem::QuasiFem, Tmax::Number, τMax::Number; τMin::Number = 1e-7, λ::Number = 2)
    qs = [Float64.(initialQ(fem.pts))]
    steps = [] # Time steps
    ts = [0.0]

    t = 0.0
    tComp = 0.0
    τ = τMax

    goodCnt = 0
    goodLim = 2

    while t < Tmax
        prev = qs[end]
        cprev = copy(prev)
        curr = missing
        firstGood = true

        while true
            curr, ok = newton2(fem, prev, cprev, τ, minErr = 1e-2, maxIters = 5, maxItersIsFail = true)
            if ok || τ <= τMin
                break
            end
            firstGood = false
            τ = max(τMin, τ / λ)
        end

        push!(qs, curr)
        push!(steps, τ)
        # t += τ
        t, tComp = kahanSum(t, τ, tComp)
        push!(ts, t)

        # Update step
        if firstGood
            goodCnt += 1
            if goodCnt >= goodLim
                goodCnt = 0
                τ = min(τMax, τ * λ)
            end
        else
            goodCnt = 0
        end
    end

    return ts, qs, steps
end
