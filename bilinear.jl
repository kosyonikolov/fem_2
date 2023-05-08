function localMassMatrix(h::Number)
    return (h^2 / 9) * [1 1/2 1/2 1/4; 1/2 1 1/4 1/2; 1/2 1/4 1 1/2; 1/4 1/2 1/2 1]
end

# Doesn't depend on square size
function localStiffnessMatrix()
    return (1 / 6) * [4 -1 -1 -2; -1 4 -2 -1; -1 -2 4 -1; -2 -1 -1 4]
end