function kahanSum(sum, x, comp)
    add = x - comp
    s1 = sum + add
    comp = (s1 - sum) - add
    sum = s1
    return sum, comp
end