# file with auxiliary functions to calculate the persistence image
using Statistics

"""
    transformdiagram(diagram)

Return the given persistence diagram to landscape.
"""
function tolandscape(diagram::Array{Float64, 2})
    newdiagram = copy(diagram)
    newdiagram[:, 2] -= newdiagram[:, 1]
    return newdiagram
end

function tolandscape2(diagram::Array{Float64, 2}; scaling::Float64=√2)
    newdiagram = copy(diagram)

    newdiagram[:, 1] = mean(diagram, dims=2)
    newdiagram[:, 2] = diagram[:, 2] - diagram[:, 1]
    newdiagram .*= scaling

    return newdiagram
end

"""
    weighting(x, y)

Return the weight of a point in the persistence landscape.
"""
function weighting(x, y)
    if y <= 0
        return 0
    elseif 0 < y && y < x
        return y/x
    else
        return 1
    end
end
