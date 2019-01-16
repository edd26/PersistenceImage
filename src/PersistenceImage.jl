module PersistenceImage


export transformdiagram, tolandscape, weighting

using Distributions

include("auxiliaryfunctions.jl")

"""
    transformdiagram(diagram; pixels, σ)

Return the persistence image of a persistence diagram, given the
persistence diagram, pixels and standard deviation.
"""
function transformdiagram(diagram::Array{Float64, 2};
                          pixels::Tuple{Int64,Int64}=(10,10), σ::Float64 = 1.)

    landscape = tolandscape(diagram)
    # min and max in the x axis of the landscape
    xmin = minimum([minimum(landscape[:,1]), 0])
    xmax = maximum(landscape[:,1])
    # min and max in the y axis of the landscape
    ymin = minimum([minimum(landscape[:,2]), 0])
    ymax = maximum(landscape[:,2])
    # create a grid over our landscape
    xaxis = range(xmin, stop=xmax, length=pixels[1]+1)
    yaxis = range(ymin, stop=ymax, length=pixels[2]+1)

    img = Array{Float64}(undef, pixels[1], pixels[2])

    for i in 1:size(landscape, 1)
        birth = landscape[i,1]
        persistence = landscape[i,2]
        weights = weighting(birth, persistence)
        dB = Normal(birth, σ)
        dP = Normal(persistence, σ)
        for k in 1:size(img, 1)
            for j in 1:size(img, 2)
                # here we get the center of each pixel in the grid
                xcenter = mean([xaxis[k], xaxis[k+1]])
                ycenter = mean([yaxis[j], yaxis[j+1]])
                # calculate the pdf in the center of the pixel in row k
                # column j
                xpdf = pdf(dB, xcenter)
                ypdf = pdf(dP, ycenter)
                # sum to the corresponding pixel
                img[size(img,1)+1-k,j] += weights * xpdf * xpdf
            end
        end
    end
    return img
end

end # module
