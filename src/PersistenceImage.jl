module PersistenceImage


export transformdiagram, tolandscape, weighting
export transformdiagram2

using Distributions

include("auxiliaryfunctions.jl")

"""
    transformdiagram(diagram; pixels, σ)

Return the persistence image of a persistence diagram, given the
persistence diagram, pixels and standard deviation.
"""
function transformdiagram(diagram::Array{Float64,2}; pixels::Tuple{Int64,Int64}=(10, 10), σ::Float64=1.0)

    landscape = tolandscape(diagram)
    # min and max in the x axis of the landscape
    xmin = minimum([minimum(landscape[:, 1]), 0])
    xmax = maximum(landscape[:, 1])
    # min and max in the y axis of the landscape
    ymin = minimum([minimum(landscape[:, 2]), 0])
    ymax = maximum(landscape[:, 2])
    # create a grid over our landscape
    xaxis = range(xmin, stop=xmax, length=pixels[1] + 1)
    yaxis = range(ymin, stop=ymax, length=pixels[2] + 1)

    # img = Array{Float64}(undef, pixels[1], pixels[2])
    img = zeros(pixels[1], pixels[2])

    for i in 1:size(landscape, 1)
        birth = landscape[i, 1]
        persistence = landscape[i, 2]
        weights = weighting(birth, persistence) # Change this to maxiaml persistence for all data, not local persistence of current point
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
                img[size(img, 1)+1-k, j] += weights * xpdf * ypdf
            end
        end
    end
    return img
end

"""
    transformdiagram2(diagram; pixels, σ, birth_range, death_range)
"""
function transformdiagram2(diagram::Array{Float64,2}; pixels::Tuple{Int64,Int64}=(10, 10), σ::Float64=1.0,
    birth_range::Tuple{Float64,Float64}=(min(0, diagram[:, 1]...), max(diagram[:, 1]...)),
    persistence_range::Tuple{Float64,Float64}=(min(0, diagram[:, 2]...), max(diagram[:, 2]...)),
    weight_with_maxp::Bool=false
)
    landscape = tolandscape(diagram)
    total_points = size(landscape, 1)
    if weight_with_maxp
        max_persistence = maximum(landscape[:, 2])[1]
    end

    # create a grid over our landscape
    x_start = birth_range[1]
    x_end = birth_range[end]
    y_start = persistence_range[1]
    y_end = persistence_range[end]

    xaxis = range(x_start, x_end, length=pixels[1] + 1)
    yaxis = range(y_start, y_end, length=pixels[2] + 1)

    xcentres = [mean([xaxis[k], xaxis[k+1]]) for k in 1:pixels[1]]
    ycentres = [mean([yaxis[k], yaxis[k+1]]) for k in 1:pixels[1]]

    img = zeros(Float64, pixels[1], pixels[2])
    for i in 1:total_points
        birth = landscape[i, 1]
        persistence = landscape[i, 2]

        if weight_with_maxp
            weights = weighting(birth, max_persistence)
        else
            weights = weighting(birth, persistence) # Change this to maxiaml persistence for all data, not local persistence of current point
        end

        xpdf(x) = Distributions.pdf(Normal(birth, σ), x)
        ypdf(y) = Distributions.pdf(Normal(persistence, σ), y)

        img_slice = (ypdf.(ycentres)*(xpdf.(xcentres)')*weights)[end:-1:1, :]
        img .+= img_slice
    end

    return img
end

end # module
