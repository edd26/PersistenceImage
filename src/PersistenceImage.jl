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
                img[size(img,1)+1-k,j] += weights * xpdf * ypdf
            end
        end
    end
    return img
end

"""
    transformdiagram2(diagram; pixels, σ, birth_range, death_range)
"""
function transformdiagram2(diagram::Array{Float64, 2};
                          pixels::Tuple{Int64,Int64}=(10,10),
                          σ::Float64 = 1.,
                          birth_range::Tuple{Float64, Float64}=(min(diagram[:,1]...),max(diagram[:,1]...)),
                          death_range::Tuple{Float64, Float64}=(min(diagram[:,2]...),max(diagram[:,2]...)),
                                                                                     )
    landscape = tolandscape2(diagram)
    births = diagram[:,1]
    persistances = diagram[:,2] - diagram[:,1]

    # create a grid over our landscape
    landscape_space = tolandscape2(diagram)
    space_diagram = [birth_range[1] death_range[1]
                        birth_range[1] death_range[2]
                        birth_range[2] death_range[2]]
    landscape_space = tolandscape2(space_diagram)
    x_start = min(landscape_space[:,1]...)
    x_end= max(landscape_space[:,1]...)
    y_start = min(landscape_space[:, 2]...)
    y_end = max(landscape_space[:, 2]...)

    xaxis = range(x_start, x_end, length=pixels[1])
    yaxis = range(y_start, y_end, length=pixels[1])
    img = zeros(Float64, pixels[1], pixels[2] )

    for i in 1:size(landscape, 1)
        birth = births[i]
        persistence = persistances[i]
        x_coord = landscape[i,1]
        y_coord = landscape[i,2]
        weights = weighting(birth, persistence)

        dB = Normal(x_coord, σ)
        dP = Normal(y_coord, σ)
        # get a vector of x values
        x_vals = pdf(dB, xaxis)
        y_vals = pdf(dP, yaxis)
        img_slice = (y_vals * (x_vals') * weights)[end:-1:1,:]
        img .+= img_slice
    end

    to_remove = findall(x -> x < eps(), img)
    for coord in to_remove
        img[coord] = 0
    end

    return img
end

end # module
