using Plots

using Distributions
include("../src/PersistenceImage.jl")
include("../src/auxiliaryfunctions.jl")
# include("auxiliaryfunctions.jl")
import .PersistenceImage: transformdiagram as transformdiagram
import .PersistenceImage: transformdiagram2 as transformdiagram2
import .PersistenceImage: toalndscape2 as toalndscape2
import .PersistenceImage: toalndscape as toalndscape

barcodes1 = [
    0.1 0.2 # early born, short lived
    0.1 0.9 # early born, long lived
    0.5 0.6 # mid born, short lived
    0.8 0.9 # alte born, short lived
] # triangle-like shape 


function plot_bd_diagram(barcodes::Vector; dims=1:size(barcodes, 2),
    normilised_diagonal::Bool=true,
    alpha=0.4,
    kwargs...)

    plot_ref = plot(; xlims=(0, 1), ylims=(0, 1))# kwargs...)
    # Add diagonal
    if normilised_diagonal
        max_coord = 1
    else
        max_x = max([k for k in vcat([barcodes[d][:, 1] for (d, dim) in dims |> enumerate]...) if !isinf(k)]...)
        max_y = max([k for k in vcat([barcodes[d][:, 2] for (d, dim) in dims |> enumerate]...) if !isinf(k)]...)
        max_coord = max(max_x, max_y)
    end

    scaling_factor = 1.05
    min_val = -0.05
    plot!([0, scaling_factor * max_coord], [0, scaling_factor * max_coord], label="")
    xlims!(min_val, scaling_factor * max_coord)
    ylims!(min_val, scaling_factor * max_coord)

    for (p, dim) in enumerate(dims)
        my_vec = barcodes[p]

        args = (
            label="β$(dim)",
            aspect_ratio=1,
            size=(600, 600),
            legend=:bottomright,
            framestyle=:origin,
            alpha=alpha,
            kwargs...)

        plot!(my_vec[:, 1], my_vec[:, 2], seriestype=:scatter; args...)
    end

    xlabel!("birth")
    ylabel!("death")

    return plot_ref
end

function get_comparison_plt(barcodes1, grid_size, σ, birth_range, death_range)
    reformatted_barcodes = [vcat([hcat([barcodes1[k, 1], barcodes1[k, 2]]...) for k in 1:size(barcodes1, 1)]...)]
    bd_diagram = plot_bd_diagram(reformatted_barcodes)
    xticks!(0:0.1:1)
    yticks!(0:0.1:1)

    # ===-
    pers_image1 = transformdiagram(barcodes1, pixels=grid_size, σ=σ,)
    pers_image2_v1 = transformdiagram2(barcodes1, pixels=grid_size, σ=σ,)
    pers_image2_v2 = PersistenceImage.transformdiagram2(
        barcodes1,
        pixels=grid_size,
        σ=σ,
        birth_range=birth_range,
        death_range=death_range,
    )

    # ===- 
    heatmapargs = (aspectratio=1, cbar=:outerbottom)
    wid = 600
    hei = 500
    all_plots = [
        bd_diagram,
        Plots.heatmap(pers_image1; title="original", heatmapargs...),
        # plot(Gray.(pers_image1 ./ findmax(pers_image1)[1])),
        Plots.heatmap(pers_image2_v1; title="mod, data range", heatmapargs...),
        Plots.heatmap(pers_image2_v2; title="mod, range chagned", heatmapargs...),
        plot(Gray.(pers_image2_v2 ./ findmax(pers_image2_v2)[1]))
    ]
    total_plots = length(all_plots)

    plot(all_plots...;
        layout=(1, total_plots),
        size=(total_plots * wid, hei),
        plot_title="σ=$(σ), grid_size=($(grid_size[1]),$(grid_size[2]))"
    )
end

@testset "Transformation arguments set nr 1" begin
    grid_size = (100, 100)
    σ = 0.04 # standard deviation
    min_death = 0.0
    max_death = 1.0
    min_birth = 0.0
    max_birth = 1.0
    birth_range = (min_birth, max_birth)
    death_range = (min_death, max_death)

    pl1 = get_comparison_plt(barcodes1, grid_size, σ, birth_range, death_range)

    # ===-
    grid_size = (100, 100)
    σ = 0.02 # standard deviation
    min_death = 0.0
    max_death = 1.0
    min_birth = 0.0
    max_birth = 1.0
    birth_range = (min_birth, max_birth)
    death_range = (min_death, max_death)

    pl2 = get_comparison_plt(barcodes1, grid_size, σ, birth_range, death_range)
    # ===-
    grid_size = (1000, 1000)
    σ = 0.02 # standard deviation
    min_death = 0.0
    max_death = 1.0
    min_birth = 0.0
    max_birth = 1.0
    birth_range = (min_birth, max_birth)
    death_range = (min_death, max_death)

    pl3 = get_comparison_plt(barcodes1, grid_size, σ, birth_range, death_range)

    # ===-
    grid_size = (100, 100)
    σ = 0.05 # standard deviation
    min_death = 0.0
    max_death = 1.0
    min_birth = 0.0
    max_birth = 1.0
    birth_range = (min_birth, max_birth)
    death_range = (min_death, max_death)

    pl4 = get_comparison_plt(barcodes1, grid_size, σ, birth_range, death_range)

    # ===-
    grid_size = (10, 10)
    σ = 0.1 # standard deviation
    min_death = 0.0
    max_death = 1.0
    min_birth = 0.0
    max_birth = 1.0
    birth_range = (min_birth, max_birth)
    death_range = (min_death, max_death)

    pl5 = get_comparison_plt(barcodes1, grid_size, σ, birth_range, death_range)

    wid = 600
    hei = 500

    pl_collection = [pl1, pl2, pl3, pl4, pl5]
    tot = pl_collection |> length
    final_comparison_plt = plot(pl_collection...;
        layout=(tot, 1),
        size=(5 * wid, tot * hei)
    )
    savefig(final_comparison_plt, "results_comparison.png")
end