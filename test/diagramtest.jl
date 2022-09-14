# file test/diagramtest.jl

srcroot = "$(dirname(@__FILE__))"

using PersistenceImage
using Random
# using Test
# include("src/PersistenceImage.jl")
# import .PersistenceImage


# tests the ProteinPersistent module with proteins from the PDB archive
@testset "Set 1" begin
    a = rand(100)
    diagram = [a a .^ (2 / 3)]

    persimg = transformdiagram(diagram)
    @test typeof(persimg) == Array{Float64,2}
end


@testset "Set 2" begin
    a = rand(MersenneTwister(0),100)
    diagram = [a a .^ (2 / 3)]

    σ = 2.0
    pixels = (25, 25)
    persimg = transformdiagram(diagram, pixels = pixels, σ = σ)
    @test typeof(persimg) == Array{Float64,2}
end


@testset "NaN presence test" begin
    a = rand(MersenneTwister(0),100)
    diagram = [a a .^ (2 / 3)]

    σ = 2.0
    pixels = (25, 25)
    persimg = transformdiagram(diagram, pixels = pixels, σ = σ)
    @test !(any(x-> isnan(x), persimg))

    a = rand(MersenneTwister(1),100)
    diagram = [a a .^ (2 / 3)]

    σ = 2.0
    pixels = (25, 25)
    persimg = transformdiagram(diagram, pixels = pixels, σ = σ)
    @test !(any(x-> isnan(x), persimg))

    a = rand(MersenneTwister(1),100)
    diagram = [a a .^ (2 / 3)]

    σ = 0.04
    pixels = (25, 25)
    persimg = transformdiagram(diagram, pixels = pixels, σ = σ)
    @test !(any(x-> isnan(x), persimg))
end
