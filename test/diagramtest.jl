# file test/diagramtest.jl

srcroot = "$(dirname(@__FILE__))"

using PersistenceImage
using Test

# tests the ProteinPersistent module with proteins from the PDB archive

a = rand(100)
diagram = [a a.^(2/3)]

persimg = transformdiagram(diagram)
@test typeof(persimg) == Array{Float64, 2}

σ = 2.0
pixels = (25,25)
persimg = transformdiagram(diagram, pixels=pixels, σ=σ)
@test typeof(persimg) == Array{Float64, 2} 
