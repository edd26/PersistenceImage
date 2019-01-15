PersistenceImage
===

Package that provides an implementation to the persistence algorithm proposed in https://dl.acm.org/citation.cfm?id=3122017.

# Installation
In order to install the package follow the steps below

```julia
using Pkg
Pkg.add("https://github.com/chronchi/PersistenceImage.jl")
```

# Usage

Given a persistence diagram the algorithm will output a persistence image using the function ``transformdiagram``.

```julia
a = rand(100)
diagram = [a a.^(2/3)]
persimg = transformdiagram(diagram)
```

## Parameters

There are some parameters described in the paper that you have to choose, namely the grid size and distribution. The current implementation only supports the Gaussian distribution, therefore there is a choice of the
standard deviation `σ`.
The grid size is a tuple, namely `pixels` and is used as below.  

```julia
a = rand(100)
diagram = [a a.^(2/3)]
σ = 1.0 # standard deviation
pixels = (10, 10) # grid size
persimg = transformdiagram(diagram, pixels=pixels, σ=σ)
```
