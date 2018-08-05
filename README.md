# NotSoSimpleMultigrid
[![Build Status](https://travis-ci.org/PieterjanRobbe/NotSoSimpleMultigrid.jl.png)](https://travis-ci.org/PieterjanRobbe/NotSoSimpleMultigrid.jl)

Simple implementation of geometric semi-coarsened multigrid methods. 

## Installation

```julia
Pkg.clone("https://github.com/PieterjanRobbe/NotSoSimpleMultigrid.jl")
```

Load the package in Julia by

```julia
using NotSoSimpleMultigrid
```

## Details

Basic implementation of a Geometric Multiple Semicoarsened Multigrid V-cycle, W-cycle and FMG-cycle.

Implemented smoothers use the `jacobi!` and `gauss_seidel!` methods from the [IterativeSolvers](https://github.com/JuliaMath/IterativeSolvers.jl) package.

Uses a matrix-dependent prolongation strategy to combine the corrections.
