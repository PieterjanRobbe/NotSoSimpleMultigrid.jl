# NotSoSimpleMultigrid

| **Build Status** | **Coverage** |
|------------------|--------------|
| [![Build Status](https://travis-ci.org/PieterjanRobbe/NotSoSimpleMultigrid.jl.png)](https://travis-ci.org/PieterjanRobbe/NotSoSimpleMultigrid.jl) [![Build status](https://ci.appveyor.com/api/projects/status/cglp5y9k6k3j5sj2?svg=true)](https://ci.appveyor.com/project/PieterjanRobbe/notsosimplemultigrid-jl) | [![Coverage](https://codecov.io/gh/PieterjanRobbe/NotSoSimpleMultigrid.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PieterjanRobbe/NotSoSimpleMultigrid.jl) [![Coverage Status](https://coveralls.io/repos/github/PieterjanRobbe/NotSoSimpleMultigrid.jl/badge.svg)](https://coveralls.io/github/PieterjanRobbe/NotSoSimpleMultigrid.jl) |

Simple implementation of geometric semi-coarsened multigrid methods. 

## Installation

From the Julia REPL, type ] to enter Pkg mode and run

```julia
pkg> add https://github.com/PieterjanRobbe/NotSoSimpleMultigrid.jl
```

Load the package in Julia by

```julia
julia> using NotSoSimpleMultigrid
```

## Details

Basic implementation of a Geometric Multiple Semicoarsened Multigrid V-cycle, W-cycle and FMG-cycle.

Implemented smoothers use the `jacobi!` and `gauss_seidel!` methods from [IterativeSolvers.jl](https://github.com/JuliaMath/IterativeSolvers.jl).

Uses a matrix-dependent prolongation strategy to combine the corrections.
