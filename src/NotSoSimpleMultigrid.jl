module NotSoSimpleMultigrid

# dependencies 
using SimpleMultigrid, StaticArrays, LinearAlgebra, SparseArrays, Statistics, Printf

import SimpleMultigrid: stencil2mat, TransferKind, R₁, P₁, Grid, MultigridCycle, V, W, F, Smoother, zero_x, MultigridIterable, μ_cycle!, F_cycle!, smooth!, residu, factor_twos

# export statements
export laplace2d, elliptic2d, anisotropic2d, rotated_anisotropic2d, laplace3d, elliptic3d, anisotropic3d
export Injection, FullWeighting, Cubic
export GaussSeidel, Jacobi
export coarsen
export MultigridIterable, V, W, F, V_cycle, W_cycle, F_cycle, \

# include source files
include("compose_matrix.jl")
include("grid_transfer_operators.jl")
include("grids.jl")
include("multigrid.jl")

end # module
