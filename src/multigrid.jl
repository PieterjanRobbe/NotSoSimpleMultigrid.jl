# multigrid.jl : define multigrid cycle operations

"""
MultigridMethod(A, sz, cycle_type)
MultigridMethod(A, sz, cycle_type; kwargs...)
MultigridMethod(f, sz, cycle_type)
MultigridMethod(f, sz, cycle_type; kwargs...)

Geometric Multigrid method of type `cycle_type` for matrix `A`, that results from a discretization of a PDE on [0,1]^d using an `sz`-point mesh. A Galerkin approach is used to compose the coarse matrices, unless a function `f` is provided for direct discretization.

Inputs
======
A          : SparseMatrixCSC, sparse matrix with discretized PDE
f          : Function, function used for direct-discretization of the PDE
sz         : NTuple, PDE grid size, e.g., `(n,m)`
cycle_type : Multigrid cycle type, can be `V(ν₁,ν₂)` (V-cycle), `W(ν₁,ν₂)` (W-cycle) or `F(ν₀,ν₁,ν₂)` (FMG)

Options
=======
* max_iter : maximum number of iterations. Default value is `20` For FMG using the `F()`-cycle, this is set to 1.
* R_op     : restriction operator type, can be `Injection()` or `FullWeighting()` (default) 
* P_op     : interpolation operator type, can be `Injection()`, `FullWeighting()` (default), or `Cubic()` 
* ngrids   : total number of grids to use, default is `min.(⌊log₂(sz)⌋)`
* smoother : smoother, can be `GaussSeidel()` of `Jacobi()`
"""
MultigridMethod(A::Union{AbstractMatrix,Function}, sz::NTuple, cycle_type::MultigridCycle; max_iter::Int=20, R_op::TransferKind=FullWeighting(), P_op::TransferKind=FullWeighting(), ngrids=factor_twos.(sz), smoother::Smoother=GaussSeidel()) = MultigridIterable(coarsen(A, sz, R_op, P_op, ngrids), max_iter, cycle_type, smoother, Vector{Float64}(undef, 0))

"""
V_cycle(A, sz)
V_cycle(A, sz; kwargs...)

Geometric semicoarsened Multigrid V(2,1)-cycle for matrix `A`, that results from a discretization of a PDE on [0,1]^d using an `sz`-point mesh.

Inputs
======
A          : SparseMatrixCSC, sparse matrix with discretized PDE
sz         : NTuple, PDE grid size, e.g., `(n,m)`

For other options, see `MultigridMethod`.
"""
V_cycle(A::Union{AbstractMatrix,Function}, sz::NTuple; kwargs...) = MultigridMethod(A, sz, V(); kwargs...)

"""
W_cycle(A, sz)
W_cycle(A, sz; kwargs...)

Geometric semicoarsened Multigrid W(2,1)-cycle for matrix `A`, that results from a discretization of a PDE on [0,1]^d using an `sz`-point mesh.

Inputs
======
A          : SparseMatrixCSC, sparse matrix with discretized PDE
sz         : NTuple, PDE grid size, e.g., `(n,m)`

For other options, see `MultigridMethod`.
"""
W_cycle(A::Union{AbstractMatrix,Function}, sz::NTuple; kwargs...) = MultigridMethod(A, sz, W(); kwargs...)

"""
F_cycle(A, sz)
F_cycle(A, sz; kwargs...)

Geometric semicoarsened Full Multigrid for matrix `A`, that results from a discretization of a PDE on [0,1]^d using an `sz`-point mesh.

Inputs
======
A          : SparseMatrixCSC, sparse matrix with discretized PDE
sz         : NTuple, PDE grid size, e.g., `(n,m)`

For other options, see `MultigridMethod`.
"""
F_cycle(A::Union{AbstractMatrix,Function}, sz::NTuple; kwargs...) = MultigridMethod(A, sz, F(); max_iter=1, kwargs...)

# semicoarsening cycle
filter_at(R, s) = Base.Iterators.filter(i->sum(i.I)==s, R)

filter_fit(I1, Rall, R) = Base.Iterators.filter(i->i[2] ∈ Rall, enumerate(R))

child_iter(Rall, I1, I) = filter_fit(I1, Rall, Base.Iterators.reverse(filter_at(CartesianIndices(UnitRange.(Tuple(I-I1), Tuple(I))), sum(Tuple(I))-1)))

parent_iter(Rall, I1, I) = filter_fit(I1, Rall, filter_at(CartesianIndices(UnitRange.(Tuple(I), Tuple(I+I1))), sum(Tuple(I))+1))

grids_at_level(R::CartesianIndices{d}, s) where d = filter_at(R, s+d-1)

function μ_cycle!(grids::Array{G} where {G<:Grid}, μ::Int, ν₁::Int, ν₂::Int, grid_ptr::Int, smoother::Smoother)
    R = CartesianIndices(size(grids))
    I1, Iend = first(R), last(R)
    for I in grids_at_level(R, grid_ptr)
        smooth!(grids[I], ν₁, smoother)
    end
    if grid_ptr == sum(Tuple(Iend-I1)) + 1
        grids[Iend].x .= grids[Iend].A\grids[Iend].b # exact solve
    else
        for I in grids_at_level(R, grid_ptr+1)
            R_child = child_iter(R, I1, I)
            grids[I].b .= mean(map(i->grids[last(i)].R[first(i)]*residu(grids[last(i)]), R_child))
            fill!(grids[I].x, zero(eltype(grids[I].x)))
        end
        for i in 1:μ
            μ_cycle!(grids, μ, ν₁, ν₂, grid_ptr+1, smoother)
        end
        for I in grids_at_level(R, grid_ptr)
            R_parent = parent_iter(R, I1, I)
            # matrix-dependent prolongation
            λ = map(i->grids[I].A * high_freq_mode(first(i), grids[I].sz), R_parent)
            λ² = broadcast(i->broadcast(j->j^2, i), λ)
            ω = map(i->λ²[i]./sum(λ²), 1:length(λ)) # weight factors from [Naik, Van Rosendale]
            ip = map(i->grids[last(i)].P[first(i)]*grids[last(i)].x, R_parent)
            coarse_grid_correction!(grids, I, sum(map(i->ω[i].*ip[i], 1:length(ω))))
            smooth!(grids[I], ν₂, smoother)
        end
    end
end

high_freq_mode(dir,sz) = view(map(i->-iseven(i[dir])+isodd(i[dir]), CartesianIndices(sz.-1)), :)

function F_cycle!(grids::Array{G} where {G<:Grid}, ν₀::Int, ν₁::Int, ν₂::Int, grid_ptr::Int, smoother::Smoother)
    R = CartesianIndices(size(grids))
    I1, Iend = first(R), last(R)
    if grid_ptr == sum(Tuple(Iend-I1)) + 1
        fill!(grids[grid_ptr].x, zero(eltype(grids[grid_ptr].x)))
    else
        for I in grids_at_level(R, grid_ptr+1)
            R_child = child_iter(R, I1, I)
            grids[I].b .= mean(map(i->grids[last(i)].R[first(i)]*grids[last(i)].b, R_child))
        end
        F_cycle!(grids, ν₀, ν₁, ν₂, grid_ptr+1, smoother)
        for I in grids_at_level(R, grid_ptr)
            R_parent = parent_iter(R, I1, I)
            # matrix-dependent prolongation
            λ = map(i->grids[I].A * high_freq_mode(first(i), grids[I].sz), R_parent)
            λ² = broadcast(i->broadcast(j->j^2, i), λ)
            ω = map(i->λ²[i]./sum(λ²),1:length(λ)) # weight factors from [Naik, Van Rosendale]
            ip = map(i->P̃(first(i), Cubic(), grids[last(i)].sz...) * grids[last(i)].x, R_parent)
            grids[I].x .= sum(map(i->ω[i].*ip[i], 1:length(ω)))
        end
    end
    for i in 1:ν₀
        μ_cycle!(grids, 1, ν₁, ν₂, grid_ptr, smoother)
    end
end
