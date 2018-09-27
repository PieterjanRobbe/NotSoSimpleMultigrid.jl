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
MultigridMethod(A::Union{AbstractMatrix,Function}, sz::NTuple, cycle_type::MultigridCycle; max_iter::Int=20, R_op::TransferKind=FullWeighting(), P_op::TransferKind=FullWeighting(), ngrids=factor_twos.(sz), smoother::Smoother=GaussSeidel()) = MultigridIterable(coarsen(A,sz,R_op,P_op,ngrids),max_iter,cycle_type,smoother,Float64[])

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
V_cycle(A::Union{AbstractMatrix,Function}, sz::NTuple; kwargs...) = MultigridMethod(A,sz,V(); kwargs...)

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
W_cycle(A::Union{AbstractMatrix,Function}, sz::NTuple; kwargs...) = MultigridMethod(A,sz,W(); kwargs...)

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
F_cycle(A::Union{AbstractMatrix,Function}, sz::NTuple; kwargs...) = MultigridMethod(A,sz,F(); max_iter=1, kwargs...)

# semicoarsening cycle
grids_at_level(sz::Tuple,ℓ::Int) = Base.Iterators.filter(idx->sum(idx.-1)==ℓ-1,Base.product(range.(1,sz)...))

children(idx) = Base.Iterators.filter(i->all(i[2].>0),enumerate(broadcast(-,idx,δ.(1:length(idx),length(idx)))))

parents(idx,sz) = Base.Iterators.filter(i->all(i[2]<=sz),enumerate(broadcast(+,idx,δ.(1:length(idx),length(idx)))))

function μ_cycle!(grids::Array{G} where {G<:Grid}, μ::Int, ν₁::Int, ν₂::Int, grid_ptr::Int, smoother::Smoother)
	d = ndims(grids)
	for idx in grids_at_level(size(grids),grid_ptr)
		smooth!(grids[idx...],ν₁,smoother)
	end
	if grid_ptr == sum(size(grids).-1)+1
		grids[end].x .= grids[end].A\grids[end].b # exact solve
	else
		for idx in grids_at_level(size(grids),grid_ptr+1)
			child_iter = Base.Iterators.filter(i->all(i[2].>=1),enumerate([idx.-δ(i,d) for i in 1:d]))
			grids[idx...].b .= mean(map(i->grids[last(i)...].R[first(i)]*residu(grids[last(i)...]),child_iter))
			grids[idx...].x .= zeros(grids[idx...].x)
		end
		μ_cycle!(grids,μ,ν₁,ν₂,grid_ptr+1,smoother)
		for idx in grids_at_level(size(grids),grid_ptr)
			parent_iter = Base.Iterators.filter(i->all(i[2].<=size(grids)),enumerate([idx.+δ(i,d) for i in 1:d]))
			# matrix-dependent prolongation
			λ = map(i->grids[idx...].A*high_freq_mode(first(i),grids[idx...].sz),parent_iter)
			λ² = broadcast(i->broadcast(j->j^2,i),λ)
			ω = map(i->λ²[i]./sum(λ²),1:length(λ)) # weight factors from [Naik, Van Rosendale]
			ip = map(i->grids[last(i)...].P[first(i)]*grids[last(i)...].x,parent_iter)
			grids[idx...].x .+= sum(map(i->ω[i].*ip[i],1:length(ω)))
			smooth!(grids[idx...],ν₂,smoother)
		end
	end
end

high_freq_mode(dir,sz) = vec(Int[-iseven(i[dir])+isodd(i[dir]) for i in Base.product(range.(1,sz.-1)...)])

function F_cycle!(grids::Array{G} where {G<:Grid}, ν₀::Int, ν₁::Int, ν₂::Int, grid_ptr::Int, smoother::Smoother)
	d = ndims(grids)
	if grid_ptr == sum(size(grids).-1)+1
		grids[grid_ptr].x .= zeros(grids[grid_ptr].x)
	else
		for idx in grids_at_level(size(grids),grid_ptr+1)
			child_iter = Base.Iterators.filter(i->all(i[2].>=1),enumerate([idx.-δ(i,d) for i in 1:d]))
			grids[idx...].b .= mean(map(i->grids[last(i)...].R[first(i)]*grids[last(i)...].b,child_iter))
		end
		F_cycle!(grids,ν₀,ν₁,ν₂,grid_ptr+1,smoother)
		for idx in grids_at_level(size(grids),grid_ptr)
			parent_iter = Base.Iterators.filter(i->all(i[2].<=size(grids)),enumerate([idx.+δ(i,d) for i in 1:d]))
			# matrix-dependent prolongation
			λ = map(i->grids[idx...].A*high_freq_mode(first(i),grids[idx...].sz),parent_iter)
			λ² = broadcast(i->broadcast(j->j^2,i),λ)
			ω = map(i->λ²[i]./sum(λ²),1:length(λ)) # weight factors from [Naik, Van Rosendale]
			ip = map(i->P̃(first(i),Cubic(),grids[last(i)...].sz...)*grids[last(i)...].x,parent_iter) # cubic
			grids[idx...].x .= sum(map(i->ω[i].*ip[i],1:length(ω)))
		end
	end
	for i in 1:ν₀
		μ_cycle!(grids,1,ν₁,ν₂,grid_ptr,smoother)
	end
end
