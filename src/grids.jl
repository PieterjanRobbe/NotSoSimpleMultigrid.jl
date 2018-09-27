# grids.jl : introduces the Grid type for easy use of multigrid algorithms

# coarsen the problem A into a sequence of coarser grids
"""
    coarsen(A, sz, R_op, P_op, ngrids)

Coarsen the matrix `A` using `R_op` as restriction operator and `P_op` as interpolation operator. The matrix `A` is the discrete version of a PDE defined on an `sz`-point mesh, and `ngrids` is the number of levels. For MSG, `ngrids` is a tuple with the number of grids in each direction. A Galerkin approach is used to compose the coarse matrices, unless a function `f` is provided for direct discretization.
"""
function coarsen(A::SparseMatrixCSC, sz::NTuple, R_op::TransferKind, P_op::TransferKind, ngrids::NTuple)
    @assert length(sz) == length(ngrids)
    d = length(ngrids)
    grids = Array{Grid{typeof(A),Vector{eltype(A)},NTuple{d,typeof(A)},typeof(sz)},length(sz)}(ngrids...)
    grids[1] = Grid(A,zero_x(A),zero_x(A),tuple(R̃.(1:d,R_op,sz...)...),ntuple(i->spzeros(0,0),d),sz)
    for (idx,i) in enumerate(Base.Iterators.drop(Base.product(range.(1,ngrids)...),1))
        sz_c = sz.>>(i.-1)
        R_mats = tuple(R̃.(1:d,R_op,sz_c...)...)
        P_mats = tuple(P̃.(1:d,P_op,sz_c...)...)
        grid_num = i.-δ(indmax(i),d)
        A_c = grids[grid_num...].R[indmax(i)]*grids[grid_num...].A*P_mats[indmax(i)]
        grids[i...] = Grid(A_c,zero_x(A_c),zero_x(A_c),R_mats,P_mats,sz_c)
    end
    grids
end

function coarsen(f::Function, sz::NTuple, R_op::TransferKind, P_op::TransferKind, ngrids::NTuple)
    @assert length(sz) == length(ngrids)
    d = length(ngrids)
	A0 = f(sz...) 
    grids = Array{Grid{typeof(A0),Vector{eltype(A0)},NTuple{d,typeof(A0)},typeof(sz)},length(sz)}(ngrids...)
    grids[1] = Grid(A0,zero_x(A0),zero_x(A0),tuple(R̃.(1:d,R_op,sz...)...),ntuple(i->spzeros(0,0),d),sz)
    for (idx,i) in enumerate(Base.Iterators.drop(Base.product(range.(1,ngrids)...),1))
        sz_c = sz.>>(i.-1)
        R_mats = tuple(R̃.(1:d,R_op,sz_c...)...)
        P_mats = tuple(P̃.(1:d,P_op,sz_c...)...)
        grid_num = i.-δ(indmax(i),d)
		A_c = f(sz_c...)
        grids[i...] = Grid(A_c,zero_x(A_c),zero_x(A_c),R_mats,P_mats,sz_c)
    end
    grids
end

δ(i,n) = ntuple(j->j==i,n)
