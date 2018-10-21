# transfer_operators.jl : tranform quantities from fine to coarse grid and back
Base.broadcastable(scaling::UniformScaling) = Ref(scaling)

# restriction
R̃(dir::Int,op::TransferKind,n::Int...) = kron(sparse.(I,n[end:-1:dir+1].-1,n[end:-1:dir+1].-1)...,R₁(op,n[dir]),sparse.(I,n[dir-1:-1:1].-1,n[dir-1:-1:1].-1)...)

# interpolation
P̃(dir::Int,op::TransferKind,n::Int...) = kron(sparse.(I,n[end:-1:dir+1].-1,n[end:-1:dir+1].-1)...,P₁(op,n[dir]),sparse.(I,n[dir-1:-1:1].-1,n[dir-1:-1:1].-1)...)
