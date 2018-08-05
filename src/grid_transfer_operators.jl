# transfer_operators.jl : tranform quantities from fine to coarse grid and back

# restriction
R̃(dir::Int,op::TransferKind,n::Int...) = R̃_reversed(dir,op,reverse(n)...)
R̃_reversed(dir::Int,op::TransferKind,n::Int...) = kron(R₁.(op,n[1:dir-1])...,speye(n[dir]-1,n[dir]-1),R₁.(op,n[dir+1:end])...)

# interpolation
P̃(dir::Int,op::TransferKind,n::Int...) = P̃_reversed(dir,op,reverse(n)...)
P̃_reversed(dir::Int,op::TransferKind,n::Int...) = kron(P₁.(op,2.*n[1:dir-1])...,speye(n[dir]-1,n[dir]-1),P₁.(op,2.*n[dir+1:end])...)
