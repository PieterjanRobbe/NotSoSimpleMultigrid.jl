# transfer_operators.jl : tranform quantities from fine to coarse grid and back

# restriction
R̃(dir::Int,op::TransferKind,n::Int...) = kron(speye.(n[end:-1:dir+1].-1,n[end:-1:dir+1].-1)...,R₁(op,n[dir]),speye.(n[dir-1:-1:1].-1,n[dir-1:-1:1].-1)...)

# interpolation
P̃(dir::Int,op::TransferKind,n::Int...) = kron(speye.(n[end:-1:dir+1].-1,n[end:-1:dir+1].-1)...,P₁(op,n[dir]),speye.(n[dir-1:-1:1].-1,n[dir-1:-1:1].-1)...)
#P̃(dir::Int,op::TransferKind,n::Int...) = kron(speye.(n[dir+1:end].-1,n[dir+1:end].-1)...,P₁(op,n[dir]),speye.(n[1:dir-1].-1,n[1:dir-1].-1)...)
