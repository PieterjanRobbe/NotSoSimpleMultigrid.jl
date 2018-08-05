# test_semicoarsening.jl : test semicoarsening method

g(x,y) = 1 + sin.(π*x)*sin.(π*y)'
f(n,m) = g(linspace(0,1,n+1),linspace(0,1,m+1))
p₁(n,m) = laplace2d(n,m)
p₂(n,m) = elliptic2d(f(n,m))
p₃(ϵ,n,m) = anisotropic2d(ϵ,n,m)

n = m = 512
A = p₃(1e-8,n,m)
b = ones((n-1)*(m-1))
#x = A\b
A = V_cycle(A,(n,m))
x = A\b
@show length(A.resnorm)
@show A.resnorm

#msurf([reshape(A.grids[i].x,A.grids[i].sz.-1) for i in 1:length(A.grids)]...)
