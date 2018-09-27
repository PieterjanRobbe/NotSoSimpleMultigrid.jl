# test_problems.jl : definition of some test problems

g(x,y) = Float64[1 + sin(π*xᵢ)*sin(π*yᵢ) for xᵢ in x, yᵢ in y]
g(x,y,z) = Float64[1 + sin(π*xᵢ)*sin(π*yᵢ)*sin(π*zᵢ) for xᵢ in x, yᵢ in y, zᵢ in z]

f(n,m) = g(linspace(0,1,n),linspace(0,1,m))
f(n,m,l) = g(linspace(0,1,n),linspace(0,1,m),linspace(0,1,l))

p₁₂(n) = laplace2d(n,n)
p₁₃(n) = laplace3d(n,n,n)

p₂₂(n) = elliptic2d(f(n+1,n+1))
p₂₃(n) = elliptic3d(f(n+1,n+1,n+1))

p₃₂(n) = anisotropic2d(1e-4,n,n)
p₃₃(n) = anisotropic3d(1e-4,1,n,n,n)

p₄₂(n) = anisotropic2d(1e-4,10*π/180,n,n)
