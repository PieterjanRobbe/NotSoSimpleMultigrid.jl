# test_problems.jl : definition of some test problems

g(x,y) = Float64[1 + sin(π*xᵢ)*sin(π*yᵢ) for xᵢ in x, yᵢ in y]
g(x,y,z) = Float64[1 + sin(π*xᵢ)*sin(π*yᵢ)*sin(π*zᵢ) for xᵢ in x, yᵢ in y, zᵢ in z]

f(n,m) = g(range(1/n,stop=1-1/n,length=n-1),range(1/m,stop=1-1/m,length=m-1))
f(n,m,l) = g(range(1/n,stop=1-1/n,length=n-1),range(1/m,stop=1-1/m,length=m-1),range(1/l,stop=1-1/l,length=l-1))

p₁₂(n) = laplace2d(n,n)
p₁₃(n) = laplace3d(n,n,n)

p₂₂(n) = elliptic2d(f(n,n))
p₂₃(n) = elliptic3d(f(n,n,n))

p₃₂(n) = anisotropic2d(1e-4,n,n)
p₃₃(n) = anisotropic3d(1e-4,1,n,n,n)

p₄₂(n) = anisotropic2d(1e-4,10*π/180,n,n)
