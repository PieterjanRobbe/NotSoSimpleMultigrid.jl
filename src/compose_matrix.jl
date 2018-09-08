## compose_matrix.jl : methods to compose example matrices

"""
anisotropic2d(ϵ, n, m)

2d anisotropic Laplacian on [0,1]^2 using n x m points and anisotropy number ϵ
"""
anisotropic2d(ϵ, n, m) = stencil2mat(SArray{Tuple{3,3}}([0 -m^2 0 -ϵ*n^2 2(ϵ*n^2+m^2) -ϵ*n^2 0 -m^2 0]), n, m)

"""
anisotropic2d(ϵ, β, n, m)

2d rotated anisotropic Laplacian on [0,1]^2 using n x m points with anisotropy number ϵ and angle β
"""
anisotropic2d(ϵ, β, n, m) = stencil2mat(SArray{Tuple{3,3}}([0 0 0 -n^2*(ϵ*cos(β)^2+sin(β)^2) 2n^2*(ϵ*cos(β)^2+sin(β)^2) -n^2*(ϵ*cos(β)^2+sin(β)^2) 0 0 0]), n, m) + stencil2mat(SArray{Tuple{3,3}}([n*m/2*((ϵ-1)*cos(β)*sin(β)) 0 -n*m/2*((ϵ-1)*cos(β)*sin(β)) 0 0 0 -n*m/2*((ϵ-1)*cos(β)*sin(β)) 0 n*m/2*((ϵ-1)*cos(β)*sin(β))]), n, m) + stencil2mat(SArray{Tuple{3,3}}([0 -m^2*(ϵ*sin(β)^2+cos(β)^2) 0 0 2m^2*(ϵ*sin(β)^2+cos(β)^2) 0 0 -m^2*(ϵ*sin(β)^2+cos(β)^2) 0]), n, m)


"""
anisotropic3d(ϵ₁, ϵ₂, n, m, l)

3d anisotropic Laplacian on [0,1]^3 using n x m x l points and anisotropy numbers ϵ₁ and ϵ₂
"""
anisotropic3d(ϵ₁, ϵ₂, n, m, l) = SArray{Tuple{3,3,3}}([0 0 0 0 -l^2 0 0 0 0 0 -ϵ₂*m^2 0 -ϵ₁*n^2 2(ϵ₁*n^2+ϵ₂*m^2+l^2) -ϵ₁*n^2 0 -ϵ₂*m^2 0 0 0 0 0 -l^2 0 0 0 0], n, m, l)
