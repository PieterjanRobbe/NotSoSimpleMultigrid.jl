## compose_matrix.jl : methods to compose example matrices

"""
anisotropic2d(ϵ,n,m)

2d anisotropic Laplacian on [0,1]^2 using n x m points and anisotropy number ϵ
"""
function anisotropic2d(ϵ,n,m)
    stencil = @SMatrix [0 -m^2 0; -ϵ*n^2 2(ϵ*n^2+m^2) -ϵ*n^2; 0 -m^2 0]
    Is = Int64[]
    Js = Int64[]
    Vs = Float64[]
    for i in 1:n-1, j in 1:m-1
        (is,js,vs) = stencil2mat(stencil,n,m,i,j)
        push!(Is,is...)
        push!(Js,js...)
        push!(Vs,vs...)
    end
    return sparse(Is,Js,Vs)
end

"""
anisotropic2d(ϵ,β,n,m)

2d rotated anisotropic Laplacian on [0,1]^2 using n x m points with anisotropy number ϵ and angle β
"""
function anisotropic2d(ϵ,β,n,m)
    stencils = SMatrix{3,3,Float64,9}[]
    push!(stencils, @SMatrix [0 0 0; -1 2 -1; 0 0 0])
    stencils[1] *= n^2*(ϵ*cos(β)^2+sin(β)^2)
    push!(stencils, @SMatrix [1 0 -1; 0 0 0; -1 0 1])
    stencils[2] *= n*m/2*((ϵ-1)*cos(β)*sin(β))
    push!(stencils, @SMatrix [0 -1 0; 0 2 0; 0 -1 0])
    stencils[3] *= m^2*(ϵ*sin(β)^2+cos(β)^2)
    Is = Int64[]
    Js = Int64[]
    Vs = Float64[]
    for i in 1:n-1, j in 1:m-1
        for k in 1:length(stencils)
            (is,js,vs) = stencil2mat(stencils[k],n,m,i,j)
            push!(Is,is...)
            push!(Js,js...)
            push!(Vs,vs...)
        end
    end
    return sparse(Is,Js,Vs)
end

