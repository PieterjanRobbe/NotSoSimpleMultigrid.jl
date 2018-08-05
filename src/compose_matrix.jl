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
