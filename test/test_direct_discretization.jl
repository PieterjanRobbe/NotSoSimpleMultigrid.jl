# test_direct_discretization.jl

# direct discretization for elliptic and anisotropic problem in 2D/3D
function get_problem_2d_elliptic(n)
    A = NotSoSimpleMultigrid.MultigridMethod((n,m)->elliptic2d(f(n,m)), (n, n), V(3, 2))
    b = fill(1,(n-1)*(n-1))
    return (A,b)
end

@testset "Direct discretization of elliptic PDE, 2D" begin
    for n in [16 32 64 128 256 512]
        (A, b) = get_problem_2d_elliptic(n)
        A.grids[1].b .= b # copy rhs
        for item in Base.Iterators.take(A,15) end # iterate
        @test A.resnorm[end] < 1/n^2
        log(A,2)
    end
end

function get_problem_2d_anisotropic(n)
    A = NotSoSimpleMultigrid.MultigridMethod((n,m)->anisotropic2d(1e-8,n,m), (n, n), V(3, 2))
    b = fill(1,(n-1)*(n-1))
    return (A,b)
end

@testset "Direct discretization of anisotropic PDE, 2D" begin
    for n in [16 32 64 128 256 512]
        (A, b) = get_problem_2d_anisotropic(n)
        A.grids[1].b .= b # copy rhs
        for item in Base.Iterators.take(A,15) end # iterate
        @test A.resnorm[end] < 1/n^2
        log(A,2)
    end
end

function get_problem_2d_rotated_anisotropic(n)
    A = NotSoSimpleMultigrid.MultigridMethod((n,m)->anisotropic2d(1e-4,10*π/180,n,m), (n, n), V(3, 2))
    b = fill(1,(n-1)*(n-1))
    return (A,b)
end

@testset "Direct discretization of rotated anisotropic PDE, 2D" begin
    for n in [16 32 64 128 256 512]
        (A, b) = get_problem_2d_anisotropic(n)
        A.grids[1].b .= b # copy rhs
        for item in Base.Iterators.take(A,15) end # iterate
        @test A.resnorm[end] < 1/n^2
        log(A,2)
    end
end
function get_problem_3d_elliptic(n)
    A = NotSoSimpleMultigrid.MultigridMethod((n,m,k)->elliptic3d(f(n,m,k)), (n, n, n), V(3, 2))
    b = fill(1,(n-1)*(n-1)*(n-1))
    return (A,b)
end

@testset "Direct discretization of elliptic PDE, 3D" begin
    for n in [4 8 16 32 64]
        (A, b) = get_problem_3d_elliptic(n)
        A.grids[1].b .= b # copy rhs
        for item in Base.Iterators.take(A,15) end # iterate
        @test A.resnorm[end] < 1/n^3
        log(A,3)
    end
end

function get_problem_3d_anisotropic(n)
    A = NotSoSimpleMultigrid.MultigridMethod((n,m,k)->anisotropic3d(1e-8,1,n,m,k), (n, n, n), V(3, 2))
    b = fill(1,(n-1)*(n-1)*(n-1))
    return (A,b)
end

@testset "Direct discretization of anisotropic PDE, 3D" begin
    for n in [4 8 16 32 64]
        (A, b) = get_problem_3d_anisotropic(n)
        A.grids[1].b .= b # copy rhs
        for item in Base.Iterators.take(A,15) end # iterate
        @test A.resnorm[end] < 1/n^3
        log(A,3)
    end
end
