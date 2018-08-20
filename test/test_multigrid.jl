# test_multigrid.jl : test multigrid implementation

g(x,y) = 1 + sin.(π*x)*sin.(π*y)'
f(n,m) = g(linspace(0,1,n+1),linspace(0,1,m+1))
p₁(n,m) = laplace2d(n,m)
p₂(n,m) = elliptic2d(f(n,m))
p₃(n,m) = anisotropic2d(1e-4,n,m)
p₄(n,m) = anisotropic2d(1e-4*pi/180,n,m)

@testset "Default Multigrid methods for 2D anisotropic" begin
    for (problem,problem_name) in zip([p₁,p₂,p₃,p₄],["     LAPLACE 2D      ","     ELLIPTIC 2D     ","   ANISOTROPIC 2D    ","ROTATED ANISTROPIC 2D"])
        println("*************************")
        println("* $(problem_name) *")
        println("*************************")
        for (method,method_name) in zip([V_cycle, W_cycle, F_cycle],["V-CYCLE", "W-CYCLE", "  FMG  "])
            println("+-----------------------+")
            println("|         $(method_name)       |")
            println("+-----------------------+")
            for n in 2.^(2:9)
                print("n = $(n)...")
                A = method(problem(n,n),(n,n))
                b = ones((n-1)*(n-1))
                x = A\b
                @test A.resnorm[end] < 1/n^2
                println("done")
            end
        end
    end
end
