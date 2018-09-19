# test_multigrid.jl : test multigrid implementation

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

p₄₂(n) = anisotropic2d(1e-4,10*pi/180,n,n)

for (d,d_name) in zip(2:3,["₂","₃"])
    @testset "Default Multigrid methods for $(d)D anisotropic problems" begin
        for (problem,problem_name) in zip(["p₁","p₂","p₃","p₄"],["     LAPLACE $(d)D      ","     ELLIPTIC $(d)D     ","   ANISOTROPIC $(d)D    ","ROTATED ANISTROPIC $(d)D"])
            if !( d == 3 && problem == "p₄" )
                println("*************************")
                println("* $(problem_name) *")
                println("*************************")
                for (method,method_name) in zip([V_cycle, W_cycle, F_cycle],["V-CYCLE", "W-CYCLE", "  FMG  "])
                    println("+-----------------------+")
                    println("|         $(method_name)       |")
                    println("+-----------------------+")
                    for n in 2.^(1:(12-2d))
                        print("n = $(n)...")
                        A = method(eval(Symbol(problem,d_name))(n),ntuple(i->n,d))
                        b = ones(prod([n-1 for i in 1:d]))
                        x = A\b 
                        @test A.resnorm[end] < 1/n^d
                        println("done ($(length(A.resnorm-1)) iterations)")
                    end 
                end 
            end 
        end 
    end 
end
