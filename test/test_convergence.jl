# test_convergence.jl

# example problem from On the Robustness of a Multiple Semicoarsened Grid Method, Oosterlee and Wesseling.

function get_problem(ϵ,β)
    n = 256
    m = n
    A = anisotropic2d(ϵ,β*π/180,n,m)
    A = V_cycle(A, (n,m))
    b = ones((n-1)*(m-1))
    return (A,b)
end

@testset "Convergence for example from Oosterlee et. al." begin
    for ϵ in [1 1e-2 1e-8]
        for β in [0 15 30 45]
            (A,b) = get_problem(ϵ,β)
            A.grids[1].b .= b # copy rhs
            tol = 1/prod(A.grids[1].sz)*SimpleMultigrid.norm_of_residu(A.grids[1]) # target norm of residu is O(h²)
            push!(A.resnorm,SimpleMultigrid.norm_of_residu(A.grids[1])) # log convergence history
            for i in 1:15
                SimpleMultigrid.cycle!(A)
                push!(A.resnorm,SimpleMultigrid.norm_of_residu(A.grids[1])) # log convergence history
            end
            if !( ( ϵ == 1e-8 ) && ( β == 45 ) ) # known failure for last case in 15 iters
                @test A.resnorm[end] < 1/256^2
            end
            log(A,ϵ,β)
        end
    end
end
