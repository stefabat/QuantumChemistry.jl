
@testset "McMurchie-Davidson" begin

    tol = 1e-15

    E = hermite_expansion(3, 3, 1.55, 0.2, 0.8, 4.0)
    # note that E = E[t+1,i+1,j+1], e.g. E₀¹² = E[1,2,3]
    @test E[1,1,1] ≈ 1.55       rtol = tol atol = tol
    @test E[1,1,2] ≈ 1.24       rtol = tol atol = tol
    @test E[1,1,3] ≈ 1.18575    rtol = tol atol = tol
    @test E[1,1,4] ≈ 1.2586     rtol = tol atol = tol

    @test E[1,2,1] ≈ 0.31       rtol = tol atol = tol
    @test E[1,2,2] ≈ 0.44175    rtol = tol atol = tol
    @test E[1,2,3] ≈ 0.54715    rtol = tol atol = tol
    @test E[1,2,4] ≈ 0.69637625 rtol = tol atol = tol

    @test E[1,3,1] ≈ 0.25575    rtol = tol atol = tol
    @test E[1,3,2] ≈ 0.2821     rtol = tol atol = tol
    @test E[1,3,3] ≈ 0.36808625 rtol = tol atol = tol
    @test E[1,3,4] ≈ 0.5017815  rtol = tol atol = tol

    @test E[1,4,1] ≈ 0.12865       rtol = tol atol = tol
    @test E[1,4,2] ≈ 0.19882625    rtol = tol atol = tol
    @test E[1,4,3] ≈ 0.28092975    rtol = tol atol = tol
    @test E[1,4,4] ≈ 0.41248270625 rtol = tol atol = tol

end
