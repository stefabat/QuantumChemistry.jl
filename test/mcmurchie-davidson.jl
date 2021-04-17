
@testset "Hermite expansion coefficients" begin

    tol = 1e-15

    E = hermite_expansion(3, 3, 1.55, 4.0, 0.2, 0.8)
    # note that Eₜⁱʲ = E[t+1,i+1,j+1], e.g. E₀¹² = E[1,2,3]
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

@testset "Hermite Coulomb integrals" begin

    tol = 1e-15

    R = hermite_integral(5, 5, 1, 4.0, 0.8, 0.4, 0.2)
    # note that Rⁿₜᵤᵥ = R[n+1,t+1,u+1,v+1], e.g. R⁰₂₀₁ = R[1,3,1,2]
    @test R[1,1,1,1] ≈  0.47886707534919914 rtol = tol atol = tol
    @test R[1,2,1,1] ≈ -0.42298268228996244 rtol = tol atol = tol
    @test R[1,3,1,1] ≈  0.46807392296951067 rtol = tol atol = tol
    @test R[1,4,1,1] ≈  0.34633694124289693 rtol = tol atol = tol
    @test R[1,5,1,1] ≈ -6.82592015565251800 rtol = tol atol = tol
    @test R[1,6,1,1] ≈ 30.76951320573933000 rtol = tol atol = tol

    @test R[1,1,3,1] ≈ -0.27952778390446210 rtol = tol atol = tol
    @test R[1,2,3,1] ≈  0.39808494650821280 rtol = tol atol = tol
    @test R[1,3,3,1] ≈ -0.25720677706409134 rtol = tol atol = tol
    @test R[1,4,3,1] ≈ -1.40091854240456470 rtol = tol atol = tol
    @test R[1,5,3,1] ≈  5.96994606350125800 rtol = tol atol = tol
    @test R[1,6,3,1] ≈  5.35853864250552200 rtol = tol atol = tol

    @test R[1,1,1,1] ≈   0.4788670753491991 rtol = tol atol = tol
    @test R[1,1,2,1] ≈  -0.2114913411449812 rtol = tol atol = tol
    @test R[1,1,3,1] ≈  -0.2795277839044621 rtol = tol atol = tol
    @test R[1,1,4,1] ≈   1.4450453180440610 rtol = tol atol = tol
    @test R[1,1,5,1] ≈  -0.8156794363483952 rtol = tol atol = tol
    @test R[1,1,6,1] ≈ -19.9655797366111650 rtol = tol atol = tol

    @test R[1,4,1,2] ≈    0.3050844126528191 rtol = tol atol = tol
    @test R[1,4,2,2] ≈   -1.4631703028343300 rtol = tol atol = tol
    @test R[1,4,3,2] ≈    1.0367547739123495 rtol = tol atol = tol
    @test R[1,4,4,2] ≈   21.2827655740582670 rtol = tol atol = tol
    @test R[1,4,5,2] ≈  -80.3741443438472000 rtol = tol atol = tol
    @test R[1,4,6,2] ≈ -410.4436752870131500 rtol = tol atol = tol

end
