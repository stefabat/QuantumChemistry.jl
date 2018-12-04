
@testset "integrals" begin
    Ga = GTO(zeros(3), 5101.0250, (0,0,0))
    Gb = GTO(zeros(3), 510.10250, (0,0,0))

    @test overlap(Ga,Ga) == overlap(Ga)
    @test kinetic(Ga,Ga)/sqrt(overlap(Ga)*overlap(Ga)) = ...
    @test kinetic(Ga,Gb)/sqrt(overlap(Ga)*overlap(Gb))
    @test kinetic(Gb,Gb)/sqrt(overlap(Gb)*overlap(Gb))
end
