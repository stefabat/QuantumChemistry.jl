
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

@testset "Overlap integrals" begin

    α = 3.5; la = 0; Ax = 0.0; Ay = 0.0; Az = 0.0; # s-type GTO shell Ga
    β = 1.5; lb = 1; Bx = 0.0; By = 0.0; Bz = 1.0; # p-type GTO shell Gb
    γ = 0.5; lc = 2; Cx = 0.4; Cy = 0.3; Cz = 0.2; # d-type GTO shell Gc
    δ = 0.5; ld = 3; Dx = 0.1; Dy = 0.3; Dz = 0.5; # f-type GTO shell Gd

    # self-overlaps, needed for normalization
    Sa = overlap(α, la, Ax, Ay, Az, α, la, Ax, Ay, Az)
    @test Sa ≈ [0.30066145098071745]
    # Na = 1.0/sqrt(Sa[1,1])

    Sb = overlap(β, lb, Bx, By, Bz, β, lb, Bx, By, Bz)
    @test Sb ≈
    [0.17860420377260638 0.0 0.0;
     0.0 0.17860420377260638 0.0;
     0.0 0.0 0.17860420377260638]
    # Nb = 1.0/sqrt(Sb[1,1])

    Sc = overlap(γ, lc, Cx, Cy, Cz, γ, lc, Cx, Cy, Cz)
    @test Sc ≈
    [4.176245997623781  0.0 0.0 1.3920819992079267 0.0 1.3920819992079267;
     0.0 1.3920819992079267 0.0 0.0 0.0 0.0;
     0.0 0.0 1.3920819992079267 0.0 0.0 0.0;
     1.3920819992079267 0.0 0.0 4.176245997623781  0.0 1.3920819992079267;
     0.0 0.0 0.0 0.0 1.3920819992079267 0.0;
     1.3920819992079267 0.0 0.0 1.3920819992079267 0.0 4.176245997623781 ]

    Sd = overlap(δ, ld, Dx, Dy, Dz, δ, ld, Dx, Dy, Dz)
    @test Sd[1:end,1:2] ≈
    [10.44061499405945  0.0;
     0.0                2.0881229988118903;
     0.0                0.0;
     2.0881229988118903 0.0;
     0.0                0.0;
     2.0881229988118903 0.0;
     0.0                2.0881229988118903;
     0.0                0.0;
     0.0                0.6960409996039634;
     0.0                0.0]

    Sab = overlap(α, la, Ax, Ay, Az, β, lb, Bx, By, Bz)
    @test Sab ≈ [0.0 0.0 -0.12199966455329908]

    Sac = overlap(α, la, Ax, Ay, Az, γ, lc, Cx, Cy, Cz)
    @test Sac ≈ [0.15174308630054437;
                 0.05632887294489904;
                 0.03755258196326603;
                 0.11888457708268658;
                 0.02816443647244952;
                 0.09541421335564533]'

    Sad = overlap(α, la, Ax, Ay, Az, δ, ld, Dx, Dy, Dz)
    @test Sad ≈
    [-0.01999633222652011;
     -0.02079651211961892;
     -0.03466085353269818;
     -0.010132890278934855;
     -0.006001349198241031;
     -0.016534329423725284;
     -0.06959115539674599;
     -0.050664451394674265;
     -0.04960298827117586;
     -0.1479924547185288]'

    Sbc = overlap(β, lb, Bx, By, Bz, γ, lc, Cx, Cy, Cz)
    @test Sbc ≈
    [-0.163565618962189     0.03595623520289499  -0.09588329387438661;
     -0.06979739774679614  -0.098615262725587    -0.0190356539309444 ;
      0.18612639399145642  -0.019035653930944408 -0.05499188913383941;
      0.042389581207334534 -0.1268382635190358   -0.08477916241466905;
     -0.019035653930944404  0.19723052545117398  -0.04124391685037954;
      0.08601295479908215   0.0645097160993116    0.2509886222006004]'

    Sbd = overlap(β, lb, Bx, By, Bz, δ, ld, Dx, Dy, Dz)
    @test Sbd ≈
    [0.32856546266268327   -0.007338503569025084   0.012230839281708475;
     0.0120851617915045     0.10288944848633026    0.012412936144463434;
    -0.020141936319174165   0.012412936144463435   0.08964898326556925;
     0.1287876689670364     0.011648129320892585   0.004866032832507672;
    -0.03614622725686051   -0.011320354967933648   0.005918148039536357;
     0.1673436447076876    -0.0037936846407284337 -0.01795677396611459;
    -0.007775536039637      0.3659560184817027     0.038877680198184994;
     0.0048660328325076715 -0.05824064660446291    0.10543071137099953;
    -0.0037936846407284337  0.15722715233241175   -0.05387032189834376;
     0.014416001634768048   0.04324800490430414    0.4337446105899509]'

    Scd = overlap(γ, lc, Cx, Cy, Cz, δ, ld, Dx, Dy, Dz)
    @test Scd ≈
    [0.5813101742955326 0.0 0.5984678599200383 0.6078551800653832 0.0 0.6352086631683255;
     0.0 0.19064095138339596 0.0 0.0 0.20860711434099352 0.0;
    -0.5813101742955329 -0.0 0.18206210857114313 -0.2086071143409936 0.0 0.1992197941956487;
    -0.1906409513833961 0.0 0.19064095138339596 0.5988720985865844 0.0 0.20860711434099355;
     0.0 -0.19064095138339604 0.0 -0.0 0.190640951383396 0.0;
    -0.19921979419564892 0.0 -0.1820621085711433 0.20860711434099355 -0.0 0.5813101742955328;
     0.0 -0.5988720985865844 -0.0 0.0 0.5988720985865842 0.0;
    -0.2086071143409936 0.0 -0.19064095138339604 -0.5988720985865844 0.0 0.1906409513833959;
     0.0 -0.2086071143409936 0.0 0.0 -0.19064095138339612 0.0;
    -0.6352086631683255 0.0 -0.5984678599200385 -0.6078551800653833 0.0 -0.5813101742955333]'

end

@testset "Kinetic integrals" begin

    α = 3.5; la = 0; Ax = 0.0; Ay = 0.0; Az = 0.0; # s-type GTO shell Ga
    β = 1.5; lb = 1; Bx = 0.0; By = 0.0; Bz = 1.0; # p-type GTO shell Gb
    γ = 0.5; lc = 2; Cx = 0.4; Cy = 0.3; Cz = 0.2; # d-type GTO shell Gc
    δ = 0.5; ld = 3; Dx = 0.1; Dy = 0.3; Dz = 0.5; # f-type GTO shell Gd

    # on-site kinetic
    Ta = kinetic(α, la, Ax, Ay, Az, α, la, Ax, Ay, Az)
    @test Ta ≈ [1.578472617648767]

    Tb = kinetic(β, lb, Bx, By, Bz, β, lb, Bx, By, Bz)
    @test Tb ≈ [0.6697657641472743 0.0 0.0;
                0.0 0.6697657641472743 0.0;
                0.0 0.0 0.6697657641472743]

    Tc = kinetic(γ, lc, Cx, Cy, Cz, γ, lc, Cx, Cy, Cz)
    @test Tc ≈
    [4.524266497425761 0.0 0.0 -0.34802049980198163 0.0 -0.34802049980198163;
     0.0 2.4361434986138715 0.0 0.0 0.0 0.0;
     0.0 0.0 2.4361434986138715 0.0 0.0 0.0;
    -0.34802049980198163 0.0 0.0 4.524266497425761 0.0 -0.34802049980198163;
     0.0 0.0 0.0 0.0 2.4361434986138715 0.0;
    -0.34802049980198163 0.0 0.0 -0.34802049980198163 0.0 4.524266497425761]

    Td = kinetic(δ, ld, Dx, Dy, Dz, δ, ld, Dx, Dy, Dz)
    @test Td ≈
    [10.96264574376242 0.0 0.0 0.5220307497029725 0.0 0.5220307497029725 0.0 0.0 0.0 0.0;
     0.0 3.3061947481188256 0.0 0.0 0.0 0.0 0.5220307497029725 0.0 0.17401024990099082 0.0;
     0.0 0.0 3.3061947481188256 0.0 0.0 0.0 0.0 0.17401024990099082 0.0 0.5220307497029725;
     0.5220307497029725 0.0 0.0 3.3061947481188256 0.0 0.17401024990099082 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 1.5660922491089173 0.0 0.0 0.0 0.0 0.0;
     0.5220307497029725 0.0 0.0 0.17401024990099082 0.0 3.3061947481188256 0.0 0.0 0.0 0.0;
     0.0 0.5220307497029725 0.0 0.0 0.0 0.0 10.96264574376242 0.0 0.5220307497029725 0.0;
     0.0 0.0 0.17401024990099082 0.0 0.0 0.0 0.0 3.3061947481188256 0.0 0.5220307497029725;
     0.0 0.17401024990099082 0.0 0.0 0.0 0.0 0.5220307497029725 0.0 3.3061947481188256 0.0;
     0.0 0.0 0.5220307497029725 0.0 0.0 0.0 0.0 0.5220307497029725 0.0 10.96264574376242]

    # hopping kinetic
    Tab = kinetic(α, la, Ax, Ay, Az, β, lb, Bx, By, Bz)
    @test Tab ≈ [0.0 0.0 -0.3714889785647957]

    Tac = kinetic(α, la, Ax, Ay, Az, γ, lc, Cx, Cy, Cz)
    @test Tac ≈
    [-0.15565629046501359;
      0.16625378835822974;
      0.11083585890548647;
     -0.25263766700731427;
      0.08312689417911485;
     -0.3219100788232432]'

    Tad = kinetic(α, la, Ax, Ay, Az, δ, ld, Dx, Dy, Dz)
    @test Tad ≈
    [0.07826404588850908;
     0.07522054916229694;
     0.12536758193716152;
     0.01289952948258379;
    -0.02282622544659097;
    -0.011448444327113257;
     0.19827017695098154;
     0.06449764741291891;
    -0.03434533298133976;
     0.20871042586981742]'

    Tbc = kinetic(β, lb, Bx, By, Bz, γ, lc, Cx, Cy, Cz)
    @test Tbc ≈
    [-0.4846530367634919      0.013207967795210496 -0.035221247453894614;
     -0.15860863773774078    -0.22882665403558922  -0.05948046989234781;
      0.4229563673006421     -0.059480469892347836 -0.09251724386562542;
      0.00026215334167921167 -0.37650113036156996  -0.000524306683358396;
     -0.059480469892347836    0.4576533080711784   -0.06938793289921905;
      0.13657156351164304     0.10242867263373225   0.7313841939575928]'

    Tbd = kinetic(β, lb, Bx, By, Bz, δ, ld, Dx, Dy, Dz)
    @test Tbd ≈
    [-0.1301435940319088   -0.002234413863907837  0.0037240231065130648;
      0.04684069064855518  -0.04580727830680528   0.004457247505449855;
     -0.07806781774759197   0.004457247505449853 -0.05056167564595174;
      0.01962619437647875   0.04508095209110688   0.004418646764230437;
     -0.11823051377570767  -0.03647727366282208   0.018366887082796664;
      0.14573874240390022  -0.006170665173434843 -0.06926912496035052;
     -0.0039941524213561305 0.049547313923906315  0.019970762106780657;
      0.004418646764230437 -0.22540476045553443  -0.0015833100918272737;
     -0.006170665173434843  0.1292836352747407   -0.20780737488105155;
      0.012522715893754515  0.037568147681263556  0.38077331291636335]'

    Tcd = kinetic(γ, lc, Cx, Cy, Cz, δ, ld, Dx, Dy, Dz)
    @test Tcd ≈
    [1.2503572445158508 0.0 0.7607600651551513 0.15176167568307286 0.0 0.1996211757432247;
     0.0 0.6197852113293099 0.0 0.0 0.2695398550717677 0.0;
    -0.8690753417490591 -0.0 0.5876054554133645 -0.06093274073077415 0.0 0.4613240158618593;
    -0.42914425994591393 0.0 0.22952022708371908 0.9347395338772271 0.0 0.060932740730774246;
     0.0 -0.4201611784671152 0.0 -0.0 0.4201611784671152 0.0;
    -0.4613240158618594 0.0 -0.5876054554133645 0.0609327407307742 -0.0 0.8690753417490591;
     0.0 -0.7351155010150324 -0.0 0.0 0.7351155010150321 0.0;
    -0.060932740730774204 0.0 -0.22952022708371925 -0.9347395338772271 0.0 0.4291442599459139;
     0.0 -0.26953985507176775 0.0 0.0 -0.61978521132931 0.0;
    -0.1996211757432245 0.0 -0.7607600651551515 -0.1517616756830727 0.0 -1.2503572445158515]'


end