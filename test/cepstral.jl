series = [0.5783834024666892, 0.41692819343103993, 0.28467412058772745, 0.040996614637221374, 0.06804363629643206, 0.9244312331822693, 0.08725427297969246, 0.7808434807429435, 0.17007538615159856, 0.832499149162595]

@testset "cepstral.jl" begin
  @testset "ARMACepstral" begin
    @test (cc(ARMACepstral(1., (1, 1)), series, 5) .≈ [-0.6231197903689889, -0.01564051553647603, -0.005955402230438618, -0.0027120968883413143, -0.0015279908018555185]) |> all
  end

  @testset "ARCepstral" begin
    @test (cc(ARCepstral(1), series, 5) .≈ [0.5686327103717447, 0.16167157965235823, 0.06128783235186792, 0.026137699667288815, 0.011890200803754471]) |> all
    @test (cc(ARCepstral(2), series, 5) .≈ [0.5061405677355508, -0.08265600639423334, -0.06344592753639652, -0.015374742371503975, 0.0017971279975581714]) |> all
  end

  @testset "RealCepstral" begin
    @test (cc(RealCepstral(), series, 5) .≈ [0.03700275546236409, 0.0867638177224842, 0.38225344455973626, -0.03467668552262836, 0.10502390427503827]) |> all
    @test (cc(RealCepstral(), series, 5, normalize=true) .≈ [-2.482417812401184, -3.5121897532599697, -3.2167001264227184, -3.633630256505082, -3.493929666707416]) |> all
  end

  @testset "Errors" begin
    @test_throws MethodError ARCepstral(2.)
    @test_throws MethodError ARMACepstral(2, (1, 1))
    @test_throws MethodError ARMACepstral(2., (1., 1))
    @test_throws MethodError ARMACepstral(2., (1, 1.))
    @test_throws MethodError ARMACepstral(2., (1., 1.))
    @test_throws ArgumentError cc(RealCepstral(), series, 0)
  end
end
