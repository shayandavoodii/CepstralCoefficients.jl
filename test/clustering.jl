using Clustering

vals = [
  0.526449  0.253821  0.536858  0.15444   0.39956    0.278192  0.145485  0.685143   0.495732  0.596486
  0.913692  0.205756  0.403905  0.529887  0.104882   0.180711  0.651375  0.208858   0.307472  0.277686
  0.537414  0.870541  0.163864  0.374202  0.0393039  0.68875   0.828697  0.0172519  0.668273  0.765799
]

@testset "clustering.jl" begin
  cepsclustering(vals, 3)

  @testset "With invalid arguments" begin
    @test_throws ArgumentError cepsclustering(vals, 4)
  end
end
