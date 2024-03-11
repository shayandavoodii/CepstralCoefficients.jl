using CepstralCoefficients
using Test

@testset "CepstralCoefficients.jl" begin
  @info "Testing cepstral.jl unit tests"
  include("cepstral.jl")
  @info "Testing clustering.jl unit tests"
  include("clustering.jl")
end
