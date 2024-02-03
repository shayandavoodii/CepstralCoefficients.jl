module CepstralClustering

using ARCHModels
using Statistics
using Clustering
using Distances
using FFTW
using QuadGK

include("arima.jl")
include("cepstral.jl")
include("normalizer.jl")
include("clustering.jl")

export cc, clustering
export RealCepstral, ARCepstral, ARMACepstral

"""
    cc(
      model::CepstralCoeffModel,
      tseries::AbstractMatrix,
      n::Integer;
      normalize::Bool=false
    )

Calculate the cepstral coefficients of the AR(p), or ARMA(p, q) processes or Real cepstral \
coefficients for each asset according to the given time series. See [`ARCepstral`](@ref), \
[`ARMACepstral`](@ref), and [`RealCepstral`](@ref) for more information.

# Arguments
- `model::Type{<:CepstralCoeffModel}`: A subtype of `CepstralCoeffModel`. Currently, \
  the [`ARCepstral`](@ref), [`ARMACepstral`](@ref), and [`RealCepstral`](@ref) are supported.
- `tseries::AbstractMatrix`: a matrix of time series observations, with each row representing \
  an asset and each column representing a time step.
- `n::Integer`: the number of cepstral coefficients to calculate.

## Keyword Arguments
- `normalize::Bool=false`: whether to normalize the series before fitting the AR(p) process.

# Returns
- `cc_mat::AbstractMatrix`: a n×m matrix of cepstral coefficients, with each row representing \
a cepstral coefficient and each column representing an asset.
"""
function cc(
  model::CepstralCoeffModel,
  tseries::AbstractMatrix,
  n::Integer;
  normalize::Bool=false
)
  tseries  = permutedims(tseries)
  n_assets = size(tseries, 2)
  cc_mat   = similar(tseries, n, n_assets)
  for asset ∈ 1:n_assets
    cc_mat[:, asset] = cc(model, tseries[:, asset], n, normalize=normalize)
  end
  return cc_mat
end

function cc(
  model::CepstralCoeffModel,
  tseries::AbstractVector,
  n::Integer;
  normalize::Bool=false
)
  series = tseries
  if normalize
    series = copy(tseries)
    normalizer!(series)
  end
  coefs = cepscoef(model, series, n)
  return coefs
end

end #module
