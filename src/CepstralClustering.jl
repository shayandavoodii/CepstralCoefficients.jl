module CepstralClustering

using ARCHModels
using Statistics
using Clustering
using Distances
using FFTW

include("arima.jl")
include("cepstral.jl")
include("normalizer.jl")
include("clustering.jl")

export cc, clustering
export RealCepstral, ARCepstral

"""
    cc(tseries::AbstractMatrix, p::Int, n::Int)

Calculate the cepstral coefficients of the AR(p) process for each asset in the series.

# Arguments
- `model::Type{<:CepstralCoeffModel}`: A subtype of `CepstralCoeffModel`. Currently, \
  the `ARCepstral` and `RealCepstral` are supported.
- `tseries::AbstractMatrix`: a matrix of time series observations, with each row representing \
  an asset and each column representing a time step.
- `p::Int`: the order of the AR(p) process.
- `n::Int`: the number of cepstral coefficients to calculate.

## Keyword Arguments
- `normalize::Bool=false`: whether to normalize the series before fitting the AR(p) process.

# Returns
- `cc_mat::AbstractMatrix`: a n×m matrix of cepstral coefficients, with each row representing \
a cepstral coefficient and each column representing an asset.
"""
function cc(
  model::Type{<:CepstralCoeffModel},
  tseries::AbstractMatrix,
  p::Integer,
  n::Integer;
  normalize::Bool=false
)
  tseries  = permutedims(tseries)
  n_assets = size(tseries, 2)
  cc_mat   = similar(tseries, n, n_assets)
  for asset ∈ 1:n_assets
    cc_mat[:, asset] = cc(model, tseries[:, asset], p, n, normalize=normalize)
  end
  return cc_mat
end

function cc(
  model::Type{<:CepstralCoeffModel},
  tseries::AbstractVector,
  p::Integer,
  n::Integer;
  normalize::Bool=false
)
  series = tseries
  if normalize
    series = copy(tseries)
    normalizer!(series)
  end
  α = fit_arima(series, p)
  coefs = cepscoef(model, α, p, n)
  return real.(coefs)
end

end #module
