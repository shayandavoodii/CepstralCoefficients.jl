module CepstralClustering

using ARCHModels
using Statistics

include("arima.jl")
include("cepstral.jl")
include("normalizer.jl")

export cc

"""
    cc(tseries::AbstractMatrix, p::Int, n::Int)

Calculate the cepstral coefficients of the AR(p) process for each asset in the series.

# Arguments
- `tseries::AbstractMatrix`: a matrix of time series observations, with each row representing \
an asset and each column representing a time step.
- `p::Int`: the order of the AR(p) process.
- `n::Int`: the number of cepstral coefficients to calculate.

## Keyword Arguments
- `normalize::Bool=false`: whether to normalize the series before fitting the AR(p) process.

# Returns
- `cc_mat::AbstractMatrix`: a matrix of cepstral coefficients, with each row representing an \
asset and each column representing a cepstral coefficient.
"""
function cc(
  tseries::AbstractMatrix,
  p::Int,
  n::Int;
  normalize::Bool=false
)
  tseries = permutedims(tseries)
  series = copy(tseries)
  n_assets = size(series, 2)
  cc_mat = similar(series, n, n_assets)
  normalize && normalizer!(series)
  α = fit_arima(series, p)
  for asset ∈ 1:n_assets
    cc_mat[:, asset] = cepscoef(p, α[:, asset], n)
  end
  return cc_mat
end

end #module
