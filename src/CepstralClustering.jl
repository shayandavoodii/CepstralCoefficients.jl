module CepstralClustering

using ARCHModels

include("arima.jl")
include("cepstral.jl")

export cc

"""
    cc(series::AbstractMatrix, p::Int, n::Int)

Calculate the cepstral coefficients of the AR(p) process for each asset in the series.

# Arguments
- `series::AbstractMatrix`: a matrix of time series observations, with each row representing \
an asset and each column representing a time step.
- `p::Int`: the order of the AR(p) process.
- `n::Int`: the number of cepstral coefficients to calculate.

# Returns
- `cc_vec::AbstractMatrix`: a matrix of cepstral coefficients, with each row representing an \
asset and each column representing a cepstral coefficient.
"""
function cc(
  series::AbstractMatrix,
  p::Int,
  n::Int
)
  n_assets = size(series, 1)
  cc_vec = similar(series, n_assets, n)
  α = fit_arima(series, p)
  for asset ∈ 1:n_assets
    cc_vec[asset, :] = cepscoef(p, α[asset, :], n)
  end
  return cc_vec
end

end #module
