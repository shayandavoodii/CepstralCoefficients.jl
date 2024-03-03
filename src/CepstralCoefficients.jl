module CepstralCoefficients

using ARCHModels
using FFTW
using QuadGK

include("arima.jl")
include("cepstral.jl")
include("tools.jl")
include("clustering.jl")

export cc, cepsclustering
export RealCepstral, ARCepstral, ARMACepstral

"""
    cc(
      method::CepstralCoeffModel,
      tseries::AbstractMatrix,
      n::Integer;
      normalize::Bool=false
    )

    cc(
      method::CepstralCoeffModel,
      tseries::AbstractVector,
      n::Integer;
      normalize::Bool=false
    )

Calculate the cepstral coefficients of the AR(p), or ARMA(p, q) processes or Real cepstral \
coefficients for each asset according to the given time series. See [`ARCepstral`](@ref), \
[`ARMACepstral`](@ref), and [`RealCepstral`](@ref) for more information.

=========================================================================================================================================
# 1st Method
## Arguments
- `method::CepstralCoeffModel`: A subtype of `CepstralCoeffModel`. Currently, \
  the [`ARCepstral`](@ref), [`ARMACepstral`](@ref), and [`RealCepstral`](@ref) methods \
  are supported.
- `tseries::AbstractMatrix`: a matrix of time series observations, with each row representing \
  an asset and each column representing a time step.
- `n::Integer`: the number of cepstral coefficients to calculate.

### Keyword Arguments
- `normalize::Bool=false`: whether to normalize the series before fitting the AR(p) process.

## Returns
- `cc_mat::AbstractMatrix`: a n×m matrix of cepstral coefficients, with each row representing \
a cepstral coefficient and each column representing an asset.

=========================================================================================================================================

# 2nd Method
## Arguments
- `method::CepstralCoeffModel`: A subtype of `CepstralCoeffModel`. Currently, \
  the [`ARCepstral`](@ref), [`ARMACepstral`](@ref), and [`RealCepstral`](@ref) methods \
  are supported.
- `tseries::AbstractVector`: a vector of time series observations.
- `n::Integer`: the number of cepstral coefficients to calculate.

### Keyword Arguments
- `normalize::Bool=false`: whether to normalize the series before fitting the AR(p) process.

## Returns
- `cc_vec::AbstractVector`: a vector of cepstral coefficients.
"""
function cc(
  method::CepstralCoeffModel,
  tseries::AbstractMatrix,
  n::Integer;
  normalize::Bool=false
)
  tseries  = permutedims(tseries)
  n_assets = size(tseries, 2)
  cc_mat   = similar(tseries, n, n_assets)
  for asset ∈ 1:n_assets
    cc_mat[:, asset] = cc(method, tseries[:, asset], n, normalize=normalize)
  end
  return cc_mat
end

function cc(
  method::CepstralCoeffModel,
  tseries::AbstractVector,
  n::Integer;
  normalize::Bool=false
)
  series = tseries
  if normalize
    series = copy(tseries)
    normalizer!(series)
  end
  coefs = cepscoef(method, series, n)
  return coefs
end

end #module
