abstract type CepstralCoeffModel end

struct ARCepstral{T<:Integer} <: CepstralCoeffModel
  p::T
end

struct RealCepstral <: CepstralCoeffModel end

struct ARMACepstral{S<:AbstractFloat, T<:Integer} <: CepstralCoeffModel
  σ²::S
  order::NTuple{2, T}
end

function cepscoef(method::ARCepstral, tseries::AbstractVector, n::Integer)
  α = fit_arima(tseries, method.p)
  c = similar(α, n)
  for n_ ∈ 1:n
    if n_==1
      c[1] = -α[1]
    elseif 1<n_≤method.p
      secondterm = 0.
      for m ∈ 1:n_-1
        secondterm += (1-m/n_)*α[m]*c[n_-m]
      end
      c[n_] = -α[n_] - secondterm
    else
      res = 0.
      for m ∈ 1:method.p
        res += (1-m/n_)*α[m]*c[n_-m]
      end
      c[n_] = -res
    end
  end
  return c
end

function cepscoef(::RealCepstral, tseries::AbstractVector, n::Integer)
  res = tseries |> fft .|> abs .|> log |> ifft
  return res[1:n] .|> real
end

function cepscoef(
  m::ARMACepstral,
  tseries::AbstractVector,
  n::Integer
)
  p, q = m.order
  α    = fit_arima(tseries, p, q)
  ψarr = similar(α, n)
  for k ∈ 1:n
    ψarr[k] = quadgk(ω->ψfunc(α, m.σ², ω, p, q, k-1), 0, 1) |> first
  end
  return ψarr
end

"""
    ψfunc(
      coefs::AbstractArray,
      σ²::AbstractFloat,
      freq::AbstractFloat,
      order::NTuple{3, Int},
      n::Int
    )

Compute the first n cepstral coefficients for each asset and each time period

# Arguments
- `α`: Vector of coefficients of the ARMA(p, q) model.
- `σ²`: Variance of the white noise process.
- `order`: Order of the ARMA(p, q) model.
- `n`: Number of cepstral coefficients to compute.

# Returns
- `ψarr::AbstractArray`: Array of cepstral coefficients. First index: coefficient, \
second index: asset, third index: time window
"""
function ψfunc(
  coefs::AbstractVector,
  σ²::S,
  freq::S,
  p::T,
  q::T,
  k::T
) where {S<:AbstractFloat, T<:Integer}
  a            = σ²/2π
  numerator_   = 0.
  denominator_ = 0.
  for h ∈ 1:p
    x = h * freq
    numerator_ += coefs[h] * cis(x)
  end
  for h ∈ 1:q
    x = h * freq
    denominator_ += coefs[p + h] * cis(x)
  end
  numerator_   = 1 - numerator_
  denominator_ = 1 - denominator_
  ψ = (a*abs(numerator_/denominator_) |> log) * cos(2π*k*freq)
  return ψ
end
