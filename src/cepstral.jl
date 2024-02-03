abstract type CepstralCoeffModel end

struct ARCepstral <: CepstralCoeffModel
  p::Integer
end

struct RealCepstral <: CepstralCoeffModel end

struct ARMACepstral{S<:AbstractFloat, T<:Integer} <: CepstralCoeffModel
  σ²::S
  order::Tuple{T, T}
end

function cepscoef(m::ARCepstral, tseries::AbstractVector, n::Integer)
  n>0 || ArgumentError("n must be positive") |> throw
  α = fit_arima(tseries, m.p)
  m.p == length(α) || ArgumentError("Length of coefficients must be equal to p") |> throw
  c = similar(α, n)
  for n_ ∈ 1:n
    if n_==1
      c[1] = -α[1]
    elseif 1<n_≤m.p
      secondterm = 0.
      for m ∈ 1:n_-1
        secondterm += (1-m/n_)*α[m]*c[n_-m]
      end
      c[n_] = -α[n_] - secondterm
    else
      res = 0.
      for m ∈ 1:m.p
        res += (1-m/n_)*α[m]*c[n_-m]
      end
      c[n_] = -res
    end
  end
  return c
end

function cepscoef(::RealCepstral, tseries::AbstractVector, n::Integer)
  res = tseries |> fft .|> abs .|> log |> ifft
  return real.(res[1:n])
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
function cepscoef(
  m::ARMACepstral,
  tseries::AbstractVector,
  n::Integer
)
  p, q = first(m.order), last(m.order)
  α = fit_arima(tseries, p, q)
  ψarr = similar(α, n)
  for k ∈ 1:n
    ψarr[k] = quadgk(ω->ψfunc(α, m.σ², ω, p, q, k-1), 0, 1) |> first
  end
  return ψarr
end

function ψfunc(
  coefs::AbstractVector,
  σ²::AbstractFloat,
  freq::AbstractFloat,
  p::Integer,
  q::Integer,
  k::Integer
)
  p+q == length(coefs) || ArgumentError("Length of coefficients must be equal to p + q") |> throw
  a = σ²/2π
  numerator_ = 0.
  denominator_ = 0.
  for h ∈ 1:p
    x = h * freq
    numerator_ += coefs[h] * cis(x)
  end
  for h ∈ 1:q
    x = h * freq
    denominator_ += coefs[p + h] * cis(x)
  end
  numerator_ = 1 - numerator_
  denominator_ = 1 - denominator_
  λxw = (a*abs(numerator_/denominator_) |> log) * cos(2π*k*freq)
  return λxw
end
