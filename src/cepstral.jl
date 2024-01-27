abstract type CepstralCoeffModel end

struct ARCepstral <: CepstralCoeffModel end

struct RealCepstral <: CepstralCoeffModel end

function cepscoef(::Type{ARCepstral}, α::AbstractVector, p::Integer, n::Int)
  n>0 || ArgumentError("n must be positive") |> throw
  c = similar(α, n)
  for n_ ∈ 1:n
    if n_==1
      c[1] = -α[1]
    elseif 1<n_≤p
      secondterm = 0.
      for m ∈ 1:n_-1
        secondterm += (1-m/n_)*α[m]*c[n_-m]
      end
      c[n_] = -α[n_] - secondterm
    else
      res = 0.
      for m ∈ 1:p
        res += (1-m/n_)*α[m]*c[n_-m]
      end
      c[n_] = -res
    end
  end
  return c
end

function cepscoef(::Type{RealCepstral}, series::AbstractVector, p::Integer, n::Int)
  p≥n || ArgumentError("`p` must be equal to or greater than `n` when using `RealCepstral`. \
  Passed $p and $n.") |> throw
  res = series |> fft .|> abs .|> log |> ifft
  return res[1:n]
end
