function cepscoef(p::Int, α::AbstractVector, n::Int)
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
