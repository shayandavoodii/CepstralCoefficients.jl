function cepscoef(p::Int, α::AbstractVector, n::Int)
  n>0 || ArgumentError("n must be positive") |> throw
  c = similar(α, n)
  if n==1
    c[1] = -α[1]
  elseif 1<n≤p
    secondterm = 0.
    for m ∈ 1:n-1
      secondterm += (1-m/n)*α[m]*c[n-m]
    end
    c[n] = -α[n] - secondterm
  else
    res = 0.
    for m ∈ 1:p
      res += (1-m/n)*α[m]*c[n-m]
    end
    c[n] = -res
  end
  return c
end
