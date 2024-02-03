function fit_arima(series::AbstractVector{T}, p::Integer, q::Integer=0) where T
  model = fit(ARMA{p, q}, series)
  return model.meanspec.coefs[2:end]::Vector{T}
end
