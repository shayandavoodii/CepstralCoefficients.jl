function fit_arima(series::AbstractVector{T}, p::Integer) where T
  model = fit(ARMA{p, 0}, series)
  return model.meanspec.coefs[2:end]::Vector{T}
end
