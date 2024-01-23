function fit_arima(series::AbstractVector, p::Int)
  model = fit(ARMA{p, 0}, series)
  return model.meanspec.coefs[2:end]
end

function fit_arima(series::AbstractMatrix, p::Int)
  p>0 || ArgumentError("p must be positive") |> throw
  n_obs, n_assets = size(series)
  coefs = similar(series, p, n_assets)
  for asset ∈ 1:n_assets
    coefs[:, asset] = fit_arima(series[:, asset], p)
  end
  return coefs
end
