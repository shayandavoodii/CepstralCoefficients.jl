function fit_arima(series::AbstractVector, p::Int)
  model = fit(ARMA{p, 0}, series)
  return model.meanspec.coefs[2:end]
end

function fit_arima(series::AbstractMatrix, order::Int)
  order>0 || ArgumentError("Order must be positive") |> throw
  n_assets, n_obs = size(series)
  coefs = similar(series, n_assets, order)
  for asset ∈ 1:n_assets
    coefs[asset, :] = fit_arima(series[asset, :], order)
  end
  return coefs
end
