function normalizer!(series::AbstractMatrix)
  n_obs, n_assets = size(series)
  for asset ∈ 1:n_assets
    normalizer!(@view series[:, asset])
  end
  return series
end

normalizer!(series::AbstractVector) = series .= (series .- mean(series)) ./ std(series)
