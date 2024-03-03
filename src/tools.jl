normalizer!(series::AbstractVector) = series .= (series .- mean(series)) ./ std(series)

@inline mean(series::AbstractVector) = sum(series) / length(series)
