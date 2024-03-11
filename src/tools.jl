normalizer!(series::AbstractVector) = series .= (series .- mean(series)) ./ std(series)
normalizer(series::AbstractVector) = (series .- mean(series)) ./ std(series)

mean(series::AbstractVector) = sum(series) / length(series)
