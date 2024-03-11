normalizer(series::AbstractVector) = (series .- mean(series)) ./ std(series)

mean(series::AbstractVector) = sum(series) / length(series)
