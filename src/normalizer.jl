normalizer!(series::AbstractVector) = series .= (series .- mean(series)) ./ std(series)
