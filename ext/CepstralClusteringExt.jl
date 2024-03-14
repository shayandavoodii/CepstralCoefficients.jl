module CepstralClusteringExt

using CepstralCoefficients
using Clustering
using Distances

function CepstralCoefficients.cepsclustering(ccmat::AbstractMatrix, k::Integer)
  k>size(ccmat, 1) && ArgumentError("`k` must be less than the number of stocks \
    ($(size(ccmat, 1)))."
  ) |> throw
  dists   = pairwise(Euclidean(), ccmat, dims=1)
  kopt    = kclusopt(dists, k)
  fitted  = kmedoids(dists, kopt)
  indices = assignments(fitted)
  return [findall(indices .== i) for i=1:kopt]
end

function kclusopt(data::AbstractMatrix, k::Integer)
  sils = zeros(Float64, k)
  for nclus ∈ 2:k
    fitted        = kmedoids(data, nclus)
    sils[nclus-1] = silhouettes(fitted, data) |> CepstralCoefficients.mean
  end
  return argmax(sils) + 1
end

end # module
