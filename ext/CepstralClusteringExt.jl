module CepstralClusteringExt

using CepstralCoefficients
using Clustering
using Distances

function CepstralCoefficients.cepsclustering(cc::AbstractMatrix, k::Integer)
  dists   = pairwise(Euclidean(), cc, dims=1)
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
