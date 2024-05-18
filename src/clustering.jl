using StatsBase: sample
using Infiltrator

abstract type CepstralClustering end

struct Clustering{T<:AbstractMatrix}<:CepstralClustering
  centeroids::T
  umatrix::T
end

"""
    function cepsclustering(
      ccmat::AbstractMatrix,
      C::AbstractVector{<:Integer},
      m::AbstractFloat,
      maxiter::Integer,
      γ::Integer
    )
    cepsclustering(ccmat::AbstractMatrix, k::Integer)

# Method 1
```julia
function cepsclustering(
  ccmat::AbstractMatrix,
  C::AbstractVector{<:Integer},
  m::AbstractFloat,
  maxiter::Integer,
  γ::Integer
)
```

Perform PAM clustering method on the given cepstral coefficients.

## Arguments
- `ccmat::AbstractMatrix`: Cepstral coefficients matrix of size `nassets×k` where `k` is \
  the number of cepstral coefficients.
- `C::AbstractVector{<:Integer}`: A Vector of Integers including the number of clusters to \
  be examined.
- `m::AbstractFloat`: Fuzziness value.
- `maxiter::Integer`: Maximum number of iterations for the PAM to be performed.
- `γ::Integer`: Coefficient in calculating the Fuzzy Silhouette metric.

## Returns
- `::Clustering`: An object of `Clustering`.


# Method 2
Perform clustering of cepstral coefficients using k-medoids algorithm.

!!! note
    This function can be used only after importing the `Clustering` package.

## Arguments
- `ccmat::AbstractMatrix`: A matrix of cepstral coefficients.
- `k::Integer`: Maximum number of clusters to be examined (to find the optimal number of clusters).

!!! warning
    `ccmat` must be a matrix of cepstral coefficients, where each row represents a stock and \
    each column represents a cepstral coefficient.

## Returns
- A vector of vectors, where each vector contains the indices of the cepstral coefficients in \
  the corresponding cluster.
"""
function cepsclustering(
  ccmat::AbstractMatrix,
  C::AbstractVector{<:Integer},
  m::AbstractFloat,
  maxiter::Integer,
  γ::Integer
)
  minimum(C)>1 || ArgumentError("Minimum number of clusters in `C` should be more than one. \
  C contains $(minimum(C))."
  ) |> throw
  maxiter>0 || ArgumentError("Maximum number of iterations `maxiter` should be more than \
    zero. maxiter=$maxiter is passed."
  ) |> throw
  nassets, K = size(ccmat)
  Wₖ         = fill(1/K, K)
  limitfs    = -Inf
  global result = Clustering(rand(2, 2), rand(2, 2))
  for c ∈ C
    randidx      = sample(1:nassets, c, replace=false)
    centroids    = ccmat[randidx, :]
    centroidsₒₗd = similar(centroids)
    umat         = zeros(nassets, c)
    iter         = 0
    while iter < maxiter || centroids != centroidsₒₗd
      centroidsₒₗd = copy(centroids)
      uᵢcfunc!(ccmat, umat, centroids, Wₖ, m)
      # all(≈(1.), sum(umat, dims=2)) || @infiltrate
      wₖfunc!(Wₖ, ccmat, umat, centroids, m)
      centroids = updatecenter(centroidsₒₗd, ccmat, umat, Wₖ, m)
      fs        = fuzzysillhouette(umat, ccmat, γ)
      if fs>limitfs
        result = Clustering(centroids, umat)
        limitfs   = fs
      end
      iter += 1
    end
  end
  return result
end

function uᵢcfunc!(ccmat::AbstractMatrix, umat::AbstractMatrix, centroids::AbstractMatrix, Wₖ::AbstractVector, m::AbstractFloat)
  C = size(centroids, 1)
  for i ∈ axes(ccmat, 1)
    for c ∈ 1:C
      summation = 0.
      numer = sum((Wₖ.*(ccmat[i, :].-centroids[c, :])).^2)
      iszero(numer) && continue
      for c′ ∈ 1:C
        denumer = sum((Wₖ.*(ccmat[i, :].-centroids[c′, :])).^2)
        iszero(denumer) && continue
        summation += (numer/denumer)^(1/(m-1))
      end
      umat[i, c] = 1/summation
    end
  end
end

function wₖfunc!(Wₖ::AbstractVector, ccmat::AbstractMatrix, umat::AbstractMatrix, centroids::AbstractMatrix, m::AbstractFloat)
  I, K = size(ccmat)
  C = size(centroids, 1)
  for k ∈ 1:K
    summation = 0.
    for k′ ∈ 1:K
      numer   = 0.
      denumer = 0.
      for i ∈ 1:I
        for c ∈ 1:C
          numer   += umat[i, c]^m * (ccmat[i, k]-centroids[c, k])^2
          denumer += umat[i, c]^m * (ccmat[i, k′]-centroids[c, k′])^2
        end
      end
      summation += numer/denumer
    end
    Wₖ[k] = 1/summation
  end
end

function updatecenter(centroidsold::AbstractMatrix, ccmat::AbstractMatrix, umat::AbstractMatrix, Wₖ::AbstractVector, m::AbstractFloat)
  I, _      = size(ccmat)
  C         = size(centroidsold, 1)
  funcval   = zeros(I)
  centroids = similar(centroidsold)
  for c ∈ 1:C
    for i′ ∈ 1:I
      for i″ ∈ 1:I
        expr_ = @. (Wₖ * (ccmat[i″, :] - ccmat[i′, :]))^2
        funcval[i′] += umat[i″, c]^m * sum(expr_)
      end
    end
    @infiltrate
    centroids[c, :] = ccmat[argmin(funcval), :]
    zeros!(funcval)
  end
  return centroids
end

zeros!(x::AbstractVector{T}) where T = x.=zero(T)

membership(umat::AbstractMatrix) = last.(getproperty.(argmax(umat, dims=2), :I)) |> vec

function aᵢfunc(memvec::AbstractVector, ccmat::AbstractMatrix)
  a = similar(memvec, Float64)
  I = length(memvec)
  for i ∈ 1:I
    idxclusᵢ = memvec[i]
    idxinᵢ   = findall(memvec.==idxclusᵢ)
    a[i]     = distance(ccmat[i, :], ccmat[idxinᵢ, :])
  end
  return a
end

function distance(vec::AbstractVector, mat::AbstractMatrix)
  I    = size(mat, 1)
  vals = similar(mat, I)
  for idxrow ∈ axes(mat, 1)
    vals[idxrow] = √sum((vec.-mat[idxrow, :]).^2)
  end
  return sum(vals)/(I-1)
end

function bᵢfunc(memvec::AbstractVector, ccmat::AbstractMatrix)
  b           = similar(memvec, Float64)
  idxclusters = Set(memvec)
  nclusters   = length(idxclusters)
  I           = length(memvec)
  vals = similar(ccmat, nclusters)
  for i ∈ 1:I
    idxclusᵢ = memvec[i]
    for (j,idx) ∈ enumerate(idxclusters)
      idx==idxclusᵢ && continue
      idxinⱼ  = findall(memvec.==idx)
      vals[j] = distance(ccmat[i, :], ccmat[idxinⱼ, :])
    end
    b[i] = minimum(vals)
  end
  return b
end

function sᵢfunc(a::AbstractVector, b::AbstractVector)
  s = similar(a)
  for i ∈ eachindex(a)
    s[i] = (b[i]-a[i])/max(a[i], b[i])
  end
  return s
end

function _fs(umat::AbstractMatrix, S::AbstractVector, γ::Integer)
  I       = size(umat, 1)
  numer   = 0.
  denumer = 0.
  for i ∈ 1:I
    uᵢ, uⱼ   = sort(umat[i, :], rev=true)[1:2]
    val      = (uᵢ-uⱼ)^γ
    numer   += val * S[i]
    denumer += val
  end
  return numer/denumer
end

function fuzzysillhouette(umat::AbstractMatrix, ccmat::AbstractMatrix, γ::Integer)
  memvec = membership(umat)
  a      = aᵢfunc(memvec, ccmat)
  b      = bᵢfunc(memvec, ccmat)
  s      = sᵢfunc(a, b)
  return _fs(umat, s, γ)
end
