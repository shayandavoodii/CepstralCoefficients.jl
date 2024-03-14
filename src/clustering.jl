"""
    cepsclustering(cc::AbstractMatrix, k::Integer)

Perform clustering of cepstral coefficients using k-medoids algorithm.

!!! note
    This function can be used only after importing the `Clustering` package.

# Arguments
- `cc::AbstractMatrix`: A matrix of cepstral coefficients.
- `k::Integer`: Maximum number of clusters to be examined (to find the optimal number of clusters).

# Returns
- A vector of vectors, where each vector contains the indices of the cepstral coefficients in \
  the corresponding cluster.
"""
function cepsclustering end
