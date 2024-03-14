"""
    cepsclustering(ccmat::AbstractMatrix, k::Integer)

Perform clustering of cepstral coefficients using k-medoids algorithm.

!!! note
    This function can be used only after importing the `Clustering` package.

# Arguments
- `ccmat::AbstractMatrix`: A matrix of cepstral coefficients.
- `k::Integer`: Maximum number of clusters to be examined (to find the optimal number of clusters).

!!! warning
    `ccmat` must be a matrix of cepstral coefficients, where each row represents a stock and \
    each column represents a cepstral coefficient.

# Returns
- A vector of vectors, where each vector contains the indices of the cepstral coefficients in \
  the corresponding cluster.
"""
function cepsclustering end
