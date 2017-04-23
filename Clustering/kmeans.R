kmeans_initialize <- function(X, K = nrow(X)) {
# initialize centroids by randomly selecting K samples from X. K is the total number of clusters.
#
# Inputs:
#   X[m, n] training set
#   K number of centroids
#  
# Outputs:
#   centroids[K, n] initial centroids
  
  if (K <= nrow(X)) {
    # generate random permutations of X  
    Xpermutations = X[sample(1:nrow(X)),]
    
    # return only the first K rows
    return(Xpermutations[1:K,])
  }
}

kmeans_run <- function(X, initial_centroids, max_iter = 100) {
# Run k-means algorithm for X. X is an m x n matrix, where m is the number of examples, and n is the number of features.
#
# Inputs:
#   X[m, n] training set
#   initial_centroids[K,n] initial centroids
#   max_iter maximum number of iterations
#
# Outputs:
#   centroids[K, n] centroids
#   clusters[m] cluster assignment of each sample in X
  
  # initial conditions
  K = nrow(initial_centroids);
  centroids = initial_centroids;
  clusters = array(0, c(nrow(X), 1))
  
  for (i in 1:max_iter) {
    
    # Assign each example in X to the closest centroid
    clusters = kmeans_assign(X, centroids)
    
    # Compute new centroid locations
    centroids = kmeans_compute(X, clusters, K)
  }
  
  return(list('centroids' = centroids, 'clusters' = clusters))
}

kmeans_assign <- function(X, centroids) {
# Finds closest centroid to each example in X and assign cluster number
#
# Inputs:
#   X[m, n] training set
#   centroids[K,n] current centroids
#
# Output:
#   clusters[m] cluster assignment of each sample in X

  # for repmat
  require(pracma)
  
  # get data dimensions
  m = nrow(X)
  K = nrow(centroids)
  
  # initialize distance matrix
  distanceToCentroid = array(0, c(m, K))
  
  for (i in 1:K) {
    # compute distance to each cluster
    centroidMatrix = repmat(centroids[i,], m, 1)
    
    # colSums/rowSums output is not an array, so we need to recover the correct shape
    distanceToCentroid[, i] = sqrt(array(rowSums((X-centroidMatrix)^2), c(m, 1)))
  }
  
  # select cluster number with the minimum distance and assign each sample to that cluster
  clusters = array(as.integer(apply(distanceToCentroid, 1, which.min)), c(m, 1))
  
  return(clusters)
}

kmeans_compute <- function(X, clusters, K = nrow(X)) {
# Computes new centroid locations
#
# Inputs:
#   X[m, n] training set
#   clusters[m] cluster assignment
#
# Output:
#   centroids[K, n] new locations of K-centroids

  # initialize centroids
  n = ncol(X)
  centroids = array(0, c(K, n))
  
  for (k in 1:K) {
    # compute mean from samples assigned to the k-th cluster
    pts = which(clusters == k)
    
    # colSums/rowSums output is not an array, so we need to recover the correct shape
    Xsum = array(colSums(X[pts,]), c(1, n))
    centroids[k, ] = Xsum/length(pts)
  }
  
  return(centroids)
}