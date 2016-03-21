pca  <- function(X, center = TRUE) {
# Computes the principal component analysis using singular value decomposition function (svd)
#
# Inputs:
#   X[m, n] training set
#   center TRUE/FALSE, i.e. normalize X
#  
# Outputs:
#   PCA of X with: 
#
#   u, v = eigenvectors 
# 	d = singular values (standard deviation)
	
	if (center) {
		# center X so that the mean is ~ 0
		Xn = X - repmat(apply(X, 2, mean), nrow(X), 1)
	} else {
		Xn = X
	}
	
	# compute for the covariance matrix
	#
	# options for computing the covariance matrix: cov(), var()
	Sigma = t(Xn) %*% Xn / (nrow(Xn)-1)
		
	# returns the eigenvectors u,v and the singular values d
	#
	# options for computing the PCA: eigen(), although svd() is a more stable implementation
	# and eigen() becomes computationally-intensive as n-increases due to O(n^3) complexity
	#
	# prcomp() can also be used to compute the PCA using eigen() and cov()
	return(svd(Sigma))
}

pca_project <- function(X, U, K = ncol(U)) {
# Project X onto the U (principal components) using the first K-columns
#
# Inputs:
#   X[m, n] training set
#   U[n, n] eigenvectors of X
# 	K columns of reduced data, K <= n
#  
# Outputs:
#   Z projection of X on U using the first K-components
#

	m = nrow(X)
	Z = array(0, c(m, K))

	# project X onto U using only the first K-components
	Z = t(t(U[,1:K]) %*% t(X))
	
	return(Z)
}

pca_estimate <- function(Z, U, K = ncol(U)) {
# Compute approximation of the original data using the projection data
# Project X onto the U (principal components) using the first K-columns
#
# Inputs:
#   Z[m, K] projection
#   U[n, n] eigenvectors of X
# 	K columns of reduced data, K <= n
#  
# Outputs:
#   Xa[m, n] approximation of X using projection data
#
    Xa = t(U[,1:K] %*% t(Z))

    return(Xa)
}
