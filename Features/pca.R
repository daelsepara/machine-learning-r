pca  <- function(X) {
# Computes the principal component analysis using singular value decomposition function (svd)
	
	# center X so that the mean is ~ 0
	Xn = X - repmat(apply(X, 2, mean), nrow(X), 1)
	
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
