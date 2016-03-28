euclidean_classifier <- function(X, y) {
# Euclidean classifier
#
# Inputs:
#   X[m, n]	data to be classified [m data, n features]
#   y[n, k] class feature vectors [n features (mean), k classes]
#
# Outputs:
#      z[m]	classification of m-samples (labels) 

	# for repmat
	require(pracma)
	
	n = nrow(y)
	k = ncol(y)
	m = nrow(X)
	
	z = array(0, c(m))
	d_e = array(0, c(k, 1))
	
	for (i in 1:m) {
		# replicate i-th column of X into into k-columns
		d = repmat(array(X[i, ], c(n, 1)), 1, k) - y
		
		# compute Euclidean distance to each class and discard off-diagonals
		d_e = sqrt(diag(t(d) %*% d))
		
		# classify according to minimum Euclidean distance
		z[i] = as.integer(which.min(d_e))
	}
	
	return(z)
}
