estimateGaussian <- function(X) {
# estimates Gaussian distribution parameters of X
#
# Inputs:
#   X[m, n] data set (m-examples, 2-features)
#  
# Outputs:
#   mu[n] mean of each feature
#   var[n] variance of each feature (sigma^2)

	m = nrow(X)
	
	# Compute the mean of each feature
	#
	mu = apply(X, 2, mean)
	
	# Compute standard deviation of each feature
	#
	sigma2 = apply(X, 2, var)*(m - 1)/m
	
	return(list('mu' = mu, 'var' = sigma2));
}