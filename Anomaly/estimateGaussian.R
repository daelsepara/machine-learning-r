estimateGaussian <- function(X) {
# estimates Gaussian distribution parameters of X
#
# Inputs:
#   X[m, n] data set (m-examples, 2-features)
#  
# Outputs:
#   mu[n] mean of each feature
#   var[n] variance of each feature (sigma^2)

	# for repmat
	require(pracma)
	
	m = nrow(X)
	
	# Compute the mean of each feature. This is for illustration purposes only. Best way to rewrite this is to use:
	#
	# mu = apply(X, 2, mean) ... instead
	#
	mu = apply(X, 2, sum)/m
	
	# Compute standard deviation of each feature
	#
	# this can be rewritten as:
	#
	# sigma2 = apply(X, 2, var)*(m - 1)/m
	#
	# or:
	#
	# sigma2 = apply(X, 2, sd))^2*(m - 1)/m
	#
	sigma2 = apply((X - repmat(mu, m, 1))^2, 2, sum)/m
	
	return(list('mu' = mu, 'var' = sigma2));
}