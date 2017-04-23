euclidean_classifier <- function(X, y) {
# Euclidean classifier
#
# Inputs:
#   X[m, n]	data to be classified [m samples, n features]
#   y[n, k]	class feature vectors [n features (mean), k classes]
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

mahalanobis_classifier <- function(X, y, S) {
# Mahalanobis classifier
#
# Inputs:
#   X[m, n]	data to be classified [m samples, n features]
#   y[n, k]	class feature vectors [n features (mean), k classes]
#   S[n, n]	covariance matrix of n-features
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
		
		# compute Mahalanobis distance to each class and discard off-diagonals
		d_e = sqrt(diag(t(d) %*% solve(S) %*% d))
		
		# classify according to minimum Mahalanobis distance
		z[i] = as.integer(which.min(d_e))
	}
	
	return(z)
}

perceptron_classifier <- function(X, y, w_init = array(1, c(ncol(X), 1)), alpha = 1, max_iter = 100000) {
# Perceptron classifier
#
# Inputs:
#   X[m, n]	data to be classified [m samples, n features]
#      y[m]	class labels of X
# w_init[n]	initial estimate of the parameter vector
#     alpha	learning rate
#
# Outputs:
#w_final[n]	final estimate of the parameter vector
#      iter	number of iterations ran until convergence
#        mc	number of misclassified samples
#      z[m] predicted classes

	m = nrow(X)
	n = ncol(X)
	
	iter = 0
	w_final = w_init
	mc = m
	
	while (mc > 0 && iter < max_iter) {
    
		iter = iter + 1

		# vectorized implementation of perceptron algorithm
		cost = (X %*% w_final) * y
		idx = which(cost < 0)
		
		# count misclassified samples
		mc = length(idx)
		
		# compute gradient
		gradient = alpha * t(-y[idx] %*% X[idx, ])
		
		if (iter == 1) {
			cat(paste('\n First Iteration: # Misclassified points = ', mc, '\n'))       
		}
		
		# Update estimates of parameter vector
		w_final = w_final - alpha * gradient
	}

	z = array(0, dim(y))
	p = X %*% w_final
	z[which(p < 0)] = -1
	z[which(p >= 0)] = 1
	
	return(list('z' = z, 'w_final' = w_final, 'mc' = mc, 'iter' = iter))
}

online_perceptron_classifier <- function(X, y, w_init = array(1, c(ncol(X), 1)), alpha = 1, max_iter = 100000) {
# Online form of the perceptron classifier
#
# Inputs:
#   X[m, n]	data to be classified [m samples, n features]
#      y[m]	class labels of X
# w_init[n]	initial estimate of the parameter vector
#     alpha	learning rate
#
# Outputs:
#w_final[n]	final estimate of the parameter vector
#      iter	number of iterations ran until convergence
#        mc	number of misclassified samples
#      z[m] predicted classes

	m = nrow(X)
	n = ncol(X)
	
	iter = 0
	w_final = w_init
	mc = m
	
	while (mc > 0 && iter < max_iter) {
		
		mc = 0
		
		for (i in 1:m) {
		
			if ((X[i, ] %*% w_final) * y[i] < 0) {
				# count misclassified samples
				mc = mc + 1
				# update estimates of parameter vector
				w_final = w_final + alpha * y[i] * X[i, ]
			}
			
			iter = iter + 1
		}
		
		if (iter == m) {
			cat(paste('\n First-pass: # Misclassified points = ', mc, '\n'))       
		}
	}

	z = array(0, dim(y))
	p = X %*% w_final
	z[which(p < 0)] = -1
	z[which(p >= 0)] = 1
	
	return(list('z' = z, 'w_final' = w_final, 'mc' = mc, 'iter' = iter))
}

sse_classifier <- function(X, y, C) {
# Sum-squared error classifier
#
# Inputs:
#   X[m, n]	data to be classified [m samples, n features]
#      y[m]	class labels of X
#         C	small positive constant that guarantees the
#           inversion of X' %*% X, when it is near singular
#
# Outputs:
#w_final[n]	final estimate of the parameter vector
#      z[m] predicted classes

	w_final = solve(t(X) %*% X + C*diag(ncol(X))) %*% (t(X) %*% y)

	z = array(0, dim(y))
	p = X %*% w_final
	z[which(p < 0)] = -1
	z[which(p >= 0)] = 1
	
	return(list('z' = z, 'w_final' = w_final))
}
