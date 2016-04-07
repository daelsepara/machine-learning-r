lr_cost_op <- function(Xp, y, theta) {
# Computes X*theta - y, an intermediate result used in computing
# for the cost function and the gradient
#
# Inputs:
#  Xp[m, n]	training set
#      y[m]	true values for X
#  theta[n]	model parameters
#  
# Outputs:
#   Computed X*theta-y

	result = Xp %*% theta - y
	
	return(result)
}

lr_cost <- function(X, y, theta) {
# Compute cost and gradient for linear regression with multiple variables
#
# Inputs:
# X[m, n-1]	training set
#      y[m]	true values for X
#  theta[n]	model parameters
#  
# Outputs:
#        J	linear regression cost function
# gradient	gradient of cost function with respect to theta

	m = nrow(y)
	n = nrow(theta)
	
	if (is.null(n)) {
		n = length(theta)
	}
	
	# add a bias column
	if (ncol(X) != n) {
		Xp = cbind(array(1, c(m, 1)), X)
	} else {
		Xp = X
	}

	result = lr_cost_op(Xp, y, theta) 
	
	# compute cost
	J = sum(result ^ 2)/(2 * m)

	# compute gradient
	gradient = t(Xp) %*% result / m
	
	return(list('J' = J, 'gradient' = gradient))
}

lr_gradientdescent <- function(X, y, theta, alpha, num_iters) {
# Performs gradient descent to optimize theta. Runs gradient descent 
# algorithm for num_iters steps with using learning rate alpha
#
# Inputs:
#   X[m, n]	training set
#      y[m]	true values for X
#  theta[n]	model parameters
#     alpha	learning rate
# num_iters	maximum number of iterations
#  
# Outputs:
#     theta	optimum theta
#         J	cost function history

	m = nrow(X)
	J = array(0, c(num_iters, 1))
	
	theta_optimum = theta
	
	for (i in 1:num_iters) {
		# compute cost and gradient
		result = lr_cost(X, y, theta_optimum)

		# record cost function history
		J[i] = result$J
		
		# perform gradient descent
		theta_optimum = theta_optimum - alpha * result$gradient
	}
	
	return(list('J' = J, 'theta' = theta_optimum))
}

sigmoid <- function(z) {
# Compute sigmoid function for z (scalar, vector, matrix)
#
# Inputs:
#      z	(scalar, vector, matrix)
#
# Outputs:
#   h(z)	sigmoid function of z

	h = 1/(1 + exp(-z))
	
	return(h)
}

softplus <- function(z) {
# Compute softplus function for z (scalar, vector, matrix)
#
# Inputs:
#      z	(scalar, vector, matrix)
#
# Outputs:
#   h(z)	softplus function of z

	h = log(1 + exp(z))
	
	return(h)
}

logr_cost <- function(X, y, theta, lambda = 0) {
# Compute cost and gradient for logistic regression with regularization parameter lambda
#
# Inputs:
#  X[m, n-1]	training set
#       y[m]	true values for X
#   theta[n] 	model parameters
#     lambda	regularization paramter
#  
# Outputs:
#          J	linear regression cost function
#   gradient	gradient of cost function with respect to theta

	m = nrow(y)
	n = nrow(theta)
	
	if (is.null(n)) {
		n = length(theta)
	}
	
	# add a bias column
	if (ncol(X) != n) {
		Xp = cbind(array(1, c(m, 1)), X)
	} else {
		Xp = X
	}
	
	# do not regularize theta[1]
	reg_theta = theta
	reg_theta[1] = 0
	
	# compute sigmoid factors
	h = sigmoid(Xp %*% theta)

	# compute cost function
	J = sum(-y * log(h) - (1 - y) * log(1 - h)) / m  + lambda * sum(reg_theta ^ 2)/(2 * m)
	
	# compute gradient
	gradient = (t(Xp) %*% (h - y) + lambda * reg_theta) / m
	
	return(list('J' = J, 'gradient' = gradient))
}

logr_optimize <- function(X, y, theta, lambda = 0, num_iters = 100, method = 'L-BFGS-B') {
# Compute optimum regularized model parameters using R's optimizer and
# logistic regression cost function
#
# Inputs:
#    X[m, n]	training set
#       y[m]	true values for X	
#     lambda	regularization parameter
#   theta[n]	model parameters
#  num_iters	maximum number of iterations
#     method	Optimization method to use: 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'SANN', 'Brent'
#
# Outputs:
#	   theta	optimum theta
#
# See: ?optim

	return(regression_optimize(logr_cost, X, y, theta, lambda, num_iters, method))
}

softmax_cost <- function(X, y, theta, lambda = 0) {
# Compute cost and gradient for softmax regression with regularization parameter lambda
#
# Inputs:
#    X[m, n]	training set
#       y[m]	true values for X	
#theta[n, k]	model parameters (n features, k-classes)
#     lambda	regularization parameter
#  
# Outputs:
#          J	softmax regression cost function
#   gradient	gradient of cost function with respect to theta
#

  # for repmat
  require(pracma)
  
  X = t(X)
  
  m = ncol(X)
  n = nrow(X)

  # reshape theta into a matrix
  theta = array(theta, c(n, length(theta) / n))
  
  k = ncol(theta)
  num_classes = k + 1
  
  # numerically stable way of computing softmax
  tx = rbind(t(theta) %*% X, array(0, c(1, m)))
  z = exp(tx - repmat(array(apply(tx, 2, max), c(1, m)), num_classes, 1))
  h = z / repmat(apply(z, 2, sum), num_classes, 1)
  
  # determine index where y == k
  ind = (1:m - 1) * num_classes + y
  
  # compute regularized cost function
  J = -sum(log(h[ind])) + lambda * sum(theta ^ 2) / 2
  
  yk = array(0, c(num_classes, m))
  yk[ind] = 1

  # ignore last row
  rows = 1:(num_classes - 1)
  
  yk = yk[rows, ]
  h = h[rows, ]

  # compute gradient with regularization
  gradient = -as.vector(X %*% t(yk - h) - lambda*theta)

  return(list('J' = J, 'gradient' = gradient))
}

softmax_optimize <- function(X, y, theta, lambda = 0, num_iters = 100, method = 'L-BFGS-B') {
# Compute regularized optimum model parameters using R's optimizer and
# softmax regression cost function
#
# Inputs:
#    X[m, n]	training set
#       y[m]	true values for X
#theta[n, k]	model parameters (n-features, k-classes)
#     lambda	regularization parameter
#  num_iters	maximum number of iterations
#     method	Optimization method to use: 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'SANN', 'Brent'
#
# Outputs:
#	   theta	optimum theta
#
# See: ?optim
  
  return(regression_optimize(softmax_cost, X, y, theta, lambda, num_iters, method))
}

softplus_cost <- function(X, y, theta, lambda = 0) {
# Compute cost and gradient for softplus regression
#
# Inputs:
#    X[m, n]	training set
#       y[m]	true values for X	
#theta[n, k]	model parameters (n features, k-classes)
#     lambda	regularization parameter
#  
# Outputs:
#          J	softmax regression cost function
#   gradient	gradient of cost function with respect to theta
#

	m = nrow(X)
	
	xt = X %*% theta

	result = softplus(xt) - y
	
	# compute cost
	J = sum(result ^ 2)/(2 * m)  + lambda * sum(theta ^ 2)/(2 * m)

	gradient = as.vector(t(X) %*% (result * sigmoid(xt)) / m + lambda * theta / m)
	
	return(list('J' = J, 'gradient' = gradient))
}

softplus_optimize <- function(X, y, theta, lambda = 0, num_iters = 100, method = 'L-BFGS-B') {
# Compute regularized optimum model parameters using R's optimizer and
# softplus regression cost function
#
# Inputs:
#    X[m, n]	training set
#       y[m]	true values for X
#theta[n, k]	model parameters (n-features, k-classes)
#     lambda	regularization parameter
#  num_iters	maximum number of iterations
#     method	Optimization method to use: 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'SANN', 'Brent'
#
# Outputs:
#	   theta	optimum theta
#
# See: ?optim
  
  return(regression_optimize(softplus_cost, X, y, theta, lambda, num_iters, method))

}


regression_optimize <- function(f, X, y, theta, lambda = 0, num_iters = 100, method = 'L-BFGS-B') {
# Compute optimum regularized model parameters using R's optimizer
#
# Inputs:
#          f	cost/gradient function
#    X[m, n]	training set
#       y[m]	true values for X	
#     lambda	regularization parameter
#  num_iters	maximum number of iterations
#     method	Optimization method to use: 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'SANN', 'Brent'
#
# Outputs:
#	   theta	optimum theta
#
# See: ?optim

	# optim works with functions with one argument/parameter. We define anonymous functions (which are just wrappers to our cost function) to acheive the desired effect
	result = optim(par = theta, fn = function(theta) { return(f(X, y, theta, lambda)$J) }, gr = function(theta) { return(f(X, y, theta, lambda)$gradient) }, control = list('maxit' = num_iters), method = method)
	
	cat('\nResults:\n')
	cat(paste('J =', result$value, '\n'))
	cat(paste('iterations =', result$counts[1], '\n'))
	
	theta_optimum = result$par
	
	return(theta_optimum)
}

logr_predict <- function(theta, X, threshold = 0.5) {
# Predicts whether the X is 0 or 1 using learned logistic 
# regression parameters theta
#

	# number of training examples
	m = nrow(X)

	# compute sigmoid values
	h = sigmoid(X %*% theta)
	
	pos = which(h >= threshold)

	# set 1 on positive labels
	p = array(0, c(m, 1))
	p[pos] = 1
	
	return(p)
}

softmax_predict <- function(theta, X) {
# Multiclass prediction for X using learned softmax
# regression parameters theta

	# for repmat
	require(pracma)
	
	# check number of classes
	k = ncol(theta)
	m = nrow(X)
	
	# numerically stable way of computing softmax
	xt = X %*% theta
	z = exp(xt - repmat(array(apply(xt, 1, max), c(m, 1)), 1, k))
	zs = array(apply(z, 1, sum), c(m, 1))
	h = z / repmat(zs, 1, k)
	
	# compute predicted class
	p = array(as.integer(apply(h, 1, which.max)), c(m, 1))

	return(p)
}

softplus_predict <- function(theta, X, threshold = 0.5) {
# Predicts whether the X is 0 or 1 using learned softplus 
# regression parameters theta
#

	# number of training examples
	m = nrow(X)

	# compute softplus values
	h = softplus(X %*% theta)
	
	pos = which(h >= threshold)

	# set 1 on positive labels
	p = array(0, c(m, 1))
	p[pos] = 1
	
	return(p)
}

mapFeature <- function(x1, x2, degree = 6) {
# Feature mapping function to polynomial features
#
# Note: x1 and x2 must be the same size
#
#   Returns a new feature array with more features, comprising of 
#   x1, x2, x1^2, x2^2, x1*x2, x1*x2^2, etc..
#
	polyFeatures = array(1, c(nrow(x1), 1))
	
	for (i in 1:degree){
		for (j in 0:i){
			polyFeatures = cbind(polyFeatures, (x1 ^ (i - j))*(x2 ^ j))
		}
	}
	
	return(polyFeatures)
}

regression_boundary <- function(X, y, theta) {
# Plots decision boundary learned by logistic regression
#
# Inputs:
#   X[m, n]	training set
#      y[m]	true values for X
#  theta[n]	model parameters

	# Plot the data
	regression_plot(X[,2:3], y)

	if (ncol(X) <= 3) {
		# Compute endpoints of decision boundary
		plot_x = c(min(X[,2])-2,  max(X[,2])+2)
		plot_y = (-1 / theta[3]) * (theta[2] * plot_x + theta[1])

		# Plot decision boundary
		plot(plot_x, plot_y, col = 'green')
		
	} else {

		u = seq(min(X[,2:3])-1, max(X[,2:3])+1, length=100)
		v = seq(min(X[,2:3])-1, max(X[,2:3])+1, length=100)
		
		z = array(0, c(length(u), length(v)))

		for (i in 1:length(u)) {
			for (j in 1:length(v)) {
			   z[i, j] = mapFeature(u[i], v[j]) %*% theta
			}
		}
		
		contour(u, v, t(z), col = 'green', add = TRUE, lw = 1)
	}
}

regression_plot <- function(X, y) {
# Plots the data points X and y into a new figure 
#   regression_plot(x,y) plots the data points with + for the positive examples
#   and o for the negative examples. X is assumed to be a mx2 matrix.
#
# Inputs:
#   X[m, n]	training set
#      y[m]	true values for X

	# Find Indices of Positive and Negative Examples
	pos = which(y == 1)
	neg = which(y == 0)

	xlim = c(min(X), max(X))
	ylim = c(min(X), max(X))
	
	# Plot Examples
	plot(X[pos, 1], X[pos, 2], pch = 3, xlab = '', ylab = '', col = 'blue', ylim = ylim, xlim = xlim)
	points(X[neg, 1], X[neg, 2], pch = 1, col = 'red')
}
