lr_cost_op <- function(Xp, y, theta) {
# Computes X*theta - y, an intermediate result used in computing
# for the cost function and the gradient
#
# Inputs:
#   Xp[m, n] training set
#   y[m] true values for X
#   theta[n] model parameters
#  
# Outputs:
# 	X*theta-y

	result = Xp %*% theta - y
	
	return(result)
}

lr_cost <- function(X, y, theta) {
# Compute cost and gradient for linear regression with multiple variables
#
# Inputs:
#   X[m, n-1]	training set
#   y[m] 		true values for X
#   theta[n] 	model parameters
#  
# Outputs:
# 	J			linear regression cost function
#  gradient		gradient of cost function with respect to theta

	# number of training examples
	m = nrow(y)

	# add a bias column
	Xp = cbind(array(1, c(m, 1)), X)

	result = lr_cost_op(Xp, y, theta) 
	
	# compute cost
	J = sum(result ^ 2)/(2 * m)

	# compute gradient
	gradient = t(Xp) %*% result / m
	
	return(list('J' = J, 'gradient' = gradient))
}

lr_gradientdescent <- function(X, y, theta, alpha, num_iters) {
# Performs gradient descent to optimize theta
# Runs gradient descent algorithm for num_iters steps with using learning
# rate alpha
#
# Inputs:
#   X[m, n] 	training set
#   y[m] 		true values for X
#   theta[n] 	model parameters
#	alpha		learning rate
#	num_iters	maximum number of iterations
#  
# Outputs:
#	theta		optimum theta
# 	J 			cost function history

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
