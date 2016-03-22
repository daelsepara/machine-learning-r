cf_costfunction <- function(X, Theta, Y, R, lambda) {
# Collaborative filtering cost function
#
# Inputs:
#     X[m, n] product features (m - products, n - features)
# Theta[m, n] model parameters, user features (m - users, n - features)
#     Y[m, n] user ratings (m - products, n - users)
#     R[m, n] = 1 if user rated product, 0 otherwise (m - products, n - users)
#
# Outputs:
#   J cost function of collaborative filtering
#   gradX, gradTheta gradients for use in training/learning/regression

	# compute cost function with regularization (only includes products that have been rated by user)
	J = sum(R*(t(Theta %*% t(X)) - Y)^2)/2 + lambda*(sum(Theta^2) + sum(X^2))/2
	
	# compute gradients of X and Theta: dJ/dX, dJ/dTheta
	# modifies X, Theta if product has been rated by user
	gradX = (R*(t(Theta %*% t(X)) - Y)) %*% Theta + lambda*X
	gradTheta = t(R*(t(Theta %*% t(X)) - Y)) %*% X  + lambda*Theta
	
	return(list('J' = J, 'X_grad' = gradX, 'Theta_grad' = gradTheta))
}

cf_fmincg_cost <- function(theta, P1, P2, P3 , P4, P5, P6) {
# Cost function for use with advanced optimization method (fmincg)
  
	# theta (X and Theta)
	# P1 Y
	# P2 R
	# P3 lambda
	# P4 number of users
	# P5 number of products
	# P6 number of features
	
	# roll up vectors into arrays
	offs = P5*P6
	X = array(theta[1:offs], c(P5, P6))
	Theta = array(theta[(offs+1):length(theta)], c(P4, P6))
	result = cf_costfunction(X, Theta, P1, P2, P3)
	
	# unroll gradient matrices into one vector
	grad = c(as.vector(result$X_grad), as.vector(result$Theta_grad))

	return(list('J' = result$J, 'grad' = grad))
}

cf_optimize <- function(X, Theta, Y, R, lambda, iterations = 100) {
# Algorithm to run fmincg with CFCostFunction
# Collaborative filtering cost function
#
# Inputs:
#     X[m, n] product features (m - products, n - features)
# Theta[m, n] model parameters, user features (m - users, n - features)
#     Y[m, n] user ratings (m - products, n - users)
#     R[m, n] = 1 if user rated product, 0 otherwise (m - products, n - users)
#      lambda regularization parameter
#  iterations maximum number of iterations to run optimizer
#
# Outputs:
#    X, Theta  optimzed product and model features/parameters
#    	 cost  final cost
#  iterations  actual number of iterations used
	
	# determine data dimensions
	num_users = ncol(X)
	num_products = nrow(X)
	num_features = ncol(Theta)

	# initialize parameters with random values
	randX = array(rnorm(n = num_products * num_users, mean = 0, sd = 1), c(num_products, num_users))
	randTheta = array(rnorm(n = num_users * num_features, mean = 0, sd = 1), c(num_users, num_features))
	
	# unroll into a single vector
	initial_params = c(as.vector(randX), as.vector(randTheta))
	
	# optimize cost function
	result = fmincg(cf_fmincg_cost, initial_params, iterations, Y, R, lambda, num_users, num_products, num_features)
	
	# roll up result into arrays
	offs = num_products*num_features
	optX = array(result$X[1:offs], c(num_products, num_features))
	optTheta = array(result$X[(offs+1):length(result$X)], c(num_users, num_features))
	
	return(list('X' = optX, 'Theta' = optTheta, 'iterations' = result$iterations, 'cost' = result$cost))
}