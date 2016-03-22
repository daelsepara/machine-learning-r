CFCostFunction <- function(X, Theta, Y, R, lambda) {
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
	gradTheta = t(R*(t(Theta %*% t(X)) - Y)) %*% X  + lambda*Theta;
	
	return(list('J' = J, 'X_grad' = gradX, 'Theta_grad' = gradTheta))
}

CFOptimizeCost <- function(theta, P1, P2, P3 , P4, P5, P6) {
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
	result = CFCostFunction(X, Theta, P1, P2, P3);
	
	# unroll gradient matrices into one vector
	grad = c(as.vector(result$X_grad), as.vector(result$Theta_grad))

	return(list('J' = result$J, 'grad' = grad))
}