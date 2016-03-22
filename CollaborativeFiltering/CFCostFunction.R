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

	gradX = array(0, dim(X))
	gradTheta = array(0, dim(Theta))

	# compute cost function with regularization (only includes products that have been rated by user)
	J = sum(R*(t(Theta %*% t(X)) - Y)^2)/2 + lambda*(sum(Theta^2) + sum(X^2))/2
	
	# compute gradients of X and Theta: dJ/dX, dJ/dTheta
	# modifies X, Theta if product has been rated by user
	gradX = (R*(t(Theta*t(X)) - Y)) %*% Theta + lambda*X
	gradTheta = (R*t(t(Theta*t(X)) - Y)) %*% X  + lambda*Theta
	
	return(list('J' = J, 'gradX' = gradX, 'gradTheta' = gradTheta))
}
