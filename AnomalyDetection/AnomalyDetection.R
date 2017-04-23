estimate_gaussian <- function(X) {
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
	
	return(list('mean' = mu, 'variance' = sigma2))
}

multivariate_gaussian <- function(X, mu = 0, variance = 1) {
# Computes the multivariate Gaussian distribution probability density function
#
# Inputs:
# X[m, n] data set (m-examples, 2-features)
#   mu[n] mean of each feature
#  var[n] variance of each feature
#
# Outputs:
#   p[m] multivariate gaussian distribution probability density function

  # for repmat and inv (if used instead of solve)
  require(pracma)
  
  # create diagonal matrix whose primary diagonal elements contains the variance
  if ((is.null(dim(variance)) && length(variance) > 0) || is.vector(variance)) {
    sigma2 = diag(variance)
  } else {
    sigma2 = variance
  }
  
  if (is.null(dim(mu))) {
    k = length(mu)
  } else {
    k = ncol(mu)
  }
  
  m = nrow(X)
  
  # center the PDF at the mean
  Xn = X - repmat(mu, m, 1)
  
  # inv() may be used instead of solve()
  p = (2*pi) ^ (-k/2) * det(sigma2) ^ (-0.5) * exp(-0.5*apply(Xn %*% solve(sigma2) * Xn, 1,sum))
  
  # reshape p into a column vector
  return(array(p, c(nrow(X), 1)))
}

select_threshold <- function(y, p) {
# Determine best threshold to use to detect outliers
#
# Inputs:
# y[m] test set (m sample)
# p[m] predicted values for test set X from estimated multivariate Gaussian PDF (using estimated parameters)
#
# Outputs:
#   epsilon: Best prediction threshold to use
#   F1: F1 score on epsilon
  
	epsilon_Best = 0
	F1_Best = 0

	# setup scan range
	stepsize = (max(p) - min(p)) / 1000
	
	for (epsilon in seq(min(p), max(p), by = stepsize)) {

		predictions = p < epsilon

		# determine true/false positives and negatives
		tp = sum((predictions == TRUE) & (y == 1))
		fp = sum((predictions == TRUE) & (y == 0))
		fn = sum((predictions == FALSE) & (y == 1))

		# compute precision and recall
		prec = tp/(tp + fp)
		rec = tp/(tp + fn)

		# compute F1 score based on precision and recall
		F1 =  2*prec*rec/(prec + rec)

		if (!is.nan(F1) && is.numeric(F1) && F1 > F1_Best) {
			F1_Best = F1
			epsilon_Best = epsilon
		}
	}
	
	return(list('F1' = F1_Best, 'epsilon' = epsilon_Best))
}
