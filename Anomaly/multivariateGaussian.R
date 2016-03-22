multivariateGaussian <- function(X, mu = 0, variance = 1) {
# Computes the multivariate gaussian distribution probability density function
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