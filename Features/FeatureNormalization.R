featureNormalize <- function(X) {
# Normalize features of X. X is a matrix of m-examples and n-features. Each feature will have a mean of 0
# and a standard deviation of 1.
  
  # bsxfun and repmat are implemented in pracma library
  require(pracma)
  
  m = nrow(X)
  
  # get mean of each feature and subtract it from all elements of X
  X_norm = bsxfun('-', X, repmat(apply(X, 2, mean), m, 1))
  
  # get standard deviation and dividate it from all elements of X
  X_norm = bsxfun('/', X, repmat(apply(X, 2, sd), m, 1))
  
  return(X_norm)
}
