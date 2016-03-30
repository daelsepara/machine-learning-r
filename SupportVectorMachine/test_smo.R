test_smo <- function() {
  
  l = 2	# Dimensionality
  N = 150	# Number of vectors
  
  # Generating the training set
  set.seed(0)
  
  X1 = 10 * rand(l, N) - 5
  y1 = array(y1, c(1, N))
  
  for (i in 1:N) {
    
    tt=0.05 * (X1[1, i]^3 + X1[1, i]^2 + X1[1, i] + 1)
    
    if (tt > X1[2, i]) {
      y1[i] = 1
    } else {
      y1[i] = -1
    }
  }
  
  krnel = 'linear'
  kpar1 = 0
  kpar2 = 0
  C = 2
  tol = 0.001
  steps = 100000
  eps = 10 ^ (-10)
  method = 1
  
  X1 = t(X1)
  y1 = t(y1)
  
  result = smo2(X = X1, Y = y1, krnel = krnel, kpar1 = kpar1, kpar2 = kpar2, C = C, tol = tol, steps = steps, eps = eps, method = method)

  return(result)
}