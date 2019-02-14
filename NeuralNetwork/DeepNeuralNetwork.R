# Extends basic Neural Network, now has 5 layers (input, 3-hidden layers, output)

# Sigmoid activation function
nnet_sigmoid <- function(x) {
  
  return(1/(1 + exp(-x)))
}

# 1st-derivative of sigmoid activation function
nnet_dsigmoid <- function(x) {
  
  z = nnet_sigmoid(x)
  
  return(z * (1 - z))
}

# Forward propagation
nnet_forward <- function(input, w1, w2, w3, w4) {
  
  # add bias column to input layer
  x1 = cbind(array(1, c(nrow(input), 1)), input)
  
  # compute 1st hidden layer activation
  z1 = x1 %*% t(w1)
  s1 = nnet_sigmoid(z1)
  
  # add bias column to 1st hidden activation
  x2 = cbind(array(1, c(nrow(s1), 1)), s1)

  # compute 2nd hidden layer activation
  z2 = x2 %*% t(w2)
  s2 = nnet_sigmoid(z2)

  # add bias column to 2nd hidden activation
  x3 = cbind(array(1, c(nrow(s2), 1)), s2)
  
  # compute 3rd hidden layer activation
  z3 = x3 %*% t(w3)
  s3 = nnet_sigmoid(z3)

  # add bias column to 3rd hidden layer activation
  x4 = cbind(array(1, c(nrow(s3), 1)), s3)
  
  # compute output layer activation
  z4 = x4 %*% t(w4)
  yk = nnet_sigmoid(z4)
  
  return(list('yk' = yk, 'z1' = z1, 'z2' = z2, 'z3' = z3, 'x1' = x1, 'x2' = x2, 'x3' = x3, 'x4' = x4))
}

# Backward propagation
nnet_backprop <- function(yk, z1, z2, z3, x1, x2, x3, x4, w1, w2, w3, w4, y_matrix) {
  
  m = nrow(x1)
  
  # compute intermediate delta values per layer
  d4 = yk - y_matrix
  d3 = d4 %*% w4[, 2:ncol(w4)] * nnet_dsigmoid(z3)
  d2 = d3 %*% w3[, 2:ncol(w3)] * nnet_dsigmoid(z2)
  d1 = d2 %*% w2[, 2:ncol(w2)] * nnet_dsigmoid(z1)
  
  dw4 = t(d4) %*% x4
  dw3 = t(d3) %*% x3
  dw2 = t(d2) %*% x2
  dw1 = t(d1) %*% x1
  
  cost = sum(-y_matrix * log(yk) - (1 - y_matrix) * log(1 - yk))
    
  cost = cost / m
  dw4 = dw4 / m
  dw3 = dw3 / m
  dw2 = dw2 / m
  dw1 = dw1 / m
  
  return(list('dw1' = dw1, 'dw2' = dw2, 'dw3' = dw3, 'dw4' = dw4, 'Error' = cost))	
}

# create labels for multi-class classification
nnet_labels <- function(output, num_labels) {
  
  # For multi-classification problem, format expected output
  # i.e. matrix, each row corresponds to a training pattern.
  # Each element in the row-vector is a 0 or 1 indicating whether
  # or not it belongs to that particular class  
  if (num_labels > 1) {
    eye_matrix = diag(num_labels)
    y_matrix = eye_matrix[output, ]
  } else {
    # binary classification
    y_matrix = output
  }
  
  return(y_matrix)
}

# intialize interconnection weights with random values (-min_max, min_max) or Gaussian (mean = 0, sd = min_max)
nnet_weights <- function(min_max = 1, m = 1, n = 1, isGaussian = FALSE) {
  
  if (!isGaussian) {
    
    return(array(runif(n = m * n, min = -min_max, max = min_max), c(m, n)))
    
  } else {
    
    return(array(rnorm(n = m * n, mean = 0, sd = abs(min_max)), c(m, n)))
  }  
}

nnet_predict <- function(test_set, w1, w2, w3, w4, threshold = 0.5) {
  # Predict using neural network parameters (multi-class classification)
  
  prediction_output = nnet_forward(test_set, w1, w2, w3, w4)$yk
  
  m = nrow(test_set)
  
  prediction = array(0, c(m, 1))
  
  if (ncol(prediction_output) > 1) {
    # for multi-class neural network classifier, each column in
    # the output correspond to a different class. The node (in the output layer)
    # with the highest output value corresponds to its predicted class
    prediction = array(as.integer(apply(prediction_output, 1, which.max)), c(m, 1))
    
  } else {
    # for binary classifier, use threshold to set the output to 0 or 1
    prediction[which(prediction_output > threshold)] = 1
  }
  
  return(prediction)
}

# Network training using gradient descent
nnet_train <- function(maxiter = 100, learning_rate = 0.1, tol = 10^(-3), training_set = array(0) , output = array(0), hidden_units = 0, num_labels = 1, min_max = 1, isGaussian = FALSE) {
  
  # determine network dimensions from user input
  j = hidden_units
  inputs = ncol(training_set)
  
  # initialize weights with random values
  w1 = nnet_weights(min_max, j, inputs + 1, isGaussian)
  w2 = nnet_weights(min_max, j, j + 1, isGaussian)
  w3 = nnet_weights(min_max, j, j + 1, isGaussian)
  w4 = nnet_weights(min_max, num_labels, j + 1, isGaussian)
  
  iter = 0
  Error = 1.0
  
  yk = numeric(0)
  
  m = nrow(training_set)
  y_matrix = nnet_labels(output, num_labels)
  
  while (!is.nan(Error) && (iter < maxiter && Error > tol)) {
    
    # for training, perform forward and backpropagation each iteration
    forward = nnet_forward(training_set, w1, w2, w3, w4)
    backward = nnet_backprop(forward$yk, forward$z1, forward$z2, forward$z3, forward$x1, forward$x2, forward$x3, forward$x4, w1, w2, w3, w4, y_matrix)
    
    dw4 = learning_rate * backward$dw4
    dw3 = learning_rate * backward$dw3
    dw2 = learning_rate * backward$dw2
    dw1 = learning_rate * backward$dw1
    
    # update weights (using learning rate and gradient descent)
    w4 = w4 - dw4
    w3 = w3 - dw3
    w2 = w2 - dw2
    w1 = w1 - dw1
    
    # save current performance
    Error = backward$Error
    
    yk = forward$yk
    
    iter = iter + 1
    
    if (iter %% 1000 == 0) {
      cat(paste('iteration = ', iter, ' Error = ', Error, '\n'))
    }
  }
  
  if  (is.nan(Error)) {
    print(paste0('Error: ', Error))
  }
  
  # add prediction
  prediction = nnet_predict(test_set = training_set, w1 = w1, w2 = w2, w3 = w3, w4 = w4)
  
  return(list('y_k' = yk, 'Error' = Error, 'iterations' = iter, 'w1' = w1, 'w2' = w2, 'w3' = w3, 'w4' = w4, 'prediction' = prediction))
}