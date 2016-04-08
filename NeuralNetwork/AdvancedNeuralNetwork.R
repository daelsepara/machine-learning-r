sigmoid <- function(x) {
# Sigmoid activation function
  
  return(1/(1 + exp(-x)))
}

dsigmoid <- function(x) {
# 1st-derivative of sigmoid activation function
  
  z = sigmoid(x)
  return(z * (1 - z))
}

nnet_forward <- function(training_set, w_ji, w_kj, softmax = FALSE) {
# Forward propagation
  
  # add bias column to input layer
  x = cbind(array(1, c(nrow(training_set), 1)), training_set)
  
  # compute hidden layer activation
  z_2 = x %*% t(w_ji)
  z_j = h_func(z_2)
  
  # add bias column
  a_2 = cbind(array(1, c(nrow(z_j), 1)), z_j)
  
  # compute output layer
  if (!softmax) {
    y_k = sigmoid(a_2 %*% t(w_kj))
  } else {
    y_k = nnet_softmax(a_2 %*% t(w_kj))
  }
  
  return(list('y_k' = y_k, 'z_2' = z_2, 'a_2' = a_2))
}

nnet_backprop <- function(training_set, y_k, z_2, a_2, w_ji, w_kj, y_matrix, lambda = 0, softmax = FALSE) {
# Backward propagation
  
  # add bias column
  x = cbind(array(1, c(nrow(training_set), 1)), training_set)
  m = nrow(x)
  
  # compute intermediate delta values per layer
  d3 = y_k - y_matrix
  d2 = d3 %*% w_kj[, 2:ncol(w_kj)] * dsigmoid(z_2)
  
  # compute cost function and gradient
  if (!softmax) {
    
    dWji = (t(d2) %*% x) / m
    dWkj = (t(d3) %*% a_2) / m
    cost = sum(-y_matrix * log(y_k) - (1 - y_matrix) * log(1 - y_k)) / m
  	
  } else {
    
    # softmax activation cost function
    cost = - sum(log(y_k[which(y_matrix == 1)]))
    dWji = (t(d2) %*% x)
    dWkj = (t(d3) %*% a_2)
  }
  
  # regularization on lambda != 0
  if (lambda != 0) {
    
    rWji = w_ji
    rWkj = w_kj
    
    # do not regularize bias column
    rWji[, 1] = array(0, nrow(w_ji))
    rWkj[, 1] = array(0, nrow(w_kj))
    
    if (!softmax)	{
      cost = cost + lambda * (sum(rWji ^ 2) + sum(rWkj ^ 2)) / (2 * m)
      dWji = dWji + lambda * rWji / m
      dWkj = dWkj + lambda * rWkj / m
    } else {
      cost = cost + lambda * (sum(rWji ^ 2) + sum(rWkj ^ 2)) / 2
      dWji = dWji + lambda * rWji
      dWkj = dWkj + lambda * rWkj
    }
  }
  
  return(list('dWkj' = dWkj, 'dWji' = dWji, 'Error' = cost))	
}

nnet_cost <- function(X, P1, P2, P3 , P4, P5, P6, softmax = FALSE) {
# Neutral network cost function for use with advanced optimization method (fmincg/optim)
  
  # P1 training_set
  # P2 y_matrix (expected output)
  # P3 number of input units
  # P4 hidden layer units
  # P5 number of labels
  # P6 lambda (regularization)
  
  # roll up vectors into arrays
  offs = P4 * (P3 + 1)
  w_ji = array(X[1:offs], c(P4, P3 + 1))
  w_kj = array(X[(1 + offs):length(X)], c(P5, P4 + 1))
  
  # compute cost function (J) and its gradients (partial derivatives)
  # using forward and backpropagation
  forward = nnet_forward(P1, w_ji, w_kj, softmax)
  result = nnet_backprop(P1, forward$y_k, forward$z_2, forward$a_2, w_ji, w_kj, P2, P6, softmax)
  
  # unroll gradient matrices into one vector
  grad = c(as.vector(result$dWji), as.vector(result$dWkj))
  
  return(list('J' = result$Error, 'grad' = grad))
}

nnet_labels <- function(output, num_labels) {
# create labels for multi-class classification
  
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

nnet_train <- function(maxiter = 100, learning_rate = 0.1, tol = 10^(-3), training_set = array(0) , output = array(0), hidden_units = 0, num_labels = 1, min_max = 1, isGaussian = FALSE, softmax = FALSE) {
# Network training
  
  y_matrix = nnet_labels(output, num_labels)
  
  # determine network dimensions from user input
  j = hidden_units
  inputs = ncol(training_set)
  
  # initialize weights with random values
  w_ji = nnet_weights(min_max, j, inputs + 1, isGaussian)
  w_kj = nnet_weights(min_max, num_labels, j + 1, isGaussian)
  
  iter = 0
  Error = 1.0
  
  y_k = numeric(0)
	
  while (iter < maxiter && Error > tol) {
    # for training, perform forward and backpropagation each iteration, no regularization
    forward = nnet_forward(training_set, w_ji, w_kj, softmax)
    backward = nnet_backprop(training_set, forward$y_k, forward$z_2, forward$a_2, w_ji, w_kj, y_matrix, 0, softmax)
    
    # update weights (using learning rate and gradient descent)
    w_ji = w_ji - learning_rate * backward$dWji
    w_kj = w_kj - learning_rate * backward$dWkj
  	
    # save current performance
    Error = backward$Error
    y_k = forward$y_k
    
    iter = iter + 1
    
    if (iter %% 1000 == 0) {
      cat(paste('iteration = ', iter, ' Error = ', Error, '\n'))
    }
  }
  
  # add prediction
  prediction = nnet_predict(training_set, w_ji, w_kj)
  
  return(list('y_k' = y_k, 'Error' = Error, 'iterations' = iter, 'w_kj' = w_kj, 'w_ji' = w_ji, 'prediction' = prediction))
}

nnet_weights <- function(min_max = 1, m = 1, n = 1, isGaussian = FALSE) {
# intialize interconnection weights with random values (-min_max, min_max) or Gaussian (mean = 0, sd = min_max)
  
  if (!isGaussian) {
    
    return(array(runif(n = m * n, min = -min_max, max = min_max), c(m, n)))
    
  } else {
    
    return(array(rnorm(n = m * n, mean = 0, sd = abs(min_max)), c(m, n)))
  }  
}

nnet_optimize <- function(maxiter = 100, training_set = array(0) , output = array(0), hidden_units = 0, num_labels = 1, min_max = 1, isGaussian = FALSE, lambda = 0) {
# Network training using advanced optimization algorithm fmincg
  
  y_matrix = nnet_labels(output, num_labels)
	
  # determine network dimensions from user input
  j = hidden_units
  inputs = ncol(training_set)
  
  # initialize weights with random values
  w_ji = nnet_weights(min_max, j, inputs + 1, isGaussian)
  w_kj = nnet_weights(min_max, num_labels, j + 1, isGaussian)
  	
  theta = c(as.vector(w_ji), as.vector(w_kj))
  result = fmincg(nnet_cost, theta, maxiter, training_set, y_matrix, inputs, j, num_labels, lambda)
  
  offs = j * (inputs + 1)
  w_ji = array(result$X[1:offs], c(j, inputs + 1))
  w_kj = array(result$X[(1 + offs):length(result$X)], c(num_labels, j + 1))
  	
  # performance
  Error = result$cost
  y_k = nnet_forward(training_set, w_ji, w_kj)$y_k
    
  # add prediction
  prediction = nnet_predict(training_set, w_ji, w_kj)
  
  return(list('y_k' = y_k, 'Error' = Error, 'iterations' = result$i, 'w_kj' = w_kj, 'w_ji' = w_ji, 'prediction' = prediction))
}

nnet_minimize <- function(maxiter = 100, training_set = array(0) , output = array(0), hidden_units = 0, num_labels = 1, min_max = 1, isGaussian = FALSE, lambda = 0, method = 'L-BFGS-B', softmax = FALSE) {
# Network training using R's optimizer
  
  y_matrix = nnet_labels(output, num_labels)
  
  # determine network dimensions from user input
  j = hidden_units
  inputs = ncol(training_set)
  
  # initialize weights with random values
  w_ji = nnet_weights(min_max, j, inputs + 1, isGaussian)
  w_kj = nnet_weights(min_max, num_labels, j + 1, isGaussian)
  
  theta = c(as.vector(w_ji), as.vector(w_kj))
  
  # optim works with functions with one argument/parameter. We define anonymous functions (which are just wrappers to our cost function) to acheive the desired effect
  result = optim(par = theta, fn = function(theta) { return(nnet_cost(theta, training_set, y_matrix, inputs, j, num_labels, lambda, softmax)$J) }, gr = function(theta) { return(nnet_cost(theta, training_set, y_matrix, inputs, j, num_labels, lambda, softmax)$grad) }, control = list('maxit' = maxiter), method = method)
  
  offs = j * (inputs + 1)
  w_ji = array(result$par[1:offs], c(j, inputs + 1))
  w_kj = array(result$par[(1 + offs):length(result$par)], c(num_labels, j + 1))
  
  # performance
  Error = result$value
  y_k = nnet_forward(training_set, w_ji, w_kj, softmax)$y_k
  
  # add prediction
  prediction = nnet_predict(test_set = training_set, w_ji = w_ji, w_kj = w_kj, softmax = softmax)
  
  return(list('y_k' = y_k, 'fn' = result$counts[1], 'gr' = result$counts[2], 'Error' = Error, 'w_kj' = w_kj, 'w_ji' = w_ji, 'prediction' = prediction))
}

nnet_softmax <- function(x) {
# numerically stable softmax activation function for use with multi-class classification
  
  # for repmat
  require(pracma)
  
  k = ncol(x)
  m = nrow(x)
  
  # get the maximum in each row
  m_k = array(apply(x, 1, max), c(m, 1))
  ezk = exp(x - repmat(m_k, 1, k))

  # compute softmax
  smk = ezk / repmat(array(apply(ezk, 1, sum), c(m, 1)), 1, k)
  
  return(smk)
}

nnet_predict <- function(test_set, w_ji, w_kj, threshold = 0.5, softmax = FALSE) {
# Predict using neural network parameters (multi-class classification)
  
  prediction_output = nnet_forward(test_set, w_ji, w_kj, softmax)$y_k
  
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
