svmPredict <- function(model, X) {
#SVMPREDICT returns a vector of predictions using a trained SVM model
#(svmTrain). 
#   pred = SVMPREDICT(model, X) returns a vector of predictions using a 
#   trained SVM model (svmTrain). X is a mxn matrix where there each 
#   example is a row. model is a svm model returned from svmTrain.
#   predictions pred is a m x 1 column of predictions of {0, 1} values.
#
# Converted to R by: SD Separa (2016/03/18)

	# for bsxfun and strcmp
	require(pracma)
	
	# Check if we are getting a column vector, if so, then assume that we only
	# need to do prediction for a single example
	if (ncol(X) == 1) {
		# Examples should be in rows
		X = t(X)
	}

	# Dataset 
	m = nrow(X)
	p = array(0, c(m, 1))
	pred = array(0, c(m, 1))

	# [sdsepara] map model kernel function to an actual function
	kernelFunction <- match.fun(model$kernelFunction)
	
	if (strcmp(model$kernelFunction,'linearKernel')) {
		# We can use the weights and bias directly if working with the 
		# linear kernel
		p = X %*% model$w + model$b
	} else if (strcmp(model$kernelFunction, 'gaussianKernel')) {
		# Vectorized RBF Kernel
		# This is equivalent to computing the kernel on every pair of examples
		X1 = as.matrix(rowSums(X^2))
		X2 = t(as.matrix(rowSums(model$X^2)))
		rows = nrow(X1)
		K = bsxfun('+', repmat(X1, 1, ncol(X2)), bsxfun('+', repmat(X2, rows, 1), -2*X%*%t(model$X)))
		K = kernelFunction(1, 0, model$kernelParam)^K
		K = bsxfun('*', repmat(t(model$y), rows, 1), K)
		K = bsxfun('*', repmat(t(model$alphas), rows, 1), K)
		p = rowSums(K)
	} else {
		# Other Non-linear kernel
		for (i in 1:m) {
			prediction = 0
			for (j in 1:nrow(model$X)) {
				prediction = prediction + model$alphas[j] * model$y[j] * kernelFunction(t(X[i,]), t(model$X[j,]), model$kernelParam)
			}
			p[i] = prediction + model$b
		}
	}

	# Convert predictions into 0 / 1
	pred[p >= 0] =  1
	pred[p <  0] =  0
	
	return(pred)
}

