gaussian_kernel <- function(x1, x2, sigma) {
#GAUSSIANKERNEL returns a radial basis function kernel between x1 and x2
#   sim = gaussian_kernel(x1, x2, sigma) returns a gaussian kernel between x1 and x2
#   and returns the value in sim
#
# Converted to R by: SD Separa (2016/03/18)

	# Ensure that x1 and x2 are column vectors
	x1 = array(x1, c(length(x1), 1))
	x2 = array(x2, c(length(x2), 1))

	return(exp(-sum((x1 - x2)^2)/(2*sigma^2)))
}

linear_kernel <- function(x1, x2) {
#LINEARKERNEL returns a linear kernel between x1 and x2
#   sim = linear_kernel(x1, x2) returns a linear kernel between x1 and x2
#   and returns the value in sim
#
# Converted to R by: SD Separa (2016/03/18)

	# Ensure that x1 and x2 are column vectors
	x1 = array(x1, c(length(x1), 1))
	x2 = array(x2, c(length(x2), 1))

	# Compute the kernel
	return(t(x1) %*% x2)  # dot product
}

svm_train <- function(X, Y, C, kernelFunction, kernelParam, tol, max_passes) {
#SVMTRAIN Trains an SVM classifier using a simplified version of the SMO 
#algorithm. 
#   [model] = svm_train(X, Y, C, kernelFunction, tol, max_passes) trains an
#   SVM classifier and returns trained model. X is the matrix of training 
#   examples.  Each row is a training example, and the jth column holds the 
#   jth feature.  Y is a column matrix containing 1 for positive examples 
#   and 0 for negative examples.  C is the standard SVM regularization 
#   parameter.  tol is a tolerance value used for determining equality of 
#   floating point numbers. max_passes controls the number of iterations
#   over the dataset (without changes to alpha) before the algorithm quits.
#
# Note: This is a simplified version of the SMO algorithm for training
#       SVMs. In practice, if you want to train an SVM classifier, we
#       recommend using an optimized package such as:  
#
#           LIBSVM   (http://www.csie.ntu.edu.tw/~cjlin/libsvm/)
#           SVMLight (http://svmlight.joachims.org/)
#
# Converted to R by: SD Separa (2016/03/18)
	
	# for bsxfun and strcmp
	require(pracma)
	
	if (missing(tol) || is.null(tol)) {
		tol = 10^(-3)
	}

	if (missing(max_passes) || is.null(max_passes)) {
		max_passes = 5
	}

	# Data parameters
	m = nrow(X)
	n = ncol(X)

	# Map 0 to -1
	Y[which(Y == 0)] = -1

	# Variables
	alphas = array(0, c(m, 1))
	b = 0
	E = array(0, c(m, 1))
	passes = 0
	eta = 0
	L = 0
	H = 0

	kernelFunc = as.character(substitute(kernelFunction))
	
	# Pre-compute the Kernel Matrix since our dataset is small
	# (in practice, optimized SVM packages that handle large datasets
	#  gracefully will *not* do this)
	# 
	# We have implemented optimized vectorized version of the Kernels here so
	# that the svm training will run faster.
	if (strcmp(kernelFunc, 'linear_kernel')) {
		# Vectorized computation for the Linear Kernel
		# This is equivalent to computing the kernel on every pair of examples
		K = X %*% t(X)
	} else if (strcmp(kernelFunc, 'gaussian_kernel')) {
		# Vectorized RBF Kernel
		# This is equivalent to computing the kernel on every pair of examples
		X2 = as.matrix(rowSums(X^2))
		K = bsxfun('+', repmat(X2, 1, m), bsxfun('+', repmat(t(X2), m, 1), - 2 * X %*% t(X)))
		K = kernelFunction(1, 0, kernelParam)^K
	} else {
		# Pre-compute the Kernel Matrix
		# The following can be slow due to the lack of vectorization
		K = array(0, c(m, m))
		for (i in 1:m) {
			for (j in i:m) {
				 K[i,j] = kernelFunction(t(X[i,]), t(X[j,]), kernelParam)
				 K[j,i] = K[i,j] #the matrix is symmetric
			}
		}
	}
	
	# Train
	cat('Training ...')
	dots = 12

	while (passes < max_passes) {
				
		num_changed_alphas = 0
		
		for (i in 1:m) {
			
			# Calculate Ei = f(x(i)) - y(i) using (2). 
			# E(i) = b + sum (X(i, :) * (repmat(alphas.*Y,1,n).*X)') - Y(i);
			E[i] = b + sum(alphas*Y*K[,i]) - Y[i]
			
			if ((Y[i]*E[i] < -tol && alphas[i] < C) || (Y[i]*E[i] > tol && alphas[i] > 0)) {
				
				# In practice, there are many heuristics one can use to select
				# the i and j. In this simplified code, we select them randomly.
				j = ceiling(m * runif(1))
				
				while (j == i) {  # Make sure i \neq j
					j = ceiling(m * runif(1))
				}

				# Calculate Ej = f(x(j)) - y(j) using (2).
				E[j] = b + sum(alphas*Y*K[,j]) - Y[j]

				# Save old alphas
				alpha_i_old = alphas[i]
				alpha_j_old = alphas[j]
				
				# Compute L and H by (10) or (11). 
				if (Y[i] == Y[j]) {
					L = max(0, alphas[j] + alphas[i] - C)
					H = min(C, alphas[j] + alphas[i])
				} else {
					L = max(0, alphas[j] - alphas[i])
					H = min(C, C + alphas[j] - alphas[i])
				}
			   
				if (L == H) {
					# continue to next i. 
					next
				}

				# Compute eta by (14).
				eta = 2 * K[i,j] - K[i,i] - K[j,j]
				if (eta >= 0) {
					# continue to next i. 
					next
				}
				
				# Compute and clip new value for alpha j using (12) and (15).
				alphas[j] = alphas[j] - (Y[j] * (E[i] - E[j])) / eta
				
				# Clip
				alphas[j] = min (H, alphas[j])
				alphas[j] = max (L, alphas[j])
				
				# Check if change in alpha is significant
				if (abs(alphas[j] - alpha_j_old) < tol) {
					# continue to next i. 
					# replace anyway
					alphas[j] = alpha_j_old
					next
				}
				
				# Determine value for alpha i using (16). 
				alphas[i] = alphas[i] + Y[i]*Y[j]*(alpha_j_old - alphas[j])
				
				# Compute b1 and b2 using (17) and (18) respectively. 
				b1 = b - E[i] - Y[i] * (alphas[i] - alpha_i_old) %*%  t(K[i,j]) - Y[j] * (alphas[j] - alpha_j_old) %*%  t(K[i,j])
				b2 = b - E[j] - Y[i] * (alphas[i] - alpha_i_old) %*%  t(K[i,j]) - Y[j] * (alphas[j] - alpha_j_old) %*%  t(K[j,j])

				# Compute b by (19). 
				if (0 < alphas[i] && alphas[i] < C) {
					b = b1
				} else if (0 < alphas[j] && alphas[j] < C) {
					b = b2
				} else {
					b = (b1+b2)/2
				}

				num_changed_alphas = num_changed_alphas + 1
			}
		}
		
		if (num_changed_alphas == 0) {
			passes = passes + 1
		} else {
			passes = 0
		}

		cat('.')
		dots = dots + 1
		if (dots > 78) {
			dots = 0
			cat('\n')
		}
	}

	cat(' Done!\n')

	# Save the model
	idx = alphas > 0
	model<-list()
	model$X = X[idx,]
	model$y = Y[idx]
	model$kernelFunction = kernelFunc
	model$kernelParam = kernelParam
	model$b = b
	model$alphas=alphas[idx]
	model$w = t(t(alphas*Y)%*%X)

	return(model)
}

svm_predict <- function(model, X) {
#SVMPREDICT returns a vector of predictions using a trained SVM model
#(svm_train). 
#   pred = SVMPREDICT(model, X) returns a vector of predictions using a 
#   trained SVM model (svm_train). X is a mxn matrix where there each 
#   example is a row. model is a svm model returned from svm_train.
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
	
	if (strcmp(model$kernelFunction,'linear_kernel')) {
		# We can use the weights and bias directly if working with the 
		# linear kernel
		p = X %*% model$w + model$b
	} else if (strcmp(model$kernelFunction, 'gaussian_kernel')) {
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

svm_boundary <- function(X, y, model) {
#SVMBOUNDARY plots a non-linear decision boundary learned by the SVM
#   svm_boundary(X, y, model) plots a non-linear decision 
#   boundary learned by the SVM and overlays the data on it
#
# Converted to R by: SD Separa (2016/03/18)

	# for meshgrid
	require(pracma)
	
	# Plot the training data on top of the boundary
	svm_plot(X, y)

	# Make classification predictions over a grid of values
	x1plot = t(seq(min(X[,1]), max(X[,1]),length=100))
	x2plot = t(seq(min(X[,2]), max(X[,2]),length=100))

	xgrid = meshgrid(x1plot, x2plot)
	X1 = xgrid$X
	X2 = xgrid$Y

	vals = array(0, dim(X1))

	for (i in 1:ncol(X1)) {
	   this_X = cbind(X1[,i], X2[,i])
	   vals[, i] = svm_predict(model, this_X)
	}

	contour(x = x1plot, y = x2plot, z = t(vals), col = 'green', add = TRUE, lw = 1)
}

svm_plot <- function(X, y) {
#SVMPLOT Plots the data points X and y into a new figure 
#   svm_plot(x,y) plots the data points with + for the positive examples
#   and o for the negative examples. X is assumed to be a Mx2 matrix.
#
# Note: This was slightly modified such that it expects y = 1 or y = 0
#
# Converted to R by: SD Separa (2016/03/18)

	# Find Indices of Positive and Negative Examples
	pos = which(y == 1)
	neg = which(y == 0)

	# Plot Examples
	plot(X[pos, 1], X[pos, 2], pch = 3, xlab = '', ylab = '', col = 'blue')
	points(X[neg, 1], X[neg, 2], pch = 19, col = 'red')
}
