# Machine Learning in R

This is a repository for R scripts developed during my machine learning studies. Some of the codes have been adapted and converted into R from their original Matlab implementations. 

## Classification ##

- Euclidean (euclidean_classifier)
- Mahalanobis (mahalanobis_classifier)
- Perceptron (perceptron_classifier)
- Online Perceptron (online_perceptron_classifier)
- Sum-Squared Error (sse_classifier)

## Regression ##

- plot data and  (regression_plot)
- plot regression decision boundary (regression_boundary)
- generic regression wrapper function (regression_optimize)

### Linear regression ###

- Linear regression cost function and gradient (lr_cost)
- Linear regression gradient descent (lr_gradientdescent)

### Logistic regression ###

- Logistic regression cost function and gradient (logr_cost)
- Logistic regression optimizer (logr_optimize)
- Predictions (logr_predict)

### Softmax regression ###

- Softmax regression cost function and gradient (softmax_cost)
- Softmax regression optimizer (softmax_optimize)
- Predictions (softmax_predict)

### Softplus regression ###

- Softplus regression cost function and gradient (softplus_cost)
- Softplus regression optimizer (softplus_optimize)
- Predictions (softplus_predict)

## Neural Network ##

- sigmoid activation and derivative (nnet_sigmoid, nnet_dsigmoid)
- softmax activation for multi-class classification/prediction (nnet_softmax)
- forward and backward propagation (nnet_forward, nnet_backprop)
- training via gradient descent (nnet_train)
- training using stochastic gradient descent without regularization (nnet_stochastic)
- multi-class classification/prediction (nnet_predict)
- cost function for use with minimization algorithms (nnet_cost)
- fast gradient descent computation (fmincg)
- training using fast gradient descent computation (nnet_optimize)
- training using R's optim function (nnet_minimize)

## Support Vector Machines ##

- training, prediction (svm_train, svm_predict)
- decision boundary visualization (svm_plot, svm_boundary)
- linear, gaussian/rbf, polynomial kernels (linear_kernel, gaussian_kernel, polynomial_kernel)
- sequential minimal optimization using various kernel functions (smo2, CalcKernel)
- sequential minimal optimization boundary plotter (svcplot)

### Supported Kernels in SMO (CalcKernel) ###

- linear
- polynomial (poly)
- radial basis function (rbf)
- extended radial basis function (erbf)
- sigmoid
- fourier
- spline
- curved spline, ANOVA (curvedspline, anova)
- B-spline (bspline)
- ANOVA spline 1,2,3 (anovaspline1, anovaspline2, anovaspline3)
- ANOVA B-spline (anovabspline)

## Features ##

- normalize all features of X[m, n]: m-examples with n-features (featureNormalize)

## Principal Component Analysis ##

- compute principal components (pca)
- create projections of X (onto U) using K principal components (pca_project)
- create approximations of X using K principal components (pca_estimate)

## Anomaly Detection ##

- estimate Guassian distribution parameters: mean, variance (estimate_gaussian)
- compute multivariate Gaussian probability distribution function (multivariate_gaussian)
- determine best prediction threshold to use (select_threshold)

## Collaborative Filtering ##

- compute cost function and gradients with Regularization (cf_costfunction )
- optimize using fast gradient descent optimizer fmincg (cf_optimize, cf_fmincg_cost)

## Clustering ##

- basic sequential algorithmic scheme (BSAS)
- generalized agglomerative scheme with single and complete linking (agglom)
- run BSAS, GAS tests (test_clustering)

### k-means clustering ###

- initialization of K-centroids (kmeans_initialize)
- cluster assignment to nearest centroid (kmeans_assign)
- compute new centroids from all points belonging to K-clusters (kmeans_compute)
- run K-means algorithm for a maximum number of iterations (kmeans_run)
