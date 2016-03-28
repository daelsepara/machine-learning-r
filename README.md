# Machine Learning in R

This is a repository for R scripts developed during my machine learning studies. Some of the codes have been adapted and converted into R from their original Matlab implementations. 

## Classification ##

- Euclidean classifier (euclidean_classifier)

## Regression ##

- plot data and decision boundary (regression_plot, regression_boundary)

### Linear regression ###

- Linear regression cost function and gradient (lr_cost)
- Linear regression gradient descent (lr_gradientdescent)

### Logistic regression ###

- Logistic regression cost function and gradient (logr_cost)
- Logistic regression optimizer (logr_optimize)

## Neural Network ##

- sigmoid activation
- training, forward and backward propagation (nnet_train, nnet_forward, nnet_backprop)
- multi-class classification/prediction (nnet_predict)
- fast gradient descent computation (nnet_optimize, nnet_cost, fmincg)

## Support Vector Machines ##

- training, prediction (svm_train, svm_predict)
- decision boundary visualization (svm_plot, svm_boundary)
- linear and gaussian kernels (linear_kernel, gaussian_kernel)

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

### k-means clustering ###

- initialization of K-centroids (kmeans_initialize)
- cluster assignment to nearest centroid (kmeans_assign)
- compute new centroids from all points belonging to K-clusters (kmeans_compute)
- run K-means algorithm for a maximum number of iterations (kmeans_run)
