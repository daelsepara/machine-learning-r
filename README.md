# Machine Learning in R

This is a repository for R scripts developed during my machine learning studies. Some of the codes have been adapted and converted into R from their original Matlab implementations. 

## Neural Network ##

- sigmoid activation
- training, forward and backward propagation (nnet_train, nnet_forward, nnet_backprop)
- multi-class classification/prediction (nnet_predict)
- fast gradient descent computation (nnet_optimize, nnet_cost, fmincg)

## Support Vector Machines ##

- training, prediction (svmTrain, svmPredict)
- decision boundary visualization (plotData, visualizeBoundary)
- linear and gaussian kernels (linearKernel, gaussianKernel)

## Features ##

- normalize all features of X is an (m x n) matrix: m-examples with n-features (featureNormalize)

## Principal Component Analysis ##

- compute principal components (pca)
- create projections of X (onto U) using K principal components (pca_project)
- create approximations of X using K principal components (pca_estimate)

## Anomaly Detection ##

- estimate Guassian distribution parameters: mean, variance (estimateGaussian)
- compute multivariate Gaussian probability distribution function (multivariateGaussian)
- determine best prediction threshold to use (selectThreshold)

## Collaborative Filtering ##

- compute cost function and gradients with Regularization (cf_costfunction )
- optimize using fast gradient descent optimizer fmincg (cf_optimize, cf_fmincg_cost)

## Clustering ##

### k-means clustering ###

- initialization of K-centroids (kmeans_initialize)
- cluster assignment to nearest centroid (kmeans_assign)
- compute new centroids from all points belonging to K-clusters (kmeans_compute)
- run K-means algorithm for a maximum number of iterations (kmeans_run)
