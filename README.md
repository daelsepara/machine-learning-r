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

## Feature Normalization ##

- normalize all features of X is an (m x n) matrix: m-examples with n-features (featureNormalize)
- principal component analysis (pca)
- project and approximate X using K-components of PCA (pca_project, pca_estimate)

## K-means clustering ##

- initialization of K-centroids (kmeans_initialize)
- cluster assignment to nearest centroid (kmeans_assign)
- compute new centroids from all points belonging to K-clusters (kmeans_compute)
- run K-means algorithm for a maximum number of iterations (kmeans_run)
