visualizeBoundary <- function(X, y, model) {
#VISUALIZEBOUNDARY plots a non-linear decision boundary learned by the SVM
#   VISUALIZEBOUNDARYLINEAR(X, y, model) plots a non-linear decision 
#   boundary learned by the SVM and overlays the data on it
#
# Converted to R by: SD Separa (2016/03/18)

	# Plot the training data on top of the boundary
	plotData(X, y)

	# Make classification predictions over a grid of values
	x1plot = t(seq(min(X[,1]), max(X[,1]),length=100))
	x2plot = t(seq(min(X[,2]), max(X[,2]),length=100))

	xgrid = meshgrid(x1plot, x2plot)
	X1 = xgrid$X
	X2 = xgrid$Y

	vals = array(0, dim(X1))

	for (i in 1:ncol(X1)) {
	   this_X = cbind(X1[,i], X2[,i])
	   vals[, i] = svmPredict(model, this_X)
	}

	contour(x = x1plot, y = x2plot, z = t(vals), col = 'green', add = TRUE, lw = 1)
}
