plotData <- function(X, y) {
#PLOTDATA Plots the data points X and y into a new figure 
#   PLOTDATA(x,y) plots the data points with + for the positive examples
#   and o for the negative examples. X is assumed to be a Mx2 matrix.
#
# Note: This was slightly modified such that it expects y = 1 or y = 0
#
# Converted to R by: SD Separa (2016/03/18)

	# Find Indices of Positive and Negative Examples
	pos = which(y == 1)
	neg = which(y == 0)

	# Plot Examples
	plot(X[pos, 1], X[pos, 2], pch = 3, xlab = '', ylab = '', col = 'black')
	points(X[neg, 1], X[neg, 2], pch = 19, col = 'yellow')
}
