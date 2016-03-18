gaussianKernel <- function(x1, x2, sigma) {

#RBFKERNEL returns a radial basis function kernel between x1 and x2
#   sim = gaussianKernel(x1, x2) returns a gaussian kernel between x1 and x2
#   and returns the value in sim

	# Ensure that x1 and x2 are column vectors
	x1 = array(x1, c(length(x1), 1))
	x2 = array(x2, c(length(x2), 1))

	return(exp(-sum((x1 - x2)^2)/(2*sigma^2)))
}
