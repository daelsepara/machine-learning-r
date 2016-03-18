linearKernel <- function(x1, x2) {
#LINEARKERNEL returns a linear kernel between x1 and x2
#   sim = linearKernel(x1, x2) returns a linear kernel between x1 and x2
#   and returns the value in sim
#
# Converted to R by: SD Separa (2016/03/18)

	# Ensure that x1 and x2 are column vectors
	x1 = array(x1, c(length(x1), 1))
	x2 = array(x2, c(length(x2), 1))

	# Compute the kernel
	return(t(x1) %*% x2)  # dot product
}
