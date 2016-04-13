svcplot <- function(X, Y, ker, kpar1, kpar2, alpha, bias, aspect, mag, xaxis, yaxis, input) {

##########################################################################
# FUNCTION
#  svcplot_book(X,Y,ker,kpar1,kpar2,alpha,bias,aspect,mag,xaxis,yaxis,input)
# Support Vector Machine Plotting routine. It plots the decision regions, the decision
# surfaces and the margin obtained by a SVM classifier.
#
# INPUT ARGUMENTS:
#  X:       training inputs
#  Y:       training targets
#  ker:     kernel function
#  kpar1:   1st parameter of kernel
#  kpar2:   2nd parameter of kernel
#  alpha:   Lagrange Multipliers
#  bias:    bias term
#  aspect:  aspect Ratio (default: 0 (fixed), 1 (variable))
#  mag:     display magnification
#  xaxis:   xaxis input (default: 1)
#  yaxis:   yaxis input (default: 2)
#  input:   vector of input values (default: zeros(no_of_inputs))
#
#  Original Author: Steve Gunn (srg@ecs.soton.ac.uk)
#  Modified by Michael Mavroforakis
#
# Converted from Matlab into R by S.D. Separa (04/12/2016)
#
##########################################################################

	# for meshgrid
	require(pracma)
	
	color_shade = 1 # 1:color in figure and shade degradation, else 0:Black and White
	gridcellsX = 50 # num of grid cells in X-dimension (test:20, presentation:60-80)
	gridcellsY = 50 # num of grid cells in Y-dimension (test:20, presentation:60-80)
	marg = 0.1 # percent of margin around the border points

	if (nargs() < 7 | nargs() > 12) { # check correct number of arguments
		stop('incorrect number of arguments')
	} else {
		epsilon = 10^(-5)
		
		if (nargs() < 12) input = array(0, c(1, ncol(X)))
		if (nargs() < 11) yaxis = 2
		if (nargs() < 10) xaxis = 1
		if (nargs() < 9) mag = 0.1
		if (nargs() < 8) aspect = 0
		
		# Calculate values to Scale the axes
		xmin = min(X[,xaxis])
		xmax = max(X[,xaxis])
		ymin = min(X[,yaxis])
		ymax = max(X[,yaxis])
		xa = (xmax - xmin)
		ya = (ymax - ymin)
		
		if (!aspect) {
			if (0.75 * abs(xa) < abs(ya)) {
				offadd = marg * (ya * 4 / 3 - xa)
				xmin = xmin - offadd - mag * marg * ya
				xmax = xmax + offadd + mag * marg *ya
				ymin = ymin - mag * marg * ya
				ymax = ymax + mag * marg *ya
			} else {
				offadd = marg * (xa * 3/4 - ya)
				xmin = xmin - mag * marg *xa
				xmax = xmax + mag * marg *xa
				ymin = ymin - offadd - mag * marg * xa
				ymax = ymax + offadd + mag * marg * xa
			}
		} else {
			xmin = xmin - mag*marg*xa
			xmax = xmax + mag*marg*xa
			ymin = ymin - mag*marg*ya
			ymax = ymax + mag*marg*ya
		}
		
		alpha_min = min(alpha)
		alpha_max = max(alpha)
		alpha_threshold = (alpha_max - alpha_min) * 0.01
		alpha_threshold = alpha_threshold + alpha_min
		
		# Plot function value
		x = seq(xmin, xmax, length = gridcellsX)
		y = seq(ymin, ymax, length = gridcellsY)
		z = bias * array(1, c(length(y), length(x)))
		
		for (x1 in 1:length(x)) {
			for (y1 in 1:length(y)) {

				input[xaxis] = x[x1]
				input[yaxis] = y[y1]
				
				for (i in 1:length(Y)) {
					if (abs(alpha[i]) >= 0) {
						z[y1, x1] = z[y1, x1] + Y[i] * alpha[i] * CalcKernel(input, X[i, ], ker, kpar1, kpar2)
					}
				}
			}
		}
		
		l = (-min(z) + max(z)) / 2.0
		
		plot.new()
		
		#Plot Training points
		pos = which(Y == 1)
		neg = which(Y != 1)
		alp = which(alpha > alpha_threshold)
		
		z = t(z)
		
		if (color_shade == 1) {
			points(x = X[pos, xaxis], y = X[pos, yaxis], xlim = c(xmin, xmax), ylim = c(ymin, ymax), col = 'red', pch = 4)
			points(x = X[neg, xaxis], y = X[neg, yaxis], xlim = c(xmin, xmax), ylim = c(ymin, ymax), col = 'blue', pch = 4)
			points(x = X[alp, xaxis], y = X[alp, yaxis], xlim = c(xmin, xmax), ylim = c(ymin, ymax), col = 'black', pch = 1)
		} else {
			points(x = X[, xaxis], y = X[, yaxis], xlim = c(xmin, xmax), ylim = c(ymin, ymax), col = 'black', pch = 4)
			points(x = X[alp, xaxis], y = X[alp, yaxis], xlim = c(xmin, xmax), ylim = c(ymin, ymax), col = 'black', pch = 19)
		}
		
		# Plot Boundary contours
		if (color_shade == 1) {
			contour(x, y, z, levels = c(0, 0), col = 'black', add = TRUE, lty = 1, drawlabels = FALSE)
			contour(x, y, z, levels = c(-1, -1), col = 'blue', add = TRUE, lty = 3, drawlabels = FALSE)
			contour(x, y, z, levels = c(1, 1), col = 'red', add = TRUE, lty = 3, drawlabels = FALSE)
		} else {
			zones = 1 # how many zones to be present in [0,1]
			steps = seq(0, 1, by = 1 / zones)
			for (j in 1:(zones + 1)) {
				if (j %% 2 == 1) {
					clsp = 1
				} else if (j %% 2 == 0) {
					clsp = 3
				}
				if (steps[j] == 0) {
					contour(x, y, z, levels = c(steps[j], steps[j]), col = 'black', lty = clsp, lwd = 2, add = TRUE, drawlabels = FALSE)
				} else {
					contour(x, y, z, levels = c(steps[j], steps[j]), col = 'black', lty = clsp, lwd = 1, add = TRUE, drawlabels = FALSE)
					contour(x, y, z, levels = c(-steps[j], -steps[j]), col = 'black', lty = clsp, lwd = 1, add = TRUE, drawlabels = FALSE)
				}
			}
		}
	}
}
