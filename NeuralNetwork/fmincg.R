fmincg<-function(f, X, options_ = 100, P1, P2, P3, P4, P5, P6) {
# Minimize a continuous differentialble multivariate function. Starting point
# is given by "X" (D by 1), and the function named in the string "f", must
# return a function value and a vector of partial derivatives. The Polack-
# Ribiere flavour of conjugate gradients is used to compute search directions,
# and a line search using quadratic and cubic polynomial approximations and the
# Wolfe-Powell stopping criteria is used together with the slope ratio method
# for guessing initial step sizes. Additionally a bunch of checks are made to
# make sure that exploration is taking place and that extrapolation will not
# be unboundedly large. The "length" gives the length of the run: if it is
# positive, it gives the maximum number of line searches, if negative its
# absolute gives the maximum allowed number of function evaluations. You can
# (optionally) give "length" a second component, which will indicate the
# reduction in function value to be expected in the first line-search (defaults
# to 1.0). The function returns when either its length is up, or if no further
# progress can be made (ie, we are at a minimum, or so close that due to
# numerical problems, we cannot get any closer). If the function terminates
# within a few iterations, it could be an indication that the function value
# and derivatives are not consistent (ie, there may be a bug in the
# implementation of your "f" function). The function returns the found
# solution "X", a vector of function values "fX" indicating the progress made
# and "i" the number of iterations (line searches or function evaluations,
# depending on the sign of "length") used.
#
# Usage: [X, fX, i] = fmincg(f, X, options, P1, P2, P3, P4, P5, P6)
#
# See also: checkgrad 
#
# Copyright (C) 2001 and 2002 by Carl Edward Rasmussen. Date 2002-02-13
#
#
# (C) Copyright 1999, 2000 & 2001, Carl Edward Rasmussen
# 
# Permission is granted for anyone to copy, use, or modify these
# programs and accompanying documents for purposes of research or
# education, provided this copyright notice is retained, and note is
# made of any changes that have been made.
# 
# These programs and documents are distributed without any warranty,
# express or implied.  As the programs were written for research
# purposes only, they have not been tested to the degree that would be
# advisable in any important application.  All use of these programs is
# entirely at the user's own risk.
#
# [ml-class] Changes Made:
# 1) Function name and argument specifications
# 2) Output display
#
# This conversion to R by:  SD Separa (2016)
#
# [sdsepara] Changes made:
#
# 1) suppress warnings on numerical error
# 2) realmin obtained from machine parameters
# 3) add one more check if some operations are NAN; if() fails
# 4) calls f() directly passing all parameters P1-P6 instead of creating an argstr

	realmin = .Machine$double.xmin

	# if check if reduction parameter was specified
	if (length(options_) > 1) {
		runLength = options_[1]
		red = options_[2]
	} else {
		runLength = options_
		red = 1
	}
	
	S='Iteration '

	RHO = 0.01	# a bunch of constants for line searches
	SIG = 0.5		# RHO and SIG are the constants in the Wolfe-Powell conditions
	INT = 0.1		# don't reevaluate within 0.1 of the limit of the current bracket
	EXT = 3.0		# extrapolate maximum 3 times the current bracket
	MAX = 20		# max 20 function evaluations per line search
	RATIO = 100	# maximum allowed slope ratio

	i = 0 			# zero the run length counter
	ls_failed = 0	# no previous line search has failed
	fX = numeric(0)

	# get function value and gradient
	eval_ = f(X, P1, P2, P3 , P4, P5, P6)
	f1 = eval_$J
	df1 = eval_$grad

	# count epochs?!
	i = i + (runLength < 0)
	
	s = -df1          # search direction is steepest
	d1 = t(-s) %*% s  # this is the slope
	z1 = red/(1 - d1) # initial step is red/(|s|+1)

	while (i < abs(runLength)) {
	  
		# count iterations?!
		i = i + (runLength > 0)

		# make a copy of current values
		X0 = X
		f0 = f1
		df0 = df1

		# begin line search
		X = X + z1 %*% s
	  
		eval_ = f(X, P1, P2, P3 , P4, P5, P6)                                 
		f2 = eval_$J
		df2 = eval_$grad
		
		# count epochs?!
		i = i + (runLength < 0)
		
		# initialize point 3 equal to point 1
		d2 = t(df2) %*% s
		f3 = f1
		d3 = d1
		z3 = -z1
		
		if (runLength>0) {
			M = MAX
		} else {
			M = min(MAX, -runLength - i)
		}
		
		# initialize quantities
		success = 0
		limit = -1
		
		while(1) {
			while (((f2 > f1 + z1 * RHO * d1) || (d2 > -SIG * d1)) && (M > 0)) {

				# tighten the bracket
				limit = z1
				
				if (f2 > f1) {
					# quadratic fit
					z2 = z3 - (0.5 * d3 * z3 * z3) / (d3 * z3 + f2 - f3)
				} else {
					# cubic fit
					A = 6 * (f2 - f3) / z3 +3 *(d2 + d3)
					B = 3 * (f3 - f2) - z3 * (d3 + 2*d2)
					
					# numerical error possible - ok!
					z2 = suppressWarnings((sqrt(B * B - A * d2 * z3 * z3)- B)/A)
				}
			
				if (is.nan(z2) || is.infinite(z2) || !is.numeric(z2)) {
					# if we had a numerical problem then bisect
					z2 = z3 / 2
				}
				
				# don't accept too close to limits
				z2 = max(min(z2, INT * z3),(1 - INT) * z3)
				
				# update the step
				z1 = z1 + z2
				X = X + z2 %*% s
		  
				eval_ = f(X, P1, P2, P3 , P4, P5, P6)
				f2 = eval_$J
				df2 = eval_$grad
				
				M = M - 1
				
				# count epochs?!
				i = i + (runLength < 0)
				
				d2 = t(df2) %*% s
				
				# z3 is now relative to the location of z2
				z3 = z3 - z2
			}
			
			if (f2 > f1 + z1 * RHO * d1 || d2 > -SIG * d1) {
				# this is a failure
				break
			} else if (d2 > SIG * d1) {
				# success
				success = 1
				break
			} else if (M == 0) {
				# failure
				break
			}	
		
			# make cubic extrapolation
			A = 6 *(f2 -f3) / z3 +3 * (d2 + d3)
			B = 3 *(f3 -f2) - z3 * (d3 + 2*d2)
			
			# num. error possible - ok!
			z2 = suppressWarnings(-d2 * z3 * z3 / (B + sqrt(B * B - A * d2 * z3 * z3)))
			
			# num prob or wrong sign?
			if (!is.numeric(z2) || is.nan(z2) || is.infinite(z2) || z2 < 0) {
				# if we have no upper limit
				if (limit < -0.5) {
					# then extrapolate the maximum amount
					z2 = z1 * (EXT - 1)
				} else {
					# otherwise bisect
					z2 = (limit - z1) / 2
				}
			} else if ((limit > -0.5) && (z2+z1 > limit)) {
				# extrapolation beyond max?
				z2 = (limit - z1) / 2
			} else if ((limit < -0.5) && (z2 + z1 > z1 * EXT)) {
				# extrapolation beyond limit?
				
				# set to extrapolation limit
				z2 = z1 * (EXT - 1.0)
			} else  if (z2 < -z3 * INT) {
				z2 = -z3 * INT
			} else if ((limit > -0.5) && (z2 < (limit - z1) * (1.0 - INT))) {
				# too close to limit?
				z2 = (limit - z1) * (1.0 - INT)
			}
		
			# set point 3 equal to point 2
			f3 = f2
			d3 = d2
			z3 = -z2
			
			# update current estimates
			z1 = z1 + z2
			X = X + z2 %*% s
			
			eval_ = f(X, P1, P2, P3 , P4, P5, P6)
			f2 = eval_$J
			df2 = eval_$grad
			
			M = M - 1
			
			# count epochs?!
			i = i + (runLength < 0)
			
			d2 = t(df2) %*% s

		# end of line search
		} 
		
		# if line search succeeded
		if (success) {
			
			f1 = f2
			fX = f1
			
			print(paste0(S, i, "| Cost: ", f1))
			
			# Polack-Ribiere direction
			s = ((t(df2) %*% df2 - t(df1) %*% df2) / (t(df1) %*% df1)) * s - df2
			
			# swap derivatives
			tmp = df1
			df1 = df2
			df2 = tmp
			
			d2 = t(df1) %*% s

			# new slope must be negative
			if (d2 > 0)  {
				# otherwise use steepest direction
				s = -df1
				d2 = -t(s) %*% s
			}

			# slope ratio but max RATIO
			z1 = z1 * min(RATIO, d1 / (d2 - realmin))
			d1 = d2
			
			# this line search did not fail
			ls_failed = 0

		} else {
		
			# restore point from before failed line search
			X = X0
			f1 = f0
			df1 = df0

			# line search failed twice in a row
			if (ls_failed || i > abs(runLength)) {
				# or we ran out of time, so we give up
				break
			}

			# swap derivatives
			tmp = df1
			df1 = df2
			df2 = tmp
			
			# try steepest
			s = -df1
			d1 = -t(s) %*% s
			z1 = 1 / (1 - d1)
			
			# this line search failed
			ls_failed = 1
		}
	}
	
	return(list('X' = X, 'cost' = fX))
}
