smo2 <- function(X, Y, krnel, kpar1, kpar2, C, tol, steps, eps, method) {
##########################################################################
# FUNCTION
#  [alpha, b, w, evals, stp, glob] = SMO2(X, Y, kernel, kpar1, kpar2, C, tol, steps, eps, method)
# SMO2 algorithm of Platt with Keerthi modifications:
# ALL KERNEL EVALUATIONS ARE DONE IN THE BEGINING AND KEPT IN MEMORY
# Implements:
# 1. Sequential Minimal Optimization: A Fast Algorithm for Training Support Vector Machines
#    John C. Platt
# 2. Improvements to Platt's SMO Algorithm for SVM Classifier Design
#    S.S. Keerthi, S.K. Shevade, C. Bhattacharyya and K.R.K. Murthy
#    Technical Report CD-99-14
# The classifier that this second algorithm outputs is f(x)=wx-b
# Furthermore, note that in SMO Platt, due to a randomisation,
# the results might not be identical even for the same data
#
# INPUT ARGUMENTS:
# X:            Training points for both class - We assume A a column of n row-vectors
# Y:            Class values corresponding to training points - We assume a column of n values
# kernel:       Type of Kernel mapping to be used
#               'linear' : Linear (Default)
#               'poly' : Polynomial
#               'rbf' : Gauss
#               'sigmoid' : tanh
#
# kpar1:        1st parameter for kernel function (optional, default=0)
# kpar2:        1st parameter for kernel function (optional, default=0)
# C:            parameter which trades off wide margin with a small number of margin failures
# tol:          tollerance (Keerthi used 0.001
# steps:        maximum allowed steps to be taken (1st stopping condition)
# eps:          accuracy
# method:       0->Platt, 1->Keerthi modification 1, 2->Keerthi modification 2
#
# OUTPUT ARGUMENTS:
# alpha:        the column-vector of m+n Lagrange multipliers for each point. The
#               first m are for the m points of A set, and the next n for the n points of
#               B set.
# b:            the threshold value
# w:            the normal to the optimal separating hyperplane (meaningfull ONLY for
#               Linear kernel. (For test purposes only.)
# evals:        num of norm evaluations
# stp:          steps taken till the end
# flag:         a logical value which indicates whether or not the selected SVM
#               training algorihm has been terminated abnormally (1) or normally
#               (0). Abnormal termination means that no solution exists with the
#               specific kernel function choice and another kernel function should
#               be selected. The results obtained in this case are unreliable.
#
# Functions in this file have been provided by Michael Mavroforakis (c) 2003.
#
# Converted from Matlab into R by S.D. Separa (03/29/2016)
##########################################################################

	if (nargs() < 10) {
		method = 1
	}

	if (nargs() < 9) {
		eps = 0.0001
	}

	if (nargs() < 8) {
		steps = 10000 
	}

	if (nargs() < 7) {
		tol = 0.001
	}

	if (nargs() < 6) {
		C = Inf
	}

	if (nargs() < 5) {
		kpar2 = 0
	}

	if (nargs() < 4) {
		kpar1 = 0
	}

	if (nargs() < 3) {
		krnel = 0
	}

	if (nargs() < 2) {
		stop('Error: At least two arguments (training points and class values) must be supplied')
	}

	n = nrow(X)
	D = ncol(X)
	n1 = nrow(Y)
	D1 = ncol(Y)

	if (D1 != 1) {
		stop('Error: Class values cannot be vectors but real numbers')
	}

	if (n != n1) {
		stop('Error: Number of rows of X and Y must be the same (one class value for each sample)')
	}

	#here we should check if representatives of both classes are presented as training points

	# returns model(alpha, b, w, evals, stp, glob)
	model = list()
	if (method == 1) {
		model = SMO_Keerthi_modif1(X, Y, krnel, kpar1, kpar2, C, tol, steps, eps)
	} else if (method == 2) {
		model = SMO_Keerthi_modif2(X, Y, krnel, kpar1, kpar2, C, tol, steps, eps)
	} else {
		model = SMO_Platt(X, Y, krnel, kpar1, kpar2, C, tol, steps, eps)
	}

	# unpack model into objects
	for (i in 1:length(model)) assign(names(model)[i], model[[i]])
				
	if (method == 1 || method == 2) {
		flag = (glob$b_up < glob$b_low - 2 * tol) | (stp >= steps)
	} else {
		flag = (stp >= steps)
	}

	if (flag) {
		cat('The algorithm has not converged. This may be due to:\n (a) the maximum number of iterations has been reached and convergence has not, yet, been achieved or \n (b) the chosen values for the hyperparameters (C as well as the parameters that define the kernel function) \n can not lead to a solution. \n')
	}
	
	# unpack result into objects
	for (i in 1:length(model)) assign(names(model)[i], model[[i]])
	
	return(list('alpha' = alpha, 'w' = w, 'b' = b, 'stp' = stp, 'evals' = evals, 'glob' = glob))
}

SMO_Platt <- function(X, Y, krnel, kpar1, kpar2, C, tol, steps, eps) {
# returns [alpha, b, w, evals, stp, glob]

	n = nrow(X)
	D = ncol(X)
	
	# initialize alpha array to all zero
	alpha = array(0, c(n, 1))
	w = array(0, c(1,D))
	b = 0
	evals = 0

	K = array(0, c(n, n))
	
	for (i in 1:n) {
		K[, i] = CalcKernel(X, X[i, ], krnel, kpar1, kpar2)
	}

	# initialize struct for temporary variables that must be global
	glob = list('ecache' = numeric(0), 'v_1' = numeric(0), 'v_2' = numeric(0), 'I_0' = numeric(0) ,'ecache_f' = numeric(0))
	
	## initialize fcache array to all zero and its size to n
	glob$ecache = array(0, c(n,1))
	glob$ecache_f = array(0, c(n,1)) # 0 -> ecache value not-OK, 1 -> value OK
	glob$v_1 = which(Y == -1)
	glob$v_2 = which(Y ==1)

	stp = 0
	numChanged = 0
	examineAll = 1
	
	while ((numChanged > 0 || examineAll == 1) && stp <= steps) {
		
		numChanged = 0
		
		if (examineAll == 1) {
			for (i in 1:n) {
				stp = stp + 1
				
				if (stp > steps) {
					break
				}

				result = examineExampleP(i, glob, alpha, w, b, X, Y, krnel, kpar1, kpar2, C, tol, steps, stp, evals, eps, K)
				# unpack result into objects
				for (i in 1:length(result)) assign(names(result)[i], result[[i]])
				
				numChanged = numChanged + retval
			}
		} else {
			
			glob$I_0 = which(alpha > 0 && alpha < C)
			k = length(glob$I_0)
			
			for (i in 1:k) {
			
				stp = stp + 1
				
				if (stp > steps) {
					break
				}
				
				# glob$I_0 changes inside loop (in examineExampleP)
				if (i > length(glob$I_0)) {
					break
				} 
				
				# returns [retval, alpha, w, b, stp, evals, glob]
				result = examineExampleP(glob$I_0[i], glob, alpha, w, b, X, Y, krnel, kpar1, kpar2, C, tol, steps, stp, evals, eps, K)
				
				# unpack result into objects
				for (i in 1:length(result)) assign(names(result)[i], result[[i]])
				
				numChanged = numChanged + retval
			}
		}
		
		if (examineAll == 1) {
			examineAll = 0
		} else if (numChanged == 0) {
			examineAll = 1
		}
	}
	
	return(list('alpha' = alpha, 'w' = w, 'b' = b, 'stp' = stp, 'evals' = evals, 'glob' = glob))
}

examineExampleP <- function(i2, glob, alpha, w, b, X, Y, krnel, kpar1, kpar2, C, tol, steps, stp, evals, eps,K) {
# returns [retval, alpha, w, b, stp, evals, glob]

	retval = 0
	n = nrow(X)
	D = ncol(X)
	y2 = Y[i2]
	
	alph2 = alpha[i2]
	
	if (glob$ecache_f[i2] == 1) {
		E2 = glob$ecache[i2]
	} else {
		ki2 = as.vector(K[,i2])
		evals = evals + n
		E2 = -y2 + (t(ki2) %*% (Y * alpha)) - b
		glob$ecache[i2] = E2
		glob$ecache_f[i2] = 1
	}
	
	r2 = E2 * y2
	
	if ((r2 < -tol && alph2 < C) || (r2 > tol && alph2 > 0)) {
		if (length(glob$I_0) > 1) {
			#i1 = result of second choice heuristic
			v = which(glob$ecache_f==1)
			k = length(v)
			Emax = 0
			for (i in 1:k) {
				tmp = abs(glob$ecache[v[i]] - E2)
				if (tmp > Emax) {
					Emax = tmp
					i1 = v[i]
				}
			}
			
			stp = stp + 1
			
			result = takeStepP(i1, i2, glob, alpha, w, b, X, Y, krnel, kpar1, kpar2, C, tol, evals, eps, K)
			
			# unpack result into objects
			for (i in 1:length(result)) assign(names(result)[i], result[[i]])
			
			if (retval == 1) {
				return(list('retval' = retval, 'alpha' = alpha, 'w' = w, 'b' = b, 'stp' = stp, 'evals' = evals, 'glob' = glob))
			}
		}
		
		# loop over all non-zero and non-C alpha, starting at a random point
		k = length(glob$I_0)
		set.seed(2)
		r = as.integer(floor(k * rand()))
		
		for (i in 1:k) {
			i1 = mod(r + i, k) + 1
			stp =stp + 1
			
			result = takeStepP(i1, i2, glob, alpha, w, b, X, Y, krnel, kpar1, kpar2, C, tol, evals, eps, K)
			
			# unpack result into objects
			for (i in 1:length(result)) assign(names(result)[i], result[[i]])
			
			if (retval == 1) {
				return(list('retval' = retval, 'alpha' = alpha, 'w' = w, 'b' = b, 'stp' = stp, 'evals' = evals, 'glob' = glob))
			}
		}
		
		# loop over all possible i1, starting at a random point
		k = n
		r = as.integer(floor(k * rand()))
		
		for (i in 1:k) {
		
			i1 = mod(r + i, k) + 1
			stp = stp + 1
			
			result = takeStepP(i1, i2, glob, alpha, w, b, X, Y, krnel, kpar1, kpar2, C, tol, evals, eps,K)
			
			# unpack result into objects
			for (i in 1:length(result)) assign(names(result)[i], result[[i]])
			
			if (retval == 1) {
				return(list('retval' = retval, 'alpha' = alpha, 'w' = w, 'b' = b, 'stp' = stp, 'evals' = evals, 'glob' = glob))
			}
		}
	}
	
	return(list('retval' = retval, 'alpha' = alpha, 'w' = w, 'b' = b, 'stp' = stp, 'evals' = evals, 'glob' = glob))
}

takeStepP <- function(i1, i2, glob, alpha, w, b, X, Y, krnel, kpar1, kpar2, C, tol, evals, eps, K) {
# returns [retval, alpha, w, b, evals, glob]

	# for strcmpi
	require(pracma)
	
	n = nrow(X)
	D = ncol(X)

	if (i1 == i2) {
		retval = 0
		return(list('retval' = retval, 'alpha' = alpha, 'w' = w, 'b' = b, 'evals' = evals, 'glob' = glob))
	}

	alph1 = alpha[i1]
	y1 = Y[i1]
	alph2 = alpha[i2]
	y2 = Y[i2]
	s = y1 * y2
	
	# Compute L, H
	if (y1 != y2) {
		L = max(0, alph2 - alph1)
		H = min(C, C + alph2 - alph1)
	} else { # y1 == y2
		L = max(0, alph1 + alph2 - C)
		H = min(C, alph1 + alph2)
	}
	
	if (L == H) {
		retval = 0
		return(list('retval' = retval, 'alpha' = alpha, 'w' = w, 'b' = b, 'evals' = evals, 'glob' = glob))
	}
	
	# calculate E1 = SVM output in X[i1] - y1 (check in error cache)
	if (glob$ecache_f[i1] == 0) {
		ki1 = as.vector(K[, i1])
		evals = evals + n
		E1 = -y1 + (t(ki1) %*% (Y * alpha)) - b
		glob$ecache[i1] = E1
		glob$ecache_f[i1] = 1
	} else {
		E1 = glob$ecache[i1]
	}
	
	# calculate E2 = SVM output in X[i2] - y2 (check in error cache)
	if (glob$ecache_f[i2] == 0) {
		ki2 = as.vector(K[, i2])
		evals = evals + n
		E2 = -y1 + (t(ki2) %*% (Y * alpha)) - b
		glob$ecache[i2] = E2
		glob$ecache_f[i2] = 1
	} else {
		E2 = glob$ecache[i2]
	}
	## computation of the derivative eta
	k11 = K[i1, i1]
	k12 = K[i2, i1]
	k22 = K[i2, i2]
	evals = evals + 3
	eta = -(2 * k12) + k11 + k22
	
	## computation of new alpha(i2)
	if (eta > 0) {
		a2 = alph2 + (y2 * (E1 - E2) / eta)
		if (a2 < L) {
			a2 = L
		} else if (a2 > H) {
			a2 = H
		}
	} else { ## the derivative is 0 => we have to make optimization by other means
		## Lobj = objective function at a2=L (according to Platt)
		## Hobj = objective function at a2=H (according to Platt)
		L1 = alph1 + (s * (alph2 - L))
		H1 = alph1 + (s * (alph2 - H))
		f1 = y1 * (E1 + b) - (alph1 * k11) - (s * alph2 * k12)
		f2 = y2 * (E2 + b) - (alph2 * k22) - (s * alph1 * k12)
		
		Lobj = (L1 * f1) + (L * f2) + (0.5 * k11 * L1^2) + (0.5 * k22 * L^2) + (s * k12 * L * L1)
		Hobj = (H1 * f1) + (H * f2) + (0.5 * k11 * H1^2) + (0.5 * k22 * H^2) + (s * k12 * H * H1)
		
		if (Lobj < Hobj - eps) {
			a2 = L
		} else if (Lobj > Hobj + eps) {
			a2 = H
		} else {
			a2 = alph2
		}
	}
	
	if (abs(a2 - alph2) < (eps * (a2 + alph2 + eps))) {
		retval = 0
		return(list('retval' = retval, 'alpha' = alpha, 'w' = w, 'b' = b, 'evals' = evals, 'glob' = glob))
	}
	
	# computation on new a1pha1(a1)
	a1 = alph1 + (s * (alph2 - a2))
	# Update threshold to reflect change in Lagrange multipliers
	b_old = b
	
	if (a1 > L && a1 < H) {
		b = E1 + (y1 * (a1 - alph1) * k11) + (y2 * (a2 - alph2) * k12) + b
	} else if (a2 > L && a2 < H) {
		b = E2 + (y1 * (a1 - alph1) * k12) + (y2 * (a2 - alph2) * k22) + b
	} else {
		b1 = E1 + (y1 * (a1 - alph1) * k11) + (y2 * (a2 - alph2) * k12) + b
		b2 = E2 + (y1 * (a1 - alph1) * k12) + (y2 * (a2 - alph2) * k22) + b
		b = (b1 + b2) / 2
	}
	
	# Update weight vector to reflect change in a1 & a2, if linear SVM
	if (strcmpi(krnel, 'linear') == 1) {
		w = w + (y1 * (a1 - alph1) * X[i1, ]) + (y2 * (a2 - alph2) * X[i2, ])
	}
	
	# Update ecache[i] using new Lagrange multipliers
	v = which(glob$ecache_f == 1)
	
	for (i in 1:length(v)) {
		ki1 = K[v[i], i1]
		ki2 = K[v[i], i2]
		evals = evals + 2
		glob$ecache[v[i]] = glob$ecache[v[i]] + b_old - b + (y1 * (a1 - alph1) * ki1) + (y2 * (a2 - alph2) * ki2)
	}
	
	# Store a1 and a2 in the alpha array
	alpha[i1] = a1
	alpha[i2] = a2

	## Compute updated E values for i1 and i2
	# glob$ecache(i1) = E1 + b - b_old - (y1 * (a1 - alph1) * k11) - (y2 * (a2 - alph2) * k12)
	glob$ecache[i1] = E1 + b_old - b + (y1 * (a1 - alph1) * k11) + (y2 * (a2 - alph2) * k12)
	glob$ecache_f[i1] = 1
	# glob$ecache(i2) = E2 + b - b_old - (y1 * (a1 - alph1) * k12) - (y2 * (a2 - alph2) * k22)
	glob$ecache[i2] = E2 + b_old - b + (y1 * (a1 - alph1) * k12) + (y2 * (a2 - alph2) * k22)
	glob$ecache_f[i2] = 1
	# Update I_0
	glob$I_0 = which(alpha > 0 && alpha < C)

	retval = 1
	
	return(list('retval' = retval, 'alpha' = alpha, 'w' = w, 'b' = b, 'evals' = evals, 'glob' = glob))
}

SMO_Keerthi_modif1 <- function(X, Y, krnel, kpar1, kpar2, C, tol, steps, eps) {
# returns [alpha, b, w, evals, stp, glob]

	n = nrow(X)
	D = ncol(X)
	
	# initialize alpha array to all zero
	# --------------------------------------------------------------------------
	alpha = array(0, c(n, 1))
	w = array(0, c(1, D))
	b = 0
	evals = 0
	# --------------------------------------------------------------------------

	K = array(0, c(n, n))
	
	for (i in 1:n) {
		K[, i] = CalcKernel(X, X[i, ], krnel, kpar1, kpar2)
	}

	# initialize struct for temporary variables that must be global
	# --------------------------------------------------------------------------

	glob = list('fcache' = numeric(0), 'b_up' = numeric(0), 'b_low' = numeric(0), 'i_up' = numeric(0), 'i_low' = numeric(0), 'v_1' = numeric(0), 'v_2' = numeric(0), 'I_0' = numeric(0), 'I_1' = numeric(0), 'I_2' = numeric(0), 'I_3' = numeric(0), 'I_4' = numeric(0))
	
	## initialize fcache array to all zero and its size to n
	glob$fcache = array(0, c(n, 1))
	# initialize b_up = -1, i_up to any one index of class 1
	glob$b_up = -1
	glob$v_1 = which(Y == 1)
	glob$i_up = glob$v_1[1]
	# initialize b_low = 1, i_low to any one index of class 2
	glob$b_low = 1
	glob$v_2 = which(Y == -1)
	glob$i_low = glob$v_2[1]
	# set fcache[i_low] = 1 and fcache[i_up] = -1
	glob$fcache[glob$i_low] = 1
	glob$fcache[glob$i_up] = -1

	# Initialize the I_* sets
	glob$I_0 = which(alpha > 0 && alpha < C)
	glob$I_1 = which(alpha[glob$v_1] == 0)
	glob$I_1 = glob$v_1[glob$I_1]
	glob$I_2 = which(alpha[glob$v_2] == C)
	glob$I_2 = glob$v_2[glob$I_2]
	glob$I_3 = which(alpha[glob$v_1] == C)
	glob$I_3 = glob$v_1[glob$I_3]
	glob$I_4 = which(alpha[glob$v_2] == 0)
	glob$I_4 = glob$v_2[glob$I_4]
	# --------------------------------------------------------------------------

	stp = 0
	numChanged = 0
	examineAll = 1
	
	while (((numChanged > 0) || (examineAll == 1)) && (stp <= steps)) {
		numChanged = 0
		
		if (examineAll==1) {
			for (i in 1:n) {
				
				stp = stp + 1
				
				if (stp > steps) {
					break
				}
				
				result = examineExampleK(i, glob, alpha, w, b, X, Y, krnel, kpar1, kpar2, C, tol, steps, stp, evals, eps, K)
				
				# unpack result into objects
				for (i in 1:length(result)) assign(names(result)[i], result[[i]])
				
				numChanged = numChanged + retval
			}
		} else {
			k = length(glob$I_0)
			
			for (i in 1:k) {
				
				stp = stp + 1
				
				if (stp > steps) {
					break
				}
				
				# glob$I_0 changes inside loop (in examineExampleK)
				if (i > length(glob$I_0)) {
					break 
				}
				 
				result = examineExampleK(glob$I_0[i], glob, alpha, w, b, X, Y, krnel, kpar1, kpar2, C, tol, steps, stp, evals, eps,K)
				
				# unpack result into objects
				for (i in 1:length(result)) assign(names(result)[i], result[[i]])
				
				numChanged = numChanged + retval
				
				# it is easy to check if optimality on I_0 is attained...
				if ( (glob$b_up) > ( glob$b_low - (2*tol) ) ) {
					# exit the loop after setting numChanged = 0
					numChanged = 0
					break # [sdsepara] long comment (in Greek) removed here
				}
			}
		}
		
		if (examineAll == 1) {
			examineAll = 0
		} else if (numChanged == 0) {
			examineAll = 1
		}
	}
	
	b = (glob$b_up + glob$b_low) / 2
	
	return(list('alpha' = alpha, 'w' = w, 'b' = b, 'stp' = stp, 'evals' = evals, 'glob' = glob))
}

SMO_Keerthi_modif2 <- function(X, Y, krnel, kpar1, kpar2, C, tol, steps, eps) {
# returns [alpha, b, w, evals, stp, glob]

	n = nrow(X)
	D = ncol(X)

	# initialize alpha array to all zero
	# --------------------------------------------------------------------------
	alpha = array(0, c(n, 1))
	w = array(0, c(1, D))
	b = 0
	evals = 0
	# --------------------------------------------------------------------------

	K = array(0, c(n, n))
	
	for (i in 1:n) {
		K[, i]= CalcKernel(X, X[i, ], krnel, kpar1, kpar2)
	}
	
	# initialize struct for temporary variables that must be global
	# --------------------------------------------------------------------------

	glob = list('fcache' = numeric(0), 'b_up' = numeric(0), 'b_low' = numeric(0), 'i_up' = numeric(0), 'i_low' = numeric(0), 'v_1' = numeric(0), 'v_2' = numeric(0), 'I_0' = numeric(0), 'I_1' = numeric(0), 'I_2' = numeric(0), 'I_3' = numeric(0), 'I_4' = numeric(0))
	## initialize fcache array to all zero and its size to n
	glob$fcache = array(0, c(n, 1))
	# initialize b_up = -1, i_up to any one index of class 1
	glob$b_up = -1
	glob$v_1 = which(Y == 1)
	glob$i_up = glob$v_1[1]
	# initialize b_low = 1, i_low to any one index of class 2
	glob$b_low = 1
	glob$v_2 = which(Y == -1)
	glob$i_low = glob$v_2[1]
	# set fcache[i_low] = 1 and fcache[i_up] = -1
	glob$fcache[glob$i_low] = 1
	glob$fcache[glob$i_up] = -1

	# Initialize the I_* sets
	glob$I_0 = which(alpha > 0 && alpha < C)
	glob$I_1 = which(alpha[glob$v_1] == 0)
	glob$I_1 = glob$v_1[glob$I_1]
	glob$I_2 = which(alpha[glob$v_2] == C)
	glob$I_2 = glob$v_2[glob$I_2]
	glob$I_3 = which(alpha[glob$v_1] == C)
	glob$I_3 = glob$v_1[glob$I_3]
	glob$I_4 = which(alpha[glob$v_2] == 0)
	glob$I_4 = glob$v_2[glob$I_4]
	# --------------------------------------------------------------------------


	stp = 0
	numChanged = 0
	examineAll = 1
	
	while ((numChanged > 0 || examineAll == 1) && stp <= steps) {
		numChanged = 0
		
		if (examineAll == 1) {
			
			for (i in 1:n) {
				
				stp = stp + 1
				
				if (stp > steps) {
					break
				}
				
				result = examineExampleK(i, glob, alpha, w, b, X, Y, krnel, kpar1, kpar2, C, tol, steps, stp, evals, eps,K)
				
				# unpack result into objects
				for (i in 1:length(result)) assign(names(result)[i], result[[i]])
				
				numChanged = numChanged + retval
			}
			
		} else {
			
			#the following loop is the only difference between the two SMO
			#modifications. Whereas, in modification 1, the inner loop selects
			#i2 from I_0 sequentially, here i2 is always set to the current
			#i_low and i1 is set to the current i_up clearly, this corresponds
			#to choosing the worst violating pair using members of I_0 and some
			#other indices.
			
			inner_loop_success = 1
			
			while (((glob$b_up) < (glob$b_low - (2 * tol))) && (inner_loop_success != 0)) {
				i2 = glob$i_low
				y2 = Y[i2]
				alph2 = alpha[i2]
				F2 = glob$fcache[i2]
				i1 = glob$i_up
				stp = stp + 1
				
				if (stp > steps) {
					break
				}
				
				stp = stp + 1
				
				result = takeStepK(glob$i_up, glob$i_low, glob, alpha, w, b, X, Y, krnel, kpar1, kpar2, C, tol, evals, eps,K)
				
				# unpack result into objects
				for (i in 1:length(result)) assign(names(result)[i], result[[i]])
				
				numChanged = numChanged + retval
			}
			
			num_changed = 0
		}
		
		if (examineAll == 1) {
			
			examineAll = 0
			
		} else if (numChanged == 0) {
			
			examineAll = 1
		}
	}
	
	b = (glob$b_up + glob$b_low) / 2
	
	return(list('alpha' = alpha, 'w' = w, 'b' = b, 'stp' = stp, 'evals' = evals, 'glob' = glob))
}

takeStepK <- function(i1, i2, glob, alpha, w, b, X, Y, krnel, kpar1, kpar2, C, tol, evals, eps, K) {
# returns [retval, alpha, w, b, evals, glob]    
# Much of this procedure is same as that in Platt's SMO pseudo-code.

	if (i1 == i2) {
		retval = 0
		return(list('retval' = retval, 'alpha' = alpha, 'w' = w, 'b' = b, 'evals' = evals, 'glob' = glob))
	} 
	
	alph1 = alpha[i1]
	y1 = Y[i1]
	F1 = glob$fcache[i1]
	alph2 = alpha[i2]
	y2 = Y[i2]
	F2 = glob$fcache[i2]
	s = y1 * y2

	# Compute L, H - If L = H return 0
	# --------------------------------------------------------------------------
	if (y1 != y2) {
		L = max(0, alph2 - alph1)
		H = min(C, C + alph2 - alph1)
	} else { # y1 == y2
		L = max(0, alph1 + alph2 - C)
		H = min(C, alph1 + alph2)
	}
	
	if (L == H) {
		retval = 0
		return(list('retval' = retval, 'alpha' = alpha, 'w' = w, 'b' = b, 'evals' = evals, 'glob' = glob))
	}
	
	# --------------------------------------------------------------------------


	## computation of the derivative eta
	# --------------------------------------------------------------------------
	k11 = K[i1, i1]
	k12 = K[i1, i2]
	k22 = K[i2, i2]
	evals = evals + 3
	eta = (2*k12) - k11 - k22
	
	## computation of new alpha(i2)
	if (eta<0) {
		a2 = alph2 - (y2 * (F1 - F2) / eta) ## HERE it is different from Platt
		if (a2 < L) {
			a2 = L
		} else if (a2 > H) {
			a2 = H
		}
	} else {
		##the derivative is 0 => we have to make optimization by other means
		##Lobj = objective function at a2=L (according to Platt)
		##Hobj = objective function at a2=H (according to Platt)
		L1 = alph1 + (s *(alph2 - L))
		H1 = alph1 + (s *(alph2 - H))
		f1 = (-y1 * F1) + (alph1 * k11) + (s * alph2 * k12)
		f2 = (-y2 * F2) + (alph2 * k22) + (s * alph1 * k12)
		Lobj = (L1 * f1) + (L * f2) - (0.5 * k11 * L1^2) - (0.5 * k22 * L^2) - (s * k12 * L * L1)
		Hobj = (H1 * f1) + (H * f2) - (0.5 * k11 * H1^2) - (0.5 * k22 * H^2) - (s * k12 * H * H1)
		
		if (Lobj > Hobj + eps) {
			a2 = L
		} else if (Lobj < Hobj - eps) {
			a2 = H
		} else {
			a2 = alph2
		}
	}
	# --------------------------------------------------------------------------


	# Calculate the change in a - if very small  return 0
	## !!!!!NOTE!! if eps not small enought it may stop the algorithm early!!!!!!
	# --------------------------------------------------------------------------
	if (abs(a2 - alph2)< eps * (a2 + alph2 + eps)) {
		retval = 0
		return(list('retval' = retval, 'alpha' = alpha, 'w' = w, 'b' = b, 'evals' = evals, 'glob' = glob))
	}
	#--------------------------------------------------------------------------


	# computation on new a1pha1(a1)
	#--------------------------------------------------------------------------
	a1 = alph1 + (s*(alph2-a2))
	#--------------------------------------------------------------------------


	#Update weight vector to reflect change in a1 & a2, if linear SVM
	#--------------------------------------------------------------------------
	if (strcmpi(krnel, 'linear')) {
		w = w + (y1 * (a1 - alph1) * X[i1, ]) + (y2 * (a2 - alph2) * X[i2, ])
	}
	#--------------------------------------------------------------------------


	# Store a1 and a2 in the alpha array
	# --------------------------------------------------------------------------
	alpha[i1] = a1
	alpha[i2] = a2
	#--------------------------------------------------------------------------

	#Update fcache[i] for i in I_0 using new Lagrange multipliers
	#--------------------------------------------------------------------------
	k = length(glob$I_0)
	
	for (i in 1:k) {
		ki1 = K[glob$I_0[i], i1]
		ki2 = K[glob$I_0[i], i2]
		evals = evals + 2
		glob$fcache[glob$I_0[i]] = glob$fcache[glob$I_0[i]] + (y1 * (a1 - alph1) * ki1) + (y2 * (a2 - alph2) * ki2)
	}
	# --------------------------------------------------------------------------


	# The update below is simply achieved by keeping and updating information
	# about alpha_i being at 0, C or in between them. Using this together with
	# target[i] gives information as to which index set i belongs.
	# Update I_0, I_1, I_2, I_3, I_4
	#--------------------------------------------------------------------------
	glob$I_0 = which(alpha > 0 && alpha < C)
	glob$I_1 = which(alpha[glob$v_1] == 0)
	glob$I_1 = glob$v_1[glob$I_1]
	glob$I_2 = which(alpha[glob$v_2] == C)
	glob$I_2 = glob$v_2[glob$I_2]
	glob$I_3 = which(alpha[glob$v_1] == C)
	glob$I_3 = glob$v_1[glob$I_3]
	glob$I_4 = which(alpha[glob$v_2] == 0)
	glob$I_4 = glob$v_2[glob$I_4]
	#--------------------------------------------------------------------------

	# Compute updated F values for i1 and i2
	#--------------------------------------------------------------------------
	glob$fcache[i1] = F1 + (y1 * (a1 - alph1) * k11) + (y2 * (a2 - alph2) * k12)
	glob$fcache[i2] = F2 + (y1 * (a1 - alph1) * k12) + (y2 * (a2 - alph2) * k22)
	# --------------------------------------------------------------------------


	#Compute (i_low, b_low) and (i_up, b_up),
	#using only i1, i2 and indices in I_0

	#--------------------------------------------------------------------------
	#--GIA TO i1 -------------------------------------------------------
	v = which(glob$I_1 == i1)
	i1_in_I_1 = length(v)
	v = which(glob$I_2 == i1)
	i1_in_I_2 = length(v)
	v = which(glob$I_3 == i1)
	i1_in_I_3 = length(v)
	v = which(glob$I_4 == i1)
	i1_in_I_4 = length(v)

	## --------------------------------------------------------------------------
	# --GIA TO i2 -------------------------------------------------------
	v = which(glob$I_1 == i2)
	i2_in_I_1 = length(v)
	v = which(glob$I_2 == i2)
	i2_in_I_2 = length(v)
	v = which(glob$I_3 == i2)
	i2_in_I_3 = length(v)
	v = which(glob$I_4 == i2)
	i2_in_I_4 = length(v)


	# 1)First Compute i_low, i_up for I_0
	# --------------------------------------------------------------------------

	if (length(glob$I_0) != 0) {   # Trying to run the smo mod1 for the alult datasets I diskovered that for small values of C there was this problem.
		glob$b_up = min(glob$fcache[glob$I_0]) 
		glob$i_up = which.min(glob$fcache[glob$I_0])
		glob$i_up = glob$I_0[glob$i_up]
		
		if (length(glob$i_up) != 1) {
			glob$i_up = glob$i_up[1]
		}
		
		glob$b_low = max(glob$fcache[glob$I_0])
		glob$i_low = which.max(glob$fcache[glob$I_0])
		glob$i_low = glob$I_0[glob$i_low]
		
		if (length(glob$i_low) != 1) {
			glob$i_low = glob$i_low[1]
		}
	}

	# 2)Then check if i1 or i2 should replace i_up or i_low

	# 2a) For i1

	if ((glob$b_up > glob$fcache[i1]) && (i1_in_I_1 + i1_in_I_2)) {
		glob$b_up = glob$fcache[i1]
		glob$i_up = i1
	}
	
	if ((glob$b_low < glob$fcache[i1])&&(i1_in_I_3 + i1_in_I_4)) {
		glob$b_low = glob$fcache[i1]
		glob$i_low = i1
	}

	# 2b) For i2

	if ((glob$b_up > glob$fcache[i2]) && (i2_in_I_1 + i2_in_I_2)) {
		glob$b_up = glob$fcache[i2]
		glob$i_up = i2
	}
	
	if ((glob$b_low < glob$fcache[i2]) && (i2_in_I_3 + i2_in_I_4)) {
		glob$b_low = glob$fcache[i2]
		glob$i_low = i2
	}
	
	retval = 1
	
	return(list('retval' = retval, 'alpha' = alpha, 'w' = w, 'b' = b, 'evals' = evals, 'glob' = glob))
}

examineExampleK <- function(i2, glob, alpha, w, b, X, Y, krnel, kpar1, kpar2, C, tol, steps, stp, evals, eps,K) {
# returns [retval, alpha, w, b, stp, evals, glob]

	retval = 0
	n = nrow(X)
	D = ncol(X)

	y2 = Y[i2]
	
	alph2 = alpha[i2]
	
	v = which(glob$I_0 == i2)
	i2_in_I_0 = length(v)
	
	if (i2_in_I_0 > 0) {
		F2 = glob$fcache[i2]
	} else {
		ki2 = as.vector(K[, i2])
		evals = evals + n
		F2 = -y2 + (t(ki2) %*% (Y * alpha))
		glob$fcache[i2] = F2
	}

	# Update (b_low, i_low) or (b_up, i_up) using (F2,#i2)
	# --------------------------------------------------------------------------
	v = which(glob$I_1 == i2)
	i2_in_I_1 = length(v)
	v = which(glob$I_2 == i2)
	i2_in_I_2 = length(v)
	v = which(glob$I_3 == i2)
	i2_in_I_3 = length(v)
	v = which(glob$I_4 == i2)
	i2_in_I_4 = length(v)

	if ((i2_in_I_1 + i2_in_I_2 > 0) && (F2 < glob$b_up)) {
		glob$b_up = F2
		glob$i_up = i2
	} else if ((i2_in_I_3 + i2_in_I_4 > 0) && (F2 > glob$b_low)) {
		glob$b_low = F2
		glob$i_low = i2
	}
	# ------------------------------------------------------------------


	#Chech optimality using current b_low and b_up and, if
	#violated, find an index i1 to do joint optimization with i2
	#------------------------------------------------------------------
	optimality = 1
	
	if ((i2_in_I_0 + i2_in_I_1 + i2_in_I_2) > 0) {
		if ((glob$b_low - F2) > (2 * tol)) {
			optimality = 0
			i1 = glob$i_low
		}
	}
	
	if ((i2_in_I_0 + i2_in_I_3 + i2_in_I_4) > 0) {
		if ((F2 - glob$b_up) > (2 * tol)) {
			optimality = 0
			i1 = glob$i_up
		}
	}
	
	if (optimality == 1) {
		retval = 0
		return(list('retval' = retval, 'alpha' = alpha, 'w' = w, 'b' = b, 'stp' = stp, 'evals' = evals, 'glob' = glob))
	}
	# ------------------------------------------------------------------



	# For i2 in I_0 choose the better i1
	# ------------------------------------------------------------------

	if (i2_in_I_0 > 0) {
		if ((glob$b_low - F2) > (F2 - glob$b_up)) {
			i1 = glob$i_low
		} else {
			i1 = glob$i_up
		}
	}
	# ------------------------------------------------------------------

	stp = stp + 1

	result =  takeStepK(i1, i2, glob, alpha, w, b, X, Y, krnel, kpar1, kpar2, C, tol, evals, eps, K)
	
	# unpack result into objects
	for (i in 1:length(result)) assign(names(result)[i], result[[i]])
	
	return(list('retval' = retval, 'alpha' = alpha, 'w' = w, 'b' = b, 'stp' = stp, 'evals' = evals, 'glob' = glob))
}

CalcKernel <- function(u, v, ker, kpar1, kpar2) {
##########################################################################
# FUNCTION
#  k = CalcKernel(u, v, ker, kpar1, kpar2)
# Calculates the Kernel function between two points (x1, x2).
#
# INPUTS ARGUMENTS:
#  ker:     Type of Kernel mapping to be used
#             'linear' : Linear 
#             'poly' : Polynomial
#             'rbf' : Gaussian  
#             'sigmoid' : tanh   
#  u:       row vector representing 1st point, or column of row-vectors
#           representing array of 1st points.
#  v:       row vector representing 2nd point.
#  kpar1:   1st parameter for kernel function (optional, default=0).
#  kpar2:   1st parameter for kernel function (optional, default=0).
#
# OUTPUT ARGUMENTS:
#  k:       the value of kernel function for these two points. If u is a
#           matrix (column of row-vectors), k is a column of values, with
#           same rows as u (one value for each row of u).
#
# (c) 2010 S. Theodoridis, A. Pikrakis, K. Koutroumbas, D. Cavouras
#
# Converted from Matlab into R by S.D. Separa (03/29/2016)
##########################################################################

	# for strcmp
	require(pracma)
	
	if (nargs() < 3) {
		stop('CalcKernel needs at least 3 arguments')
	}
	
	if (nargs() < 5) {
		kpar2 = 0
	}
	
	if (nargs() < 4) {
		kpar1 = 0
	}

	# make sure u and v have same number of columns
	if (is.vector(v)) {
		v = array(v, c(1, length(v)))
	}
	
	if (is.vector(u)) {
		u = array(u, c(1, length(u)))
	}
	
	r1 = nrow(u)
	c1 = ncol(u)
	r2 = nrow(v)
	c2 = ncol(v)
	
	if (r1 < 1 || r2 != 1) {
		stop('CalcKernel expect u=column of row-vectors and v a row-vector')
	}
	
	if (c1 != c2) {
		stop('CalcKernel needs both x1 and x2 to have same num of columns')
	}

	ker = tolower(ker)
	
	k = numeric(0)
	
	if (strcmp(ker, 'linear')) {
		k = u %*% t(v)
	} else if (strcmp(ker, 'poly')) {
		k = (u %*% t(v) + kpar1)^kpar2
	} else if (strcmp(ker, 'rbf')) {
			
		k = array(0, c(r1, 1))
		
		for (i in 1:r1) {
			k[i] = exp(-(u[i,] - v) %*% t(u[i, ]-v)/(2*kpar1^2))
		}
	} else if (strcmp(ker, 'erbf')) {
		
		k = array(0, c(r1, 1))
		
		for (i in 1:r1) {
			k[i] = exp(-sqrt((u[i, ] - v) %*% t(u[i, ]-v))/(2*kpar1^2))
		}
	} else if (strcmp(ker, 'sigmoid')) {

		k = array(0, c(r1, 1))
		
		for (i in 1:r1) {
			k[i] = tanh(kpar1*u[i, ] %*% t(v)/length(u[i, ]) + kpar2)
		}
	} else if (strcmp(ker, 'fourier')) {
		
		k = array(0, c(r1, 1))
		
		for (i in 1:r1) {
			z = sin(kpar1 + 1 / 2)* 2 * array(1, c(length(u[i, ]),1))
			j = which(u[i,] - v)
			z[j] = sin(kpar1 + 1/2) * (u[i, j] - v[j])/sin((u[i, j] - v[j]) / 2)
			k[i] = prod(z)
		}
	} else if (strcmp(ker, 'spline')) {
	
		k = array(0, c(r1, 1))
		
		for (i in 1:r1) {
			z = 1 + u[i, ] * v + u[i, ] * v * bsxfun('min', u[i, ], v) - ((u[i, ] + v) / 2) * (bsxfun('min', u[i, ], v)) ^ 2 + (1 / 3) * (bsxfun('min', u[i, ], v)) ^ 3
			k[i] = prod(z)
		}
	}  else if (strcmp(ker, 'curvspline') || strcmp(ker, 'anova')) {
			
		k = array(0, c(r1, 1))
		
		for (i in 1:r1) {
			z = 1 + u[i, ] * v + (1 / 2) * u[i, ] *v * bsxfun('min', u[i, ], v) - (1 / 6) * bsxfun('min', u[i, ], v) ^ 3
			k[i] = prod(z)
		}
		
	} else if (strcmp(ker, 'bspline')) {
		
		k = array(0, c(r1, 1))
		
		for (i in 1:r1) {
			z = 0
			for (r in 0:(2*(kpar1+1))) {
				z = z + (-1) ^ r * choose(2*(kpar1+1), r) * (bsxfun('max', array(0, length(v)), u[i, ] - v + kpar1 + 1 - r)) ^ (2 * kpar1 + 1)
			}
			
			k[i] = prod(z)
		}
	} else if (strcmp(ker, 'anovaspline1')) {
		
		k = array(0, c(r1, 1))
		
		for (i in 1:r1) {
			z = 1 + u[i, ] * v + u[i, ] * v * bsxfun('min', u[i, ], v) - ((u[i, ] + v) / 2) * bsxfun('min', u[i, ], v) ^ 2 + (1 / 3) * bsxfun('min', u[i, ], v) ^ 3
			k[i] = prod(z)
		}
	} else if (strcmp(ker, 'anovaspline2')) {
		
		k = array(0, c(r1, 1))
		
		for (i in 1:r1) {
			z = 1 + u[i, ] * v + (u[i, ] * v)^2 + (u[i, ] * v) ^ 2 * bsxfun('min', u[i, ], v) - u[i, ] * v * (u[i, ]+ v) * bsxfun('min', u[i, ], v) ^ 2 + (1 / 3) * (u[i, ] ^ 2 + 4 * u[i, ] * v + v ^ 2) * bsxfun('min', u[i, ], v) ^3 - (1 / 2) * (u[i, ] + v) * bsxfun('min', u[i, ], v) ^ 4 + (1 / 5) * bsxfun('min', u[i, ], v) ^ 5
			k[i] = prod(z)
		}
	} else if (strcmp(ker, 'anovaspline3')) {
		
		k = array(0, c(r1, 1))
		
		for (i in 1:r1) {
			z = 1 + u[i, ] * v + (u[i, ] * v) ^ 2 + (u[i, ] * v) ^ 3 + (u[i, ] * v) ^ 3 * bsxfun('min', u[i, ], v) - (3 / 2) * (u[i, ] * v) ^ 2 *(u[i, ] + v) * bsxfun('min', u[i, ], v) ^ 2 + u[i, ] * v * (u[i, ] ^ 2 + 3 * u[i, ] * v + v ^2) * bsxfun('min', u[i, ], v) ^ 3 - (1 / 4)*(u[i, ] ^ 3 + 9 * u[i, ] ^ 2  * v + 9 * u[i, ] * v ^ 2 + v ^ 3) * bsxfun('min', u[i, ], v) ^ 4 + (3 / 5) * (u[i, ] ^ 2 + 3 * u[i, ] * v + v ^ 2) * bsxfun('min', u[i, ], v) ^ 5 - (1/2) * (u[i, ] + v) * bsxfun('min', u[i, ], v) ^ 6 + (1 / 7) * bsxfun('min', u[i, ], v) ^ 7
			k[i] = prod(z)
		}
	} else if (strcmp(ker, 'anovabspline')) {

		k = array(0, c(r1, 1))
		
		for (i in 1:r1) {
			z = 0
			for (r in 0:(2*(kpar1+1))) {
				z = z + (-1) ^ r * choose(2 * (kpar1 + 1), r) * (bsxfun('max', array(0, length(v)), u[i, ] - v + kpar1 + 1 - r)) ^ (2 * kpar1 + 1)
			}
			k[i] = prod(z)
		}
	} else {
		cat(paste('CalcKernel: wrong kernel \'', ker, '\'\n'))
	}
	
	return(k)
}
