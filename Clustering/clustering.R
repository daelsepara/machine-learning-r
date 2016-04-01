BSAS <- function(X,theta, qc, ordr) {
#########################################################################
# FUNCTION
#  [bel, repre]=BSAS(X,theta,q,order)
# This function implements the BSAS (Basic Sequential Algorithmic Scheme)
# algorithm. It performs a single pass on the data. If the currently
# considered vector lies at a significant distance (greater than a given
# dissimilarity threshold) from the clusters formed so far a new cluster
# is formed represented by this vector. Otherwise the considered vector
# is assigned to its closest cluster. The results of the algorithm are 
# influenced by the order of presentation of the data.
#
# INPUT ARGUMENTS:
#  X:       lxN matrix, each column of which corresponds to an
#           l-dimensional data vector.
#  theta:   the dissimilarity threshold.
#  qc:       the maximum allowable number of clusters.
#  ordr:   N-dimensional vector containing a permutation of the integers
#           1,2,...,N. The i-th element of this vector specifies the order of
#           presentation of the i-th vector to the algorithm.
#
# OUTPUT ARGUMENTS:
#  bel:     N-dimensional vector whose i-th element contains the
#           cluster label for the i-th data vector.
#  repre:   a matrix, each column of which contains the l-dimensional (mean) 
#           representative of each cluster.
#
# (c) 2010 S. Theodoridis, A. Pikrakis, K. Koutroumbas, D. Cavouras
#
# Converted to R by: SD Separa (2016/03/31)
#
##########################################################################

	l = nrow(X)
	N = ncol(X)
	
	if (length(ordr) == 0) {
		ord = 1:N
	}

	# Cluster determination phase
	n_clust = 1  # no. of clusters
	
	bel = rep(0, N)
	bel[ordr[1]] = n_clust
	
	repre = array(X[, ordr[1]], c(length(X[, ordr[1]]), 1))
	
	for (i in 2:N) {
	   
	   m1 = nrow(repre)
	   m2 = ncol(repre)
	   
	   # Determining the closest cluster representative}}
	   x = array(X[, ordr[i]], c(m1, 1))
	   d = sqrt(apply((repre -  x %*% array(1, c(1, m2)))^2, 2, sum))
	   s1 = min(d)
	   s2 = which.min(d)
	   
	   if ((s1 > theta) && (n_clust < qc)) {
		   n_clust = n_clust + 1
		   bel[ordr[i]] = n_clust
		   repre = cbind(repre, x)
	   } else {
			# Pattern classification phase(*4)}}
		   bel[ordr[i]] = s2
		   repre[, s2] = ((sum(bel == s2) - 1) * repre[, s2] + x)/sum(bel == s2)
	   }
	}
	
	return(list('bel' = bel, 'repre' = repre))
}

SL_step <- function(prox_mat, merge_pair) {
##########################################################################
# FUNCTION
#  prox_new=SL_step(prox_mat, merge_pair)
# This function performs a step of the single link algorithm.
# Specifically, given (a) the distances between clusters of the t-th
# level clustering and (b) the pair of clusters that has been selected
# for merging, the function computes the distances of the newly formed
# cluster from the remaining ones, when the single link algorithm is
# adopted (the distances between the other clusters remains unaltered).
#
# INPUT ARGUMENTS:
#  prox_mat:    NxN dissimilarity matrix for the N vectors of
#               the data set at hand (prox_mat(i,j) is the distance between
#               vectors xi and xj).
#  merge_pair:  2-dimensional vector containing the labels of the
#               clusters that are to be merged.
#
# OUTPUT ARGUMENTS:
#  prox_new:    matrix containing the distances between the clusters of the
#               (t+1)th level clustering.
#
# (c) 2010 S. Theodoridis, A. Pikrakis, K. Koutroumbas, D. Cavouras
#
# Converted to R by: SD Separa (2016/03/31)
#
##########################################################################

	p1 = nrow(prox_mat)
	p2 = ncol(prox_mat)
	
	# handle out-of-bounds errors
	if (merge_pair[2] + 1 > p2) {
		prox_new = cbind(numeric(0), prox_mat[, 1:(merge_pair[2] - 1)])
	} else {
		prox_new = cbind(prox_mat[, 1:(merge_pair[2] - 1)], prox_mat[, (merge_pair[2] + 1):p2])
	}
	
	tt = apply(rbind(prox_new[merge_pair[1], ], prox_new[merge_pair[2], ]), 2, min)
	prox_new[merge_pair[1], ] = tt
	
	# handle out-of-bounds errors
	if (merge_pair[2] + 1 > p1) {
		prox_new = rbind(numeric(0), prox_new[1:(merge_pair[2] - 1), ])
	} else {
		prox_new = rbind(prox_new[1:(merge_pair[2] - 1), ], prox_new[(merge_pair[2] + 1):p1, ])
	}
	
	prox_new[, merge_pair[1]] = tt

	return(prox_new)
}

CL_step <- function(prox_mat, merge_pair) {
##########################################################################
# FUNCTION
#  prox_new=CL_step(prox_mat, merge_pair)
# This function performs a step of the complete link algorithm.
# Specifically, given (a) the distances between clusters of the t-th
# level clustering and (b) the pair of clusters that has been selected
# for merging, the function computes the distances of the newly formed
# cluster from the remaining ones, when the complete link algorithm is
# adopted.
#
# INPUT ARGUMENTS:
#   prox_mat:   NxN dimensional dissimilarity matrix for the N vectors of
#               the data set at hand (prox_mat(i,j) is the distance between
#               vectors xi and xj).
#
#   merge_pair: a 2-dimensional vector containing the labels of the
#               clusters that are to be merged.
#
# OUTPUT ARGUMENTS:
#   prox_new:   a matrix containing the distances between the clusters of the
#               (t+1)th level clustering.
#
# (c) 2010 S. Theodoridis, A. Pikrakis, K. Koutroumbas, D. Cavouras
#
# Converted to R by: SD Separa (2016/03/31)
#
##########################################################################

	p1 = nrow(prox_mat)
	p2 = ncol(prox_mat)
	
	# handle out-of-bounds errors
	if (merge_pair[2] + 1 > p2) {
		prox_new = cbind(numeric(0), prox_mat[, 1:(merge_pair[2] - 1)])
	} else {
		prox_new = cbind(prox_mat[, 1:(merge_pair[2] - 1)], prox_mat[, (merge_pair[2] + 1):p2])
	}
	
	tt = apply(rbind(prox_new[merge_pair[1], ], prox_new[merge_pair[2], ]), 2, max)
	
	prox_new[merge_pair[1], ] = tt
	
	# handle out-of-bounds errors
	if (merge_pair[2] + 1 > p1) {
		prox_new = rbind(numeric(0), prox_new[1:(merge_pair[2] - 1), ])
	} else {
		prox_new = rbind(prox_new[1:(merge_pair[2] - 1), ], prox_new[(merge_pair[2] + 1):p1, ])
	}
	
	prox_new[, merge_pair[1]] = tt
	prox_new[merge_pair[1], merge_pair[1]] = 0
	
	return(prox_new)
}

agglom <- function(prox_mat,code) {
##########################################################################
# FUNCTION
#  [bel, thres]=agglom(prox_mat,code)
# Implements the generalized agglomerative scheme (GAS).
#
# INPUT ARGUMENTS:
#   prox_mat:   NxN dissimilarity matrix for the N vectors of
#               the data set at hand (prox_mat(i,j) is the distance between
#               vectors xi and xj).
#
#   code:       integer indicating the specific clustering algorithm that
#               will be used: "1" stands for single link and "2" for complete
#               link.
#
# OUTPUT ARGUMENTS:
#   bel:        NxN dimensional matrix whose i-th row corresponds to the
#               i-th clustering. The bel(i,j) element of the matrix contains
#               the cluster label for the j-th vector in the i-th clustering.
#               The first row of bel corresponds to the N-cluster clustering,
#               the 2nd row corresponds to the clustering where (N-1)-cluster
#               clustering and, finally, the N-th row corresponds to the
#               single-cluster clustering.
#   thres:      N-dimensional vector containing the dissimilarity levels
#               where each new clustering is formed.
#
# (c) 2010 S. Theodoridis, A. Pikrakis, K. Koutroumbas, D. Cavouras
#
# Converted to R by: SD Separa (2016/03/31)
#
##########################################################################

	N = nrow(prox_mat)
	thres = numeric(0)
	bel = array(0, c(N, N))
	bel[1, ] = 1:N

	for (i in 2:N) {
		p_min = prox_mat + 10 ^ 10 * (prox_mat == 0)
		q1 = apply(p_min, 1, min)
		q2 = apply(p_min, 1, which.min)
		r1 = min(q1)
		r2 = which.min(q1)
		merge_pair= c(min(q2[r2], r2), max(q2[r2], r2))
		thres = c(thres, r1)
		temp = bel[i - 1, ]
		temp[which(temp == merge_pair[2])] = merge_pair[1]
		temp = temp * (temp <= merge_pair[2])+(temp - 1)*(temp > merge_pair[2])
		bel[i, ] = temp
		
		if (code == 1) {
			prox_mat = SL_step(prox_mat, merge_pair)
		} else if(code == 2) {
			prox_mat = CL_step(prox_mat, merge_pair)
		}
	}

	cluster = array(list(), c(N, N))
	
	# Making the dendrogram
	for (i in 1:N) {
		for (j in 1:N) {
			dend = which(bel[i, ] == j)
			if (length(dend) > 0) {
				cluster[[i,j]] = dend
			}
		}
	}
	
	return(list('bel' = bel, 'thres' = thres, 'cluster' = cluster))
}

dendrogram_cut <- function(bel, prox_mat) {
##########################################################################
# FUNCTION
#  [lambda,cut_point_tot,hist_cut] = dendrogram_cut(bel,prox_mat)
# This function applies the procedure discussed in the book for
# determining the clusterings of a clustering hierarchy that best fit the
# underlying clustering structure of the data set at hand. Also, it
# provides a histogram with the frequency of selection of a clustering in
# the clustering hierarchy, as the lambda parameter varies.
#
# INPUT ARGUMENTS:
#  prox_mat:    NxN dissimilarity matrix for the N vectors
#               of the data set at hand (prox_mat(i,j) is the distance between
#               vectors xi and xj).
#
#  bel:         NxN matrix whose i-th row corresponds to the
#               i-th clustering. The bel(i,j) element of the matrix contains
#               the cluster label for the j-th vector in the i-th clustering.
#               The first row of bel corresponds to the N-cluster clustering,
#               the 2nd row corresponds to the (N-1)-cluster clustering and,
#               finally, the N-th row corresponds to the single-cluster
#               clustering.
#
# OUTPUT ARGUMENTS:
#  lambda:      a vector of the values of the lambda parameter for which a
#               clustering other than the 1-cluster and the N-cluster
#               clusterings is obtained.
#
#  cut_point_tot: the index of the clustering selected for a given value
#                 of lambda
#
#  hist_cut: a vector whose t-th component contains the number of times
#           the t-th clustering has been selected (1-cluster and N-cluster
#           clusterings are excluded).
#
#  NOTE:   The function also plots the corresponding histogram.
#
# (c) 2010 S. Theodoridis, A. Pikrakis, K. Koutroumbas, D. Cavouras
#
# Converted to R by: SD Separa (2016/03/31)
#
# Changes:
#	Removed fighandle
#
##########################################################################

	N = ncol(prox_mat)

	# Put all distances to a vector form
	dist_vec = numeric(0)
	
	for (i in 1:(N-1)) {
		dist_vec = c(dist_vec, prox_mat[i, (i + 1):N])
	}

	cut_point_tot = numeric(0)

	# The mean of non-zero distances
	mean_dist = mean(dist_vec)
	
	# The standard deviation of non-zero distances
	std_dist = sd(dist_vec)

	lambda_min = -3
	lambda_max = 3
	lambda_step = 0.01

	for (lambda in seq(lambda_min, lambda_max, by = lambda_step)) { # The parameter involved in the definition of the threshold
		
		theta = mean_dist + lambda * std_dist # The external threshold theta
		
		cut_point = N  # The level where the dendrogram will be cut
		
		clust_el_old = array(1, c(1, N))
		
		id = array(0, c(1, N))
		h = array(0, c(1, N))
		
		for (i in 2:N) {
		
			clust_el = array(0, c(1, N))
			
			for (j in 1:N) { # for each cluster determine the number of elements
				clust_el[j] = sum(bel[i, ] == j)
			}
			
			max_op = (clust_el-clust_el_old) * (clust_el > 0)
			
			# Identify the newly formed cluster
			vali = max(max_op) 
			id[i] = which.max(max_op)
			
			list_ = which(bel[i, ] == id[i]) # Identifying the vectors that belong to the newly formed cluster.
			
			# Determination of the h parameter of the newly formed cluster
			dist_list = numeric(0)
			
			for (k in 1:length(list_)) {
				for (qi in (k+1):length(list_)) {
					if (qi <= length(list_)) {
						dist_list = c(dist_list, prox_mat[list_[k], list_[qi]])
					}
				}
			}
			
			h[id[i]] = max(dist_list)   # median(dist_list)
			
			if (h[id[i]] > theta) {
				cut_point = i - 1
				break
			}
			
			clust_el_old = clust_el
		}
		
		cut_point_tot = c(cut_point_tot, cut_point)
	}

	all_clust = sum(cut_point_tot == 1)
	one_clust = sum(cut_point_tot == N)

	le = length(cut_point_tot)
	temp = cut_point_tot[(all_clust + 1):(le - one_clust)]
	cut_point_tot = temp
	lambda = seq(lambda_min, lambda_max, by = lambda_step)
	lambda = lambda[(all_clust + 1):(le - one_clust)]

	qwe = numeric(0)
	hist_cut = numeric(0)
	
	for (i in 2:(N-1)) {
		tt = sum(cut_point_tot == i)
		hist_cut = c(hist_cut, tt)
		qwe = c(qwe, tt/601)
	}
	
	te = 1:length(qwe)
	x_axis = N - te
	
	plot(x_axis, qwe, pch = 8, col = 'red', xlim = c(1, N), ylim = c(0, 0.4), xlab = 'Clusterings', ylab = 'Frequency')
	
	for (i in 1:length(qwe)) {
		lines(c(x_axis[i], x_axis[i]), c(0, qwe[i]))
	}

	return(list('lambda' = lambda, 'cut_point_tot' = cut_point_tot, 'hist_cut' = hist_cut))
}

test_clustering <- function() {
	
	X = array(c(2, 5, 6, 4, 5, 3, 2, 2, 1, 4, 5, 4, 3, 3, 2, 3, 2, 4, 8, 2, 9, 2, 10, 2, 11, 2, 10, 3, 9 ,1), c(2, 15))

	# test BSAS
	ordr = array(c(8, 6, 11, 1, 5, 2, 3, 4, 7, 10, 9, 12, 13, 14, 15), c(1, 15))
	stopifnot(BSAS(X, theta = 2.5, qc = 15, ordr)$bel == c(1, 2, 2, 1, 1, 2, 1, 1, 1, 3, 3, 3, 3, 3, 3))
	stopifnot(BSAS(X, theta = 1.4, qc = 15, ordr)$bel == c(4, 2, 2, 1, 5, 2, 1, 1, 4, 3, 3, 6, 6, 6, 3))
	stopifnot(BSAS(X, theta = 1.4, qc = 2, ordr)$bel == c(1, 2, 2, 1, 1, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2))
	ordr = array(c(7, 3, 1, 5, 9, 6, 8, 4, 2, 10, 15, 13, 14, 11, 12), c(1, 15))
	stopifnot(BSAS(X, theta = 2.5, qc = 15, ordr)$bel == c(2, 1, 1, 2, 2, 1, 1, 2, 2, 3, 3, 4, 4, 4, 3))
	cat('all BSAS tests passed\n')
	
	# test GAS
	D = rbind(c(0, 1, 4, 20, 22, 23), c(1, 0, 3, 22, 24, 25), c(4, 3, 0, 23, 25, 26), c(20, 22, 23, 0, 3.5, 3.6), c(22, 24, 25, 3.5, 0, 3.6), c(23, 25, 26, 3.6, 3.7, 0))
	result = agglom(D, 1)
	stopifnot(result$thres == c(1.0, 3.0, 3.5, 3.6, 20.0))
	cat('all GAS tests passed\n')
}

