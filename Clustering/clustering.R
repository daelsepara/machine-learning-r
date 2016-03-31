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

test_clustering <- function() {
	
	X = array(c(2, 5, 6, 4, 5, 3, 2, 2, 1, 4, 5, 4, 3, 3, 2, 3, 2, 4, 8, 2, 9, 2, 10, 2, 11, 2, 10, 3, 9 ,1), c(2, 15))

	ordr = array(c(8, 6, 11, 1, 5, 2, 3, 4, 7, 10, 9, 12, 13, 14, 15), c(1, 15))
	stopifnot(BSAS(X, theta = 2.5, qc = 15, ordr)$bel == c(1, 2, 2, 1, 1, 2, 1, 1, 1, 3, 3, 3, 3, 3, 3))
	stopifnot(BSAS(X, theta = 1.4, qc = 15, ordr)$bel == c(4, 2, 2, 1, 5, 2, 1, 1, 4, 3, 3, 6, 6, 6, 3))
	stopifnot(BSAS(X, theta = 1.4, qc = 2, ordr)$bel == c(1, 2, 2, 1, 1, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2))
	ordr = array(c(7, 3, 1, 5, 9, 6, 8, 4, 2, 10, 15, 13, 14, 11, 12), c(1, 15))
	stopifnot(BSAS(X, theta = 2.5, qc = 15, ordr)$bel == c(2, 1, 1, 2, 2, 1, 1, 2, 2, 3, 3, 4, 4, 4, 3))
	cat('all BSAS tests passed\n')
}
