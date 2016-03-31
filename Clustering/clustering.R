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
