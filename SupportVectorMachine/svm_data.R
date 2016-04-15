svm_data1 <- function() {

	l = 2
	pts = 30
	N = 9 * pts
	X = numeric(0)
	y = numeric(0)

	for (i in 0:2) {
		for (j in 0:2) {
			X = cbind(X, array(runif(l * pts), c(l, pts)) + array(c(i, j), c(2, 1)) %*% array(1, c(1, pts)))
			if ((i + j) %% 2 == 0) {
				y = cbind(y, array(1, c(1, pts)));
			} else {
				y = cbind(y, array(-1, c(1, pts)));
			}
		}
	}

	X = t(X)
	y = t(y)
	
	return(list('X' = X, 'y' = y))
}

svm_data2 <- function() {

	l = 2
	N = 150

	X = 10 * array(runif(l * N), c(l, N)) - 5
	y = array(0, c(N, 1))

	for (i in 1:N) {
		tt = 0.05 * (X[1, i] ^ 3 + X[1, i] ^ 2 + X[1, i] + 1)
		if (tt > X[2, i]) {
			y[i] = 1
		} else {
			y[i] = -1
		}
	}

	X = t(X)
	
	return(list('X' = X, 'y' = y))
}
