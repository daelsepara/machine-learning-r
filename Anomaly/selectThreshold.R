selectThreshold <- function(y, p) {
# Determine best threshold to use to detect outliers
#
# Inputs:
# y[m] test set (m sample)
# p[m] predicted values for test set X from estimated multivariate Gaussian PDF (using estimated parameters)
#
# Outputs:
#   epsilon: Best prediction threshold to use
#   F1: F1 score on epsilon
  
	epsilon_Best = 0
	F1_Best = 0

	# setup scan range
	stepsize = (max(p) - min(p)) / 1000
	
	for (epsilon in seq(min(p), max(p), by = stepsize)) {

		predictions = p < epsilon

		# determine true/false positives and negatives
		tp = sum((predictions == TRUE) & (y == 1))
		fp = sum((predictions == TRUE) & (y == 0))
		fn = sum((predictions == FALSE) & (y == 1))

		# compute precision and recall
		prec = tp/(tp + fp)
		rec = tp/(tp + fn)

		# compute F1 score based on precision and recall
		F1 =  2*prec*rec/(prec + rec)

		if (!is.nan(F1) && is.numeric(F1) && F1 > F1_Best) {
			F1_Best = F1
			epsilon_Best = epsilon
		}
	}
	
	return(list('F1' = F1_Best, 'epsilon' = epsilon_Best))

}
