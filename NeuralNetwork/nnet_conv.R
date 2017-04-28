require(pracma)

# convolution layer
nnet_conv <- function(input, feature, shape = "full") {
	
	# input dimensions
	ix = dim(input)[2]
	iy = dim(input)[1]
	
	# filter dimension
	fx = dim(feature)[2]
	fy = dim(feature)[1]
	
	# convolution dimensions
	cx = ix + fx - 1
	cy = iy + fy - 1
	
	if (iy >= fy && ix >= fx) {
		
		result = array(0, c(cy, cx))
		
		for (cj in 1:(cy + 1)) {
			for (ci in 1:(cx + 1)) {
				for (ky in 1:iy) {
					for (kx in 1:ix) {
						if (ci - 1 > 0 && ci - 1 <= cx && cj - 1 > 0 && cj - 1 <= cy && ci - kx > 0 && ci - kx <= fx && cj - ky > 0 && cj - ky <= fy) {
							result[cj - 1, ci - 1] = result[cj - 1, ci - 1] + input[ky, kx] * feature[cj - ky, ci - kx]
						}
					}
				}
			}
		}
		
		if (shape == "valid") {
		
			result = result[fy:(cy - fy + 1), fx:(cx - fx + 1)]
		
		} else if (shape == "same") {
		
			result = result[fy:cy, fx:cx]
		}
		
		return(drop(result))
	
	} else {
		
		stop('input and filter dimensions are incompatible')
	}
}

# convolution layer
nnet_conv3 <- function(input, feature, shape = "full") {
	
	# input dimensions
	ix = dim(input)[2]
	iy = dim(input)[1]
	iz = dim(input)[3]
	
	# filter dimension
	fx = dim(feature)[2]
	fy = dim(feature)[1]
	fz = dim(feature)[3]
	
	# convolution dimensions
	cx = ix + fx - 1
	cy = iy + fy - 1
	cz = iz + fz - 1
	
	if (iy >= fy && ix >= fx && iz >= fz) {
		
		result = array(0, c(cy, cx, cz))
		
		for (ck in 0:(cz + 1)) {
			for (cj in 1:(cy + 1)) {
				for (ci in 1:(cx + 1)) {
					for (kz in 1:iz) {
						for (ky in 1:iy) {
							for (kx in 1:ix) {
								if ((ci - 1) > 0 && (ci - 1) <= cx && (cj - 1) > 0 && (cj - 1) <= cy && (ck - 1) > 0 && (ck - 1) <= cz && (ci - kx) > 0 && (ci - kx) <= fx && (cj - ky) > 0 && (cj - ky) <= fy && (ck - kz) > 0 && (ck - kz) <= fz) {
									result[cj - 1, ci - 1, ck - 1] = result[cj - 1, ci - 1, ck -  1] + input[ky, kx, kz] * feature[cj - ky, ci - kx, ck - kz]
								}
							}
						}
					}
				}
			}
		}
		
		if (shape == "valid") {
		
			result = result[fy:(cy - fy + 1), fx:(cx - fx + 1), fz:(cz - fz + 1)]
		
		} else if (shape == "same") {
		
			result = result[fy:cy, fy:cx, fz:cz]
		}
		
		return(drop(result))
	
	} else {
		
		stop('input and filter dimensions are incompatible')
	}
}

# pooling layer
nnet_pool <- function(img_, window_, steps_) {
  
  ix_ = ncol(img_)
  iy_ = nrow(img_)
  
  if (ix_ >= window_ && iy_ >= window_) {
    
    colseq_ = seq(1, ix_, window_)
    rowseq_ = seq(1, iy_, window_)
    
    cols_ = length(colseq_)
    rows_ = length(rowseq_)
    
    result_ = array(0, c(rows_, cols_))
    
    for (y_ in 1:rows_) {
      for(x_ in 1:cols_) {
        
        col_ = colseq_[x_]
        row_ = rowseq_[y_]
        
        px_ = col_ + window_ - 1
        py_ = row_ + window_ - 1
        
        if (px_ > ix_) {
          px_ = ix_  
        }
        
        if (py_ > iy_) {
          py_ = iy_  
        }
        
        result_[y_, x_] = max(img_[row_:py_, col_:px_])
      }
    }
    
    return(result_)

  } else {
    
    stop('input and window dimensions are incompatible')   
    
  }
}

nnet_expand <- function(A, SZ, scale = 1.0)
{
  if (length(size(A) == length(SZ))) {
    
    return (scale * repmat(A, SZ[1], SZ[2]))
    
  } else {
    
    stop('Length of size vector must equal ndims(A)')
    
  }
}

# zero-padding function
nnet_pad <- function(img_, padsize = 0) {
  
  if (padsize >= 0) {
    
    if (padsize > 0) {
      
      # zero pad columns
      conv_c = cbind(array(0, c(nrow(img_), padsize)), img_, array(0, c(nrow(img_), padsize)))
      
      # zero pad rows
      conv_r = rbind(array(0, c(padsize, ncol(conv_c))), conv_c, array(0, c(padsize, ncol(conv_c))))
      
      return(conv_r)
      
    } else {
      
      return(img_)
      
    }
    
  } else {
    
    stop('padsize must be >= 0')
    
  }
}
