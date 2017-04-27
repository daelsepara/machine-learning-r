# convolution layer, can also add a rectified linear unit layer after convolution
nnet_conv <- function(input, feature, shape = "full") {
	
	# input dimensions
	ix = ncol(input)
	iy = nrow(input)
	
	# filter dimension
	fx = ncol(feature)
	fx = nrow(feature)
	
	# convolution dimensions
	cx = ix + fx - 1;
	cy = iy + fx - 1;
	
	rx = cx;
	ry = cy;
	
	offset = 0;
	
	if (shape == "valid") {
		
		rx = ix - fx + 1;
		ry = iy - fx + 1;
		
		offset = 1;
		
	} else if (shape == "same") {
		
		rx = ix;
		ry = iy;
		
		offset = 1;
	}
	
	if (iy >= fx && ix >= fx) {
		
		result = array(0, c(ry, rx))
		
		for (cj in 1:(cy + 1)) {
			for (ci in 1:(cx + 1)) {
				for (ky in 1:iy) {
					for (kx in 1:ix) {
						if (ci - offset - 1 > 0 && ci - offset - 1 <= rx && cj - offset - 1 > 0 && cj - offset - 1 <= ry && ci - kx > 0 && ci - kx <= fx && cj - ky > 0 && cj - ky <= fx) {
							result[cj - offset - 1, ci - offset - 1] = result[cj - offset - 1, ci - offset - 1] + input[ky, kx] * feature[cj - ky, ci - kx];
						}
					}
				}
			}
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
	cx = ix + fx - 1;
	cy = iy + fy - 1;
	cz = iz + fz - 1;
	
	rx = cx;
	ry = cy;
	rz = cz;
	
	offset = 0;
	
	if (shape == "valid") {
		
		rx = ix - fx + 1;
		ry = iy - fy + 1;
		rz = iz - fz + 1;
		
		offset = 1;
		
	} else if (shape == "same") {
		
		rx = ix;
		ry = iy;
		rz = iz;
		
		offset = 1;
	}
	
	if (iy >= fy && ix >= fx && iz >= fz) {
		
		result = array(0, c(ry, rx, rz))
		
		for (ck in 0:(cz + 1)) {
			for (cj in 1:(cy + 1)) {
				for (ci in 1:(cx + 1)) {
					for (kz in 1:iz) {
						for (ky in 1:iy) {
							for (kx in 1:ix) {
								if ((ci - offset - 1) > 0 && (ci - offset - 1) <= rx && (cj - offset - 1) > 0 && (cj - offset - 1) <= ry && (ck - 1) > 0 && (ck - 1) <= rz && (ci - kx) > 0 && (ci - kx) <= fx && (cj - ky) > 0 && (cj - ky) <= fy && (ck - kz) > 0 && (ck - kz) <= fz) {
									print(paste0(' ci=',ci,' cj=',cj,' ck=',ck,' kx=',kx,' ky=',ky,' kz=',kz,' ci-kx=',ci-kx,' cj-ky=',cj-ky,' ck-kz=',ck-kz))
									result[cj - offset - 1, ci - offset - 1, ck - 1] = result[cj - offset - 1, ci - offset - 1, ck - 1] + input[ky, kx, kz] * feature[cj - ky, ci - kx, ck - kz];
								}
							}
						}
					}
				}
			}
		}
		
		return(drop(result))
	
	} else {
		
		stop('input and filter dimensions are incompatible')
	}
}

# pooling layer
nnet_pool <- function(input, window_, steps_) {
  
  ix = ncol(input)
  iy = nrow(input)
  
  if (ix >= window_ && iy >= window_) {
    
    colseq_ = seq(1, ix, window_)
    rowseq_ = seq(1, iy, window_)
    
    cols_ = length(colseq_)
    rows_ = length(rowseq_)
    
    result = array(0, c(rows_, cols_))
    
    for (y_ in 1:rows_) {
      for(x_ in 1:cols_) {
        
        col_ = colseq_[x_]
        row_ = rowseq_[y_]
        
        px_ = col_ + window_ - 1
        py_ = row_ + window_ - 1
        
        if (px_ > ix) {
          px_ = ix  
        }
        
        if (py_ > iy) {
          py_ = iy  
        }
        
        result[y_, x_] = max(input[row_:py_, col_:px_])
      }
    }
    
    return(result)

  } else {
    
    stop('input and window dimensions are incompatible')   
    
  }
}

# zero-padding function
nnet_pad <- function(input, padsize = 0) {
  
  if (padsize >= 0) {
    
    if (padsize > 0) {
      
      # zero pad columns
      conv_c = cbind(array(0, c(nrow(input), padsize)), input, array(0, c(nrow(input), padsize)))
      
      # zero pad rows
      conv_r = rbind(array(0, c(padsize, ncol(conv_c))), conv_c, array(0, c(padsize, ncol(conv_c))))
      
      return(conv_r)
      
    } else {
      
      return(drop(result))
      
    }
    
  } else {
    
    stop('padsize must be >= 0')
    
  }
}