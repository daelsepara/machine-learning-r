# convolution layer, can also add a rectified linear unit layer after convolution
nnet_conv <- function(img_, filter_, shape = "full") {
	
	# input dimensions
	ix_ = ncol(img_)
	iy_ = nrow(img_)
	
	# filter dimension
	fx_ = ncol(filter_)
	fy_ = nrow(filter_)
	
	# convolution dimensions
	cx_ = ix_ + fx_ - 1;
	cy_ = iy_ + fy_ - 1;
	
	rx_ = cx_;
	ry_ = cy_;
	
	offset = 0;
	
	if (shape == "valid") {
		
		rx_ = ix_ - fx_ + 1;
		ry_ = iy_ - fy_ + 1;
		
		offset = 1;
		
	} else if (shape == "same") {
		
		rx_ = ix_;
		ry_ = iy_;
		
		offset = 1;
	}
	
	if (iy_ >= fy_ && ix_ >= fx_) {
		
		result_ = array(0, c(ry_, rx_))
		
		for (cj in 1:(cy_ + 1)) {
			for (ci in 1:(cx_ + 1)) {
				for (ky in 1:iy_) {
					for (kx in 1:ix_) {
						if (ci - offset - 1 > 0 && ci - offset - 1 <= rx_ && cj - offset - 1 > 0 && cj - offset - 1 <= ry_ && ci - kx > 0 && ci - kx <= fx_ && cj - ky > 0 && cj - ky <= fy_) {
							result_[cj - offset - 1, ci - offset - 1] = result_[cj - offset - 1, ci - offset - 1] + img_[ky, kx] * filter_[cj - ky, ci - kx];
						}
					}
				}
			}
		}
		
		return(result_)
	
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