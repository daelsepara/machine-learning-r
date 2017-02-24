# convolution layer, can also add a rectified linear unit layer after convolution
nnet_conv <- function(img_, filter_, rectify = FALSE, normalize = TRUE) {
  
  ix_ = ncol(img_)
  iy_ = nrow(img_)
  
  fx_ = ncol(filter_)
  fy_ = nrow(filter_)
  
  if (iy_ >= fy_ && ix_ >= fx_) {

    rows_ = iy_ - fy_ + 1
    cols_ = ix_ - fx_ + 1
    
    result_ = array(0, c(rows_, cols_))
    
    if (normalize) {
      norm_ =  (fx_ * fy_)
    } else {
      norm_ = 1
    }
    
    for (y_ in 1:rows_) {
      for (x_ in 1:cols_) {
        
        result_[y_, x_] = sum(img_[y_:(y_ + fy_ - 1), x_:(x_ + fx_ - 1)] * filter_[1:fy_, 1:fx_]) / norm_
        
      }
    }
    
    # check if RLU's are needed
    if (rectify) {
      
      result_[which(result_ < 0)] = 0
      
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