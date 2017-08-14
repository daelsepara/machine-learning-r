require(pracma)

# convolution layer
nnet_conv <- function(input, feature, shape = "full") {
	
  # input dimensions
  ix = ncol(input)
  iy = nrow(input)
  
  # filter dimension
  fx = ncol(feature)
  fy = nrow(feature)
  
  # convolution dimensions
  cx = ix + fx - 1
  cy = iy + fy - 1
  
  if (iy >= fy && ix >= fx) {
    
    result = array(0, c(cy, cx, iz))
    
    for (cj in 2:(cy + 1)) {
      for (ci in 2:(cx + 1)) {
        for (ky in 1:iy) {
          if (cj - ky > 0 && cj - ky <= fy) {
            for (kx in 1:ix) {
              if (ci - kx > 0 && ci - kx <= fx) {
                result[cj - 1, ci - 1] = result[cj - 1, ci - 1] + input[ky, kx] * feature[cj - ky, ci - kx]
              }
            }
          }
        }
      }
    }
    
    if (shape == "valid") {
      
      result = result[fy:(cy - fy + 1), fx:(cx - fx + 1)]
      
    } else if (shape == "same") {
      
      dx = (cx - ix)/2
      dy = (cy - iy)/2
      
      result = result[(1 + ceil(dy)):(cy - floor(dy)), (1 + ceil(dx)):(cx - floor(dx))]
    }
    
    return(drop(result))
    
  } else {
    
    stop('input and filter dimensions are incompatible')
  }
}

# convolution layer
nnet_conv3 <- function(input, feature, shape = "full") {
	
  # input dimensions
  ix = ncol(input)
  iy = nrow(input)
  iz = dim(input)[3]
  
  # filter dimension
  fx = ncol(feature)
  fy = nrow(feature)
  fz = dim(feature)[3]
  
  # convolution dimensions
  cx = ix + fx - 1
  cy = iy + fy - 1
  cz = iz + fz - 1
  
  if (iy >= fy && ix >= fx && iz >= fz) {
    
    result = array(0, c(cy, cx, cz))
    
    for (ck in 2:(cz + 1)) {
      for (cj in 2:(cy + 1)) {
        for (ci in 2:(cx + 1)) {
          for (kz in 1:iz) {
            if (ck - kz > 0 && ck - kz <= fz) {
              for (ky in 1:iy) {
                if (cj - ky > 0 && cj - ky <= fy)
                {
                  for (kx in 1:ix) {
                    if (ci - kx > 0 && ci - kx <= fx) {
                      result[cj - 1, ci - 1, ck - 1] = result[cj - 1, ci - 1, ck -  1] + input[ky, kx, kz] * feature[cj - ky, ci - kx, ck - kz]
                    }
                  }
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
      
      dx = (cx - ix)/2
      dy = (cy - iy)/2
      dz = (cz - iz)/2
      
      result = result[(1 + ceil(dy)):(cy - floor(dy)), (1 + ceil(dx)):(cx - floor(dx)), (1 + ceil(dz)):(cz - floor(dz))]
    }
    
    return(drop(result))
    
  } else {
    
    stop('input and filter dimensions are incompatible')
  }
}

# pooling layer
nnet_pool <- function(input, feature, steps_) {
  
  ix = ncol(input)
  iy = nrow(input)
  
  if (ix >= feature && iy >= feature) {
    
    colseq = seq(1, ix, feature)
    rowseq = seq(1, iy, feature)
    
    cols = length(colseq)
    rows = length(rowseq)
    
    result = array(0, c(rows, cols))
    
    for (y in 1:rows) {
      for(x in 1:cols) {
        
        col = colseq[x]
        row = rowseq[y]
        
        px = col + feature - 1
        py = row + feature - 1
        
        if (px > ix) {
          px = ix  
        }
        
        if (py > iy) {
          py = iy  
        }
        
        result[y, x] = max(input[row:py, col:px])
      }
    }
    
    return(result)

  } else {
    
    stop('input and window dimensions are incompatible')   
    
  }
}

nnet_expand <- function(A, SZ, scale = 1.0)
{
  if (length(dim(A) == length(SZ))) {
    
    return (scale * repmat(A, SZ[1], SZ[2]))
    
  } else {
    
    stop('Length of size vector must equal ndims(A)')
    
  }
}

# zero-padding function
nnet_pad <- function(input, padsize = 0) {
  
  if (padsize >= 0) {
    
    if (padsize > 0) {
      
      # zero pad columns
      conv_c = cbind(array(0, c(dim(input)[1], padsize)), input, array(0, c(dim(input)[1], padsize)))
      
      # zero pad rows
      conv_r = rbind(array(0, c(padsize, ncol(conv_c))), conv_c, array(0, c(padsize, ncol(conv_c))))
      
      return(conv_r)
      
    } else {
      
      return(input)
      
    }
    
  } else {
    
    stop('padsize must be >= 0')
    
  }
}
