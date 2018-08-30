fast_factorial <- function(x) {
# faster computation of factorials
  
  if (x %% 2) {
    res = x
    x = x - 1
  } else {
    res = 1
  }
  
  i = 0
  
  while (x > 0) {
    i = i + x
    res = res * i
    x = x - 2
  }
  
  return(res)
}

unpack <- function(X) {
# unpacks a list X = list('variable1' = value1, 'variable2' = value2, ...)
# into the objects variable1, variable2, ... whose values are value1, value2, ...
  
  if (is.list(X)) {
    # assign to variables in the parent/calling frame
    for (i in 1:length(X)) {
      assign(names(X)[i], X[[i]], envir = parent.frame())
    }
  }
}
