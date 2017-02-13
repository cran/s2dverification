Enlarge <- function(var, numdims) {
  if (is.numeric(var) || is.logical(var)) {
    if (is.null(dim(var))) {
      dim(var) <- length(var)
    }
    if (numdims >= length(dim(var))) {
      dim(var) <- c(dim(var), rep(1, numdims - length(dim(var))))
      var
    } else {
      stop("'numdims' must be higher or equal to length(dim(var))")
    }
  } else {
    stop("'var' must be a numeric object")
  }
}
