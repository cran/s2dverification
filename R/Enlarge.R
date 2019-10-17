#'Extends The Number Of Dimensions of A Matrix
#'
#'Extends the number of dimensions of var to numdims (the added dimensions 
#'have length 1).
#'
#'@param var Matrix to be extended.
#'@param numdims Output number of dimensions.
#'
#'@return Output number of dimensions. 
#'
#'@keywords datagen
#'@author History:\cr
#'  0.1  -  2011-03  (V. Guemas, \email{virginie.guemas@@ic3.cat})  -  Original code\cr
#'  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@@ic3.cat})  -  Formatting to R CRAN\cr
#'  1.1  -  2015-03  (N. Manubens, \email{nicolau.manubens@@ic3.cat})  -  Improved\cr
#'@examples
#'data <- array(1, c(2, 2, 3))
#'print(dim(Enlarge(data, 5)))
#'@export
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
