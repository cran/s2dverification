#'Creates A List Of Integer Ranges
#'
#'This function generates a list of arrays containing integers greater than or 
#'equal to 1. This list of arrays is used in other functions as a list of 
#'indices of the elements of the matrices.
#'
#'@param dims The dimensions of a matrix for which we need the possible 
#'  indices for each dimension. For exemple, if the dimensions sent are 
#'  c(3,2,5), the following list of arrays will be generated:\cr
#'  list(c(1:3), c(1:2), c(1:5)).
#'@param lenlist 'lenlist' is the length of the list because the list will 
#'  be complemented above length(dims) by arrays of length 1.\cr
#'  For example, if lenlist is set to 7, the previous list of arrays will be 
#'  extended to:\cr
#'  list(c(1:3), c(1:2), c(1:5), 1, 1, 1, 1).
#'
#'@return A list with lenlist elements, each with arrays with integers from 1 
#'  to the numbers in dims array and with only 1 for the dimensions above 
#'  length(dims).
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2011-04  (V. Guemas)  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens)  -  Formatting to R CRAN\cr
#'1.1  -  2015-03  (N. Manubens)  -  Improved
#'@examples
#'indices <- IniListDims(c(2, 2, 4, 3), 6)
#'print(indices)
#'@export
IniListDims <- function(dims, lenlist) {
  indices <- as.list(rep(1, lenlist))
  for (jdim in 1:min(length(dims), lenlist)) {
    indices[[jdim]] <- 1:dims[jdim]
  }
  indices
}
