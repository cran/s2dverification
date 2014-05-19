IniListDims <- function(dims, lenlist) {
  # This function generates a list of arrays where those arrays contain 
  # integers from 1 to various numbers. This list of arrays is used in the 
  # other functions as a list of indices of the elements of the matrices.
  #
  # Args:
  #   dims: The dimensions of a matrix the elements of which we want to 
  #         generate the indices to. For exemple, if the dimensions sent are 
  #         c(3,2,5), the following list of arrays will be generated:
  #         list(c(1:3), c(1:2), c(1:5))
  #   lenlist: lenlist is the length of the list because the list will be 
  #            complemented above length(dims) by arrays of length 1. 
  #            For example, if lenlist is set to 7, the previous list of arrays 
  #            will be extended to:
  #            list(c(1:3), c(1:2), c(1:5), 1, 1, 1, 1)
  #
  # Returns:
  #   A list with lenlist elements, each with arrays with integers from 1 to 
  #   the corresponding number in dims array.
  #
  # History:
  #   1.0  #  2011-04  (V. Guemas, vguemas@ic3.cat)  #  Original code
  
  u <- list()
  for (jdim in 1:length(dims)) {
    u[[jdim]] <- 1:dims[jdim]
  }
  for (jdim in (length(dims) + 1):lenlist) {
    u[[jdim]] <- 1
  }
  #
  #  Outputs
  # ~~~~~~~~~
  #
  out <- u
}
