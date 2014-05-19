SelIndices <- function(var, posdim, limits) {
  # This function allows to select a subensemble from a matrix of any 
  # dimensions, providing the dimension along which the user aims at cutting 
  # the matrix and between which indices.
  #
  # Args:
  #   var: A matrix of any rank and any dimensions.
  #   posdim: The dimension along which a submatrice should be selected.
  #   limits: The lower and upper indice of the selection along the (posdim)th 
  #           dimension.
  #
  # Returns:
  #   The sliced matrix.
  #
  # History:
  #   1.0  #  2011-04  (V. Guemas, vguemas@ic3.cat)  #  Original code   

  #
  #  A few security checks
  # ~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimsvar <- dim(var)
  if (is.null(dimsvar)) {
    dimsvar <- length(var)
  }
  if (posdim > length(dimsvar)) {
    stop("posdim does not exist")
  }
  if (length(limits) != 2) {
    stop("Need lower and upper limit")
  } 
  if (dimsvar[posdim] < limits[2]) {
    stop("Check the consistency between limits and var dimensions")
  }
  #
  #  Select the correct indices 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  enlvar <- Enlarge(var, 10)
  u <- IniListDims(dimsvar, 10)
  u[[posdim]] <- limits[1]:limits[2]
  enlvarout <- enlvar[u[[1]], u[[2]], u[[3]], u[[4]], u[[5]], u[[6]], u[[7]],
                      u[[8]], u[[9]], u[[10]]]
  #
  #  Preparing the output matrice
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  if (limits[2] == limits[1]) {
    dimsvar <- dimsvar[-posdim]
  } else {
    dimsvar[posdim] <- limits[2] - limits[1] + 1
  }
  varout <- array(dim = dimsvar)
  varout[] <- enlvarout
  
  #
  #  Output
  # ~~~~~~~~
  #
  varout
}
