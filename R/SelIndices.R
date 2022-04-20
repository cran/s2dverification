#'Slices A Matrix Along A Dimension
#'
#'This function selects a subset of ensemble members from an array containing 
#'any number of dimensions.
#'
#'@param var An array with any number of dimensions.
#'@param posdim The dimension along which the ensemble subset should be 
#'  selected.
#'@param limits The lower and upper limits for the selection of ensemble 
#'  members along the posdim dimension.
#'
#'@return The subsetted array.
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2011-04  (V. Guemas)  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens)  -  Formatting to CRAN
#'@examples
#'a <- array(rnorm(24), dim = c(2, 3, 4, 1))
#'print(a)
#'print(a[, , 2:3, ])
#'print(dim(a[, , 2:3, ]))
#'print(SelIndices(a, 3, c(2, 3)))
#'print(dim(SelIndices(a, 3, c(2, 3))))
#'
#'@export
SelIndices <- function(var, posdim, limits) {
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
