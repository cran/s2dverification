#'Averages An Array Along A Dimension
#'
#'Averages the array along the posdim dimension along the user specified 
#'dimension. The user can specify a subset of the dimension to take the mean 
#'along.
#'
#'@param var Matrix to average.
#'@param posdim Dimension to average along.
#'@param narm Ignore NA (TRUE) values or not (FALSE).
#'@param limits Limits to average between. Default is to take the mean along 
#'  the entire dimension.
#'
#'@return Array with one dimension less than the input array, containing 
#'  the average along the posdim dimension.
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2011-04  (V. Guemas, \email{virginie.guemas@@bsc.es})  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@@bsc.es})  -  Formatting to R CRAN
#'@examples
#'a <- array(rnorm(24), dim = c(2, 3, 4))
#'print(a)
#'print(Mean1Dim(a, 2))
#'@export
Mean1Dim <- function(var, posdim, narm = TRUE, limits = NULL) {
  if (is.null(limits)) {
    limits <- c(1, dim(var)[posdim])
  }
  #
  #  Initialisation of the output var with the required dimension length 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  dimsvar <- dim(var)
  outdim <- c()
  if (posdim > 1) {
    outdim <- c(dimsvar[1:(posdim - 1)], outdim)
  }
  if (posdim < length(dimsvar)) {
    outdim <- c(outdim, dimsvar[(posdim + 1):length(dimsvar)])
  }
  tmpvar <- array(0, dim = c(outdim, array(1, dim = (9 - length(outdim)))))
  #
  #  Handling NA
  # ~~~~~~~~~~~~~
  #
  weight <- array(1, dim = dimsvar)
  if (narm) {
    weight[is.na(var)] <- 0
    var[is.na(var)] <- 0 
  }
  outweight <- array(0, dim = c(outdim, array(1, dim = (9 - length(outdim)))))
  #
  #  Average the matrix over the required (posdim)th dimension
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  u <- IniListDims(dimsvar, 10)
  v <- IniListDims(outdim, 9)
  enlvar <- Enlarge(var, 10)
  weight <- Enlarge(weight, 10)
 
  for (jindex in limits[1]:limits[2]) {
    u[[posdim]] <- jindex
    tmpvar[v[[1]], v[[2]], v[[3]], v[[4]], v[[5]], v[[6]], v[[7]], v[[8]],
           v[[9]]] <- tmpvar[v[[1]], v[[2]], v[[3]], v[[4]], v[[5]], v[[6]],
                             v[[7]], v[[8]], 
                             v[[9]]] + enlvar[u[[1]], u[[2]], u[[3]], u[[4]], 
                             u[[5]], u[[6]], u[[7]], u[[8]], u[[9]], u[[10]]]
    outweight[v[[1]], v[[2]], v[[3]], v[[4]], v[[5]], v[[6]], v[[7]], v[[8]],
              v[[9]]] <- outweight[v[[1]], v[[2]], v[[3]], v[[4]], v[[5]],
                                   v[[6]], v[[7]], v[[8]], 
                                   v[[9]]] + weight[u[[1]], u[[2]], u[[3]], 
                                   u[[4]], u[[5]], u[[6]], u[[7]], u[[8]], 
                                   u[[9]], u[[10]]]
  }
  tmpvar <- tmpvar / outweight
  #
  #  Reduce the number of dimensions to the required one
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  
  outvar <- array(dim = outdim)
  outvar[] <- tmpvar
  #
  #  Outputs
  # ~~~~~~~~~
  #
  outvar  
}
