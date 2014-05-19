Enlarge <- function(var, numdims) {
  # Extends the number of dimensions of var to numdims (the added dimensions 
  # have length 1).
  #
  # Args:
  #   var: Array to extend.
  #   numdims: Desired number of dimensions.
  #
  # Returns:
  #   Extended array.
  #
  # History:
  #   1.0  #  2011-03  (V. Guemas, vguemas@ic3.cat)  #  Original code

  #
  #  Enlarge the number of dimensions to 20 --> enlvar
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimsvar <- dim(var)
  if (is.null(dimsvar)) {
    dimsvar <- length(var)
  }
  d <- c(dimsvar, array(1, dim = 20))
  enlvar <- array(dim = d[1:20])
  enlvar[, , , , , , , , , , , , , , , , , , , ] <- var
  #
  #  Reduce the number of dimensions to the required one
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  outvar <- array(dim = d[1:numdims])
  outvar[] <- enlvar
  #
  #  Outputs
  # ~~~~~~~~~
  #
  outvar
}
