MeanListDim <- function(var, dims, narm = TRUE) {
  # Averages the matrix var along a set of dimensions given as argument dims.
  #
  # Args:
  #   var: Matrix to average.
  #   dims: List of dimensions to average along.
  #   narm: Ignore NA values or not. 
  #
  # Returns:
  #   Matrix with the averages applied.
  #
  # History:
  #   1.0  #  2011-04  (V. Guemas, vguemas@ic3.cat)  #  Original code           
  
  #
  whichdims <- sort(dims, decreasing = TRUE)
  for (jdim in 1:length(whichdims)) {
    if (whichdims[jdim] <= length(dim(var))) {
      if (length(dim(var)) > 1) {
        var <- Mean1Dim(var, posdim = whichdims[jdim], narm = narm)
      } else {
        var <- mean(var, na.rm = narm)
      }
    }
  }
  # 
  #  Outputs
  # ~~~~~~~~~
  # 
  outvar <- var
}
