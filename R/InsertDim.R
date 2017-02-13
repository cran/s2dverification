InsertDim <- function(var, posdim, lendim) {
  if (is.numeric(var) || is.logical(var)) {
    dimsvar <- dim(var)
    if (is.null(dimsvar)) {
      dimsvar <- length(var)
    }
    outdim <- lendim
    if (posdim <= (length(dimsvar) + 1)) {
      if (posdim > 1) {
        outdim <- c(dimsvar[1:(posdim - 1)], outdim)
      }
      if (posdim <= length(dimsvar)) {
        outdim <- c(outdim, dimsvar[posdim:length(dimsvar)])
      }
      outvar <- array(dim = c(outdim, rep(1, 10 - length(outdim))))
      #
      #  Duplicate the matrix along the required (posdim)th dimension
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #
      u <- IniListDims(outdim, 10)
      for (jindex in 1:lendim) {
        u[[posdim]] <- jindex
        outvar[u[[1]], u[[2]], u[[3]], u[[4]], u[[5]], u[[6]], u[[7]], u[[8]], u[[9]], u[[10]]] <- var
      }
      #
      #  Outputs
      # ~~~~~~~~~
      dim(outvar) <- outdim
      outvar
    } else {
      stop("'posdim' must be smaller or equal to the number of dimensions of 'var' plus 1")
    }
  } else {
    stop("'var' must be a numeric object.")
  }
}
