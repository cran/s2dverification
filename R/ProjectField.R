ProjectField <- function(ano, eof, mode = 1) {
  # Checking ano
  if (!is.numeric(ano) || !is.array(ano)) {
    stop("Parameter 'ano' must be a numeric array.")
  }
  if (length(dim(ano)) != 6) {
    stop("'ano' must have dimensions c(n. data sets, n. members, n. start dates, n. forecast time steps, n. latitudes, n. longitudes).")
  }

  # Checking eof
  print_error <- FALSE
  if (!is.list(eof)) {
    print_error <- TRUE
  } else {
    if (!all(c('EOFs', 'wght') %in% names(eof))) {
      print_error <- TRUE
    } else {
      if (!is.numeric(eof$EOFs) || !is.array(eof$EOFs)) {
        stop("The component 'EOFs' of parameter 'eof' must be a numeric array.")
      }
      if (length(dim(eof$EOFs)) != 3) {
        stop("The component 'EOFs' of parameter 'eof' must have dimensions c(n. modes, n. latitudes, n. longitudes).")
      }
      if (!is.numeric(eof$wght) || !is.array(eof$wght)) {
        stop("The component 'wght' of parameter 'eof' must be a numeric array.")
      }
      if (length(dim(eof$wght)) != 2) {
        stop("The component 'wght' of parameter 'eof' must have dimensions c(n. latitudes, n. longitudes).")
      }
    }
  }
  if (print_error) {
    stop("Parameter 'eof' must be a list with the components 'EOFs' and 'whgt' as output from EOF().")
  }

  # Checking mode
  if (!is.numeric(mode)) {
    stop("Parameter 'mode' must be numeric.")
  }
  mode <- round(mode)
  if (mode > dim(eof$EOFs)[1]) {
    stop("Parameter 'mode' is greater than the number of available modes in 'eof'.")
  }

  nlon <- dim(ano)[6]
  nlat <- dim(ano)[5]
  nyr <- dim(ano)[3]
  lt <- dim(ano)[4]
  nfc <- dim(ano)[2]
  nmod <- dim(ano)[1]
  
  # Security checks.
  if (dim(eof$EOFs)[2] != nlat || dim(eof$wght)[1] != nlat) {
    stop("Inconsistent number of latitudes between parameter 'eof' and input field")
  }
  if (dim(eof$EOFs)[3] != nlon || dim(eof$wght)[2] != nlon) {
    stop("Inconsistent number of longitudes between parameter 'eof' and input field")
  }
  
  # Initialization of pc.ver.
  pc.ver <- array(NA, dim = c(nmod, nfc, nyr, lt))
  
  # Weigths
  e.1 <- eof$EOFs[mode, , ] * eof$wght
  ano <- ano * InsertDim(InsertDim(InsertDim(InsertDim(eof$wght, 
    1, nmod), 2, nfc), 3, nyr), 4, lt)
  
  for (imo in 1:nmod) {
    for (im in 1:nfc) {
      for (iy in 1:nyr) {
        for (il in 1:lt) {
          if (!all(is.na(ano[imo, im, iy, il, , ]))) {
          pc.ver[imo, im, iy, il] <- sum(ano[imo, im, 
            iy, il, , ] * e.1, na.rm = TRUE)
          }
        }
      }
    }
  }
  
  return(pc.ver)
}

