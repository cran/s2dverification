IniListDims <- function(dims, lenlist) {
  indices <- as.list(rep(1, lenlist))
  for (jdim in 1:min(length(dims), lenlist)) {
    indices[[jdim]] <- 1:dims[jdim]
  }
  indices
}
