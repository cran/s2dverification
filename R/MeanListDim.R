MeanListDim <- function(var, dims, narm = TRUE) {
  apply(var, setdiff(c(1:length(dim(var))), dims), mean, na.rm = narm)
}
