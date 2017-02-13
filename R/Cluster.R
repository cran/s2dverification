Cluster <- function(var, weights, nclusters = NULL, 
                     index = 'sdindex', posdates = 1) {
  # Check var
  if (!is.numeric(var) || !is.array(var)) {
    stop("The parameter 'var' must be a numeric array.")
  }

  # Check weights
  if (!is.numeric(weights)) {
    stop("The parameter 'weights' must be a numeric array or vector.")
  }
  if (is.null(dim(weights))) {
    dim(weights) <- length(weights)
  }

  # Check posdates
  if (!is.numeric(posdates) || length(posdates) != 1) {
    stop("The parameter 'posdates' must be a single numeric value.")
  }
  posdates <- round(posdates)

  # Consistency checks
  if (posdates < 1 || posdates > length(dim(var))) {
    stop("The parameter 'posdates' must point to one of the dimensions in 'var'.")
  }
  if (!identical(dim(var)[-posdates], dim(weights))) {
    stop("The parameter 'weights' must have the same dimensions as 'var' without the time dimension.")
  }

  # Check nclusters
  if (!is.null(nclusters)) {
    if (!is.numeric(nclusters) || length(nclusters) != 1) {
      stop("The parameter 'nclusters' must be a single numeric value.") 
    }
    nclusters <- round(nclusters)
  }

  # Check index
  if (!is.character(index)) {
    stop("The parameter 'index' should be a vector of character strings accepted as 'index' by the function NbClust::NbClust.")
  }

  # Reordering dimensions if needed
  if (posdates != 1) {
    new_order <- (1:length(dim(var)))[-posdates]
    new_order <- c(posdates, new_order)
    var <- aperm(var, new_order)
  }
  original_dims <- dim(var)
  dim(var) <- c(dim(var)[1], prod(dim(var)[-1]))
  dim(weights) <- length(weights)

  for (i in 1:dim(var)[1]) {
    var[i, ] <- var[i, ] * weights
  }

  if (!is.null(nclusters)) {
    kmeans.results <- kmeans(var, centers = nclusters, iter.max = 300, nstart = 30) 
  } else {
    pdf(file = NULL)
    nbclust.results <- NbClust::NbClust(var, distance = 'euclidean', min.nc = 2, 
                                        max.nc = 20, method = 'kmeans', index = index)
    dev.off()

    if (index == 'all' || index == 'alllong') {
      kmc  <- hist(nbclust.results$Best.nc[1, ], breaks = seq(0, 20), plot = FALSE)$counts
      kmc1 <- which(kmc == max(kmc))
    } else {
      kmc1 <- nbclust.results$Best.nc[1]
    }

    kmeans.results <- kmeans(var, centers = kmc1, iter.max = 300, nstart = 30)
  }
  invisible(kmeans.results)
}
