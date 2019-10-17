#'K-means Clustering
#'
#'This function computes cluster centers and their time series of occurrences, 
#'with the K-means clustering method using Euclidean distance, of an array of 
#'input data with any number of dimensions, one of them (the 'posdates'th) 
#'corresponding to time. By default the first dimension is expected to 
#'correspond to time. Specifically, it partitions the array along time axis in 
#'K groups or clusters in which each space vector/array belongs to (i.e., is a 
#'member of) the cluster with the nearest center or centroid. This function 
#'relies on the NbClust package (Charrad et al., 2014 JSS). 
#'
#'@param var An array with any number of dimensions, one of them (the 
#'  'posdates'th) corresponding to time with either area-averages over a 
#'  series of domains or the grid points for any sptial grid structure (x), 
#'  (y), (z), (x,y), (x,y,z), (y,z), ...
#'@param weights A vector/array of multiplicative weights based on the areas 
#'  covering each domain/region or grid-cell of var; the dimensions of weights 
#'  vector must be equal to the dimensions of 'var' without the 
#'  'posdates'th dimension.
#'@param nclusters This is positive integer K that must be bigger than 1. 
#'  K is the number of clusters to be computed, or K initial cluster centers 
#'  to be used in the method. Default is NULL and then user has to specify 
#'  which index from NbClust and the associated criteria for selecting the 
#'  optimal number of clusters will be used for K-means clustering of var.
#'@param index A validity index from NbClust package that can be used to 
#'  determine optimal K if K is not specified as positive integer bigger than 
#'  1 or initial/seed cluster centers in nclusters. 'sdindex' is deafult 
#'  (Halkidi et al. 2001, JIIS). Other indices also available in NBClust are 
#'  "kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", 
#'  "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2", 
#'  "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", 
#'  "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", and "sdbw". 
#'  One can also use all of them with the option 'alllong' or almost all indices 
#'  except gap, gamma, gplus and tau with 'all', when the optimal number of 
#'  clusters K is detremined by the majority rule (the maximum of histogram of 
#'  the results of all indices with finite solutions). Use of some indices on 
#'  a big and/or unstructured dataset can be computationally intense and/or 
#'  could lead to numerical singularity. 
#'@param posdates The index of the dimension that corresponds to time in the 
#'  provided array in the parameter 'var', the first by default.
#'
#'@return
#'\item{cluster}{
#'  A vector (time series) of integers indicating the occurrence 
#'  of a cluster, i.e., when 'certain data member in time is allocated to a 
#'  specific cluster (e.g., 2 1 3 1 1 1 ..).
#'}
#'\item{centers}{
#'  A matrix of cluster centres or centroids (e.g. 
#'  [1:K, 1:spatial degrees of freedom]).
#'}
#'\item{totss}{
#'  The total sum of squares.
#'}
#'\item{withinss}{
#'  A vector of within-cluster sum of squares, one component 
#'  per cluster.
#'}
#'\item{tot.withinss}{
#'  Total within-cluster sum of squares, 
#'  i.e., sum(withinss).
#'}
#'\item{betweenss}{
#'  The between-cluster sum of squares, i.e. totss-tot.withinss.
#'}
#'\item{size}{
#'The number of points in each cluster.
#'}
#'
#'@references
#'Wilks, 2011, Statistical Methods in the Atmospheric Sciences, 3rd ed., Elsevire, pp 676.
#'@keywords datagen
#'@author History:\cr
#'  1.0 # 2014-10 (N.S. Fuckar, \email{neven.fuckar@@bsc.es}) - Original code
#'@examples
#'# Generating synthetic data
#'a1 <- array(dim = c(200, 4))
#'mean1 <- 0
#'sd1 <- 0.3 
#'
#'c0 <- seq(1, 200)
#'c1 <- sort(sample(x = 1:200, size = sample(x = 50:150, size = 1), replace = FALSE))
#'x1 <- c(1, 1, 1, 1)
#'for (i1 in c1) {
#'  a1[i1, ] <- x1 + rnorm(4, mean = mean1, sd = sd1)
#'}
#'
#'c1p5 <- c0[!(c0 \%in\% c1)]
#'c2 <- c1p5[seq(1, length(c1p5), 2)] 
#'x2 <- c(2, 2, 4, 4)
#'for (i2 in c2) {
#'  a1[i2, ] <- x2 + rnorm(4, mean = mean1, sd = sd1)
#'}
#'
#'c3 <- c1p5[seq(2, length(c1p5), 2)]
#'x3 <- c(3, 3, 1, 1)
#'for (i3 in c3) {
#'  a1[i3, ] <- x3 + rnorm(4, mean = mean1, sd = sd1)
#'}
#'
#'# Computing the clusters
#'res1 <- Cluster(var = a1, weights = array(1, dim = dim(a1)[2]), nclusters = 3)
#'print(res1$cluster)
#'print(res1$centers)
#'
#'res2 <- Cluster(var = a1, weights = array(1, dim = dim(a1)[2]))
#'print(res2$cluster)
#'print(res2$centers)
#'@import NbClust
#'@importFrom stats kmeans
#'@importFrom grDevices pdf dev.off 
#'@export
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
