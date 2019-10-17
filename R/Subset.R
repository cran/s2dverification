#'Subset a Data Array
#'
#'This function allows to subset (i.e. slice, take a chunk of) an array, in a 
#'similar way as done in the function \code{take()} in the package plyr. There
#'are two main inprovements:\cr\cr
#'The input array can have dimension names, either 
#'in \code{names(dim(x))} or in the attribute 'dimensions', and the dimensions 
#'to subset along can be specified via the parameter \code{along} either with 
#'integer indices or either by their name.\cr\cr
#'There are additional ways to adjust which dimensions are dropped in the 
#'resulting array: either to drop all, to drop none, to drop only the ones that 
#'have been sliced or to drop only the ones that have not been sliced.\cr\cr
#'If an array is provided without dimension names, dimension names taken from 
#'the parameter \code{dim_names} will be added to the array.
#'
#'@param x A multidimensional array to be sliced. It can have dimension names 
#'  either in \code{names(dim(x))} or either in the attribute 'dimensions'.
#'@param along Vector with references to the dimensions to take the subset 
#'  from: either integers or dimension names.
#'@param indices List of indices to take from each dimension specified in 
#'  'along'. If a single dimension is specified in 'along' the indices can be 
#'  directly provided as a single integer or as a vector.
#'@param drop Whether to drop all the dimensions of length 1 in the resulting 
#'  array, none, only those that are specified in 'along', or only those that 
#'  are not specified in 'along'. The possible values are, respectively: 'all' 
#'  or TRUE, 'none' or FALSE, 'selected', and 'non-selected'.
#'
#'@keywords datagen
#'@examples
#'subset <- Subset(sampleMap$mod, c('dataset', 'sdate', 'ftime'), 
#'                 list(1, 1, 1), drop = 'selected')
#'PlotLayout(PlotEquiMap, c('lat', 'lon'), subset, 
#'           sampleMap$lon, sampleMap$lat, 
#'           titles = paste('Member', 1:3))
#'
#'@export
Subset <- function(x, along, indices, drop = FALSE) {
  # Check x
  if (!is.array(x)) {
    stop("Input array 'x' must be a numeric array.")
  }

  # Take the input array dimension names
  dim_names <- attr(x, 'dimensions')
  if (!is.character(dim_names)) {
    dim_names <- names(dim(x))
  }
  if (!is.character(dim_names)) {
    if (any(sapply(along, is.character))) {
      stop("The input array 'x' doesn't have labels for the dimensions but the parameter 'along' contains dimension names.")
    }
  }

  # Check along
  if (any(sapply(along, function(x) !is.numeric(x) && !is.character(x)))) {
    stop("All provided dimension indices in 'along' must be integers or character strings.")
  }
  if (any(sapply(along, is.character))) {
    req_dimnames <- along[which(sapply(along, is.character))]
    if (length(unique(req_dimnames)) < length(req_dimnames)) {
      stop("The parameter 'along' must not contain repeated dimension names.")
    }
    along[which(sapply(along, is.character))] <- match(req_dimnames, dim_names)
    if (any(is.na(along))) {
      stop("Could not match all dimension names in 'indices' with dimension names in input array 'x'.")
    }
    along <- as.numeric(along)
  }

  # Check indices
  if (!is.list(indices)) {
    indices <- list(indices)
  }

  # Check parameter drop
  dims_to_drop <- c()
  if (is.character(drop)) {
    if (drop == 'all') {
      drop <- TRUE
    } else if (any(drop %in% c('selected', 'non-selected', 'none'))) {
      if (drop == 'selected') {
        dims_to_drop <- along[which(sapply(indices, length) == 1)]
      } else if (drop == 'non-selected') {
        dims_to_drop <- dim(x) == 1
        dims_to_drop[along] <- FALSE
        dims_to_drop <- which(dims_to_drop)
      }
      drop <- FALSE
    } else {
      stop("Parameter 'drop' must be one of TRUE, FALSE, 'all', 'selected', 'non-selected', 'none'.")
    }
  } else if (!is.logical(drop)) {
    stop("Parameter 'drop' must be one of TRUE, FALSE, 'all', 'selected', 'non-selected', 'none'.")
  }

  # Take the subset
  nd <- length(dim(x))
  index <- as.list(rep(TRUE, nd))
  index[along] <- indices
  subset <- eval(as.call(c(as.name("["), as.name("x"), index, drop = drop)))
  # If dropped all dimensions, need to drop dimnames too
  if (is.character(dim_names) && drop == TRUE) {
    dim_names_to_remove <- unique(c(along[which(sapply(indices, length) == 1)],
                                    which(dim(x) == 1)))
    if (length(dim_names_to_remove) > 0) {
      dim_names <- dim_names[-dim_names_to_remove]
    }
  }

  # Amend the final dimensions and put dimnames and attributes
  metadata <- attributes(x)
  metadata[['dim']] <- dim(subset)
  if (length(dims_to_drop) > 0) {
    metadata[['dim']] <- metadata[['dim']][-dims_to_drop]
    if (is.character(dim_names)) {
      names(metadata[['dim']]) <- dim_names[-dims_to_drop]
      if ('dimensions' %in% names(attributes(x))) {
        metadata[['dimensions']] <- dim_names[-dims_to_drop]
      }
    }
  } else if (is.character(dim_names)) {
    names(metadata[['dim']]) <- dim_names
    if ('dimensions' %in% names(attributes(x))) {
      metadata[['dimensions']] <- dim_names
    }
  }
  attributes(subset) <- metadata
  subset
}
