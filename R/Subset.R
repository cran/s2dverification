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
      metadata[['dimensions']] <- dim_names[-dims_to_drop]
    }
  } else if (is.character(dim_names)) {
    names(metadata[['dim']]) <- dim_names
    metadata[['dimensions']] <- dim_names
  }
  attributes(subset) <- metadata
  subset
}
