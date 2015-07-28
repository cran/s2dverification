ConfigShowSimilarVars <- function(configuration, var_name, n_results = 5) {
  ## Simon White: http://www.catalysoft.com/articles/StrikeAMatch.html
  getBigrams <- function(str) {
    bigramlst <- list()
    for (i in 1:(nchar(str) - 1)) {
      bigramlst[[i]] <- substr(str, i, i + 1)
    }
    return(bigramlst)
  }

  strSimilarity <- function(str1, str2) {
    str1 <- tolower(str1)
    str2 <- tolower(str2)

    if (is.null(str1)) {
      str1 <- ""
    } else if (is.na(str1)) {
      str1 <- ""
    }
    if (is.null(str2)) {
      str2 <- ""
    } else if (is.na(str2)) {
      str2 <- ""
    }
    if (nchar(str1) <= 1 && nchar(str2) <= 1) {
      return (ifelse(str1 == str2, 1, 0))
    } else if (nchar(str1) == 1) {
      return (ifelse(grepl(str1, str2, fixed = TRUE), 1, 0))
    } else if (nchar(str2) == 1) {
      return (ifelse(grepl(str2, str1, fixed = TRUE), 1, 0))
    } else if (nchar(str1) == 0 || nchar(str2) == 0) {
      return (0)
    } else {
      pairs1 <- getBigrams(str1)
      pairs2 <- getBigrams(str2)
      unionlen <- length(pairs1) + length(pairs2)
      hit_count <- 0
      for (x in 1:length(pairs1)) {
        for(y in 1:length(pairs2)) {
          if (pairs1[[x]] == pairs2[[y]]) {
            hit_count <- hit_count + 1
          }
        }
      }
      return ((2.0 * hit_count) / unionlen)
    }
  }

  similarities <- list()
  similarities$two_d_vars <- unlist(as.vector(sapply(configuration$two_d_vars, strSimilarity, var_name)))
  similarities$global_mean_vars <- unlist(as.vector(sapply(configuration$global_mean_var, strSimilarity, var_name)))
  n_results <- min(n_results, length(similarities$two_d_vars) + length(similarities$global_mean_vars))
  threshold <- sort(c(similarities$two_d_vars, similarities$global_mean_vars), decreasing = TRUE)[n_results]
  n_minimums <- sum(sort(c(similarities$two_d_vars, similarities$global_mean_vars), decreasing = TRUE)[1:n_results] == threshold)

  matches <- list()
  n_picked_minimums <- 0
  for (var_type in c('two_d_vars', 'global_mean_vars')) {
    matches[[var_type]] <- configuration[[var_type]][which(similarities[[var_type]] > threshold)]
    if (n_picked_minimums < n_minimums) {
      minimums <- configuration[[var_type]][which(similarities[[var_type]] == threshold)]
      if (length(minimums) + n_picked_minimums > n_minimums) {
        minimums <- minimums[1:(n_minimums - n_picked_minimums)]
      }
      matches[[var_type]] <- c(matches[[var_type]], minimums)
      n_picked_minimums <- n_picked_minimums + length(minimums)
    }
  }

  ConfigShowVars(matches)
  invisible(matches)
}
