GenSeries <- function(n, alpha, mean, std) {
  # This functions generates AR1 processes containing n data, with alpha as 
  # autocorrelation at lag 1, and mean and standard deviation provided by the 
  # mean and std arguments.
  #
  # Params:
  #   n: Length of data to be generated.
  #   alpha: Autocorrelation at lag 1.
  #   mean: Mean of the data.
  #   std: Standard deviation of the data.
  #
  # Returns:
  #   AR1 process.
  #
  # History:
  #   1.0  #  2012-04  (L. Auger, ludovic.auger@meteo.fr)  #  Original code

  res <- vector("numeric", n)
  x <- mean
  stdterm <- std * (sqrt(1 - alpha ^ 2) / (1 - alpha))
  for (i in 1:100) {
    x <- alpha * x + (1 - alpha) * rnorm(1, mean, stdterm)
  }
  for (i in 1:n) {
    x <- alpha * x + (1 - alpha) * rnorm(1, mean, stdterm)
    res[i] <- x
  }
  
  res
}
