#'Computes composites
#'
#'Composites a 3-d field var(x, y, time) according to the indices of 
#'mode/cluster occurrences in time and computes the pvalues (t-test). x and y 
#'are typically lon and lat, but function can accept other 2-d fields such as 
#'lat and depth, lon and depth, etc.
#'
#'@param var 3-dimensional array (x, y, time) containing the variable to 
#'  composite.
#'@param occ 1-dimensional array for the occurrence time series of 
#'  mode(s)/cluster(s).
#'  (*1) When one wants to composite all modes, e.g., all K = 3 clusters then 
#'    for example occurrences could look like: 1 1 2 3 2 3 1 3 3 2 3 2 2 3 2.
#'  (*2) Otherwise for compositing only the 2nd mode or cluster of the above 
#'    example occurrences should look like 0 0 1 0 1 0 0 0 0 1 0 1 1 0 1.
#'@param lag Lag time step (an integer), e.g., for lag = 2 composite will 
#'  use +2 occurrences (i.e., shifted 2 time steps forward). Default is lag = 0.
#'@param eno For using the effective sample size (TRUE) or the total sample 
#'  size (FALSE that is the default) for the number of degrees of freedom.
#'@param fileout Name of the .sav output file (NULL is the default).
#'
#'@return 
#'\item{$composite}{ 
#'  3-d array (x, y, k) containing the composites k=1,..,K for case (*1)
#'  or only k=1 for any specific cluster, i.e., case (*2).
#'}
#'\item{$pvalue}{
#'  3-d array (x, y, k) containing the pvalue of the 
#'  composites obtained through a t-test that accounts for the serial
#'  dependence of the data with the same structure as Composite.
#'}
#'@keywords datagen
#'@author History:
#'  0.1  #  2014-08  (N.S. Fuckar, \email{neven.fuckar@@bsc.es}) # Original code
#'
#'@examples
#'blank <- array(0, dim=c(20, 10, 30))
#'
#'x1 <- blank
#'t1 <- blank
#'f1 <- blank
#'
#'for (i in 1:20) {
#'  x1[i,,] <- i
#'}
#'
#'for (i in 1:30) {
#'  t1[,,i] <- i
#'}
#'
#'# This is 2D propagating sin wave example, where we use (x,y,t) structure of 
#'# f1 wave field. Compositing (like using stroboscopicc light) at different time 
#'# steps can lead to modification or cancelation of wave pattern.
#'
#'for (i in 1:20) {
#'  for (j in 1:30) {
#'    f1[i,,j] <- 3*sin(2*pi*x1[i,,j]/5. - 2*pi*t1[i,,j]/6.)
#'  }
#'}
#'
#'occ1 <- rep(0, 30)
#'occ1[c(2, 5, 8, 11, 14, 17, 20, 23)] <- 1
#'
#'filled.contour(Composite(var=f1, occ=occ1)$composite[,,1])
#'
#'occ2 <- rep(0, 30)
#'occ2[c(3, 9, 15, 21)] <- 1
#'
#'filled.contour(Composite(var=f1, occ=occ2)$composite[,,1])
#'@importFrom stats sd pt
#'@export
Composite <- function(var, occ, lag = 0, eno = FALSE, fileout = NULL) {

  if ( dim(var)[3] != length(occ) ) {
     stop("Temporal dimension of var is not equal to length of occ.")
  }
  K         <- max(occ)
  composite <- array(dim = c(dim(var)[1:2], K))
  tvalue    <- array(dim = dim(var)[1:2])
  dof       <- array(dim = dim(var)[1:2])
  pvalue    <- array(dim = c(dim(var)[1:2], K))

  if (eno == TRUE) { 
    n_tot <- Eno(var, posdim = 3)
  } else {
    n_tot <- length(occ)
  }
  mean_tot <- Mean1Dim(var, posdim = 3, narm = TRUE)
  stdv_tot <- apply(var, c(1, 2), sd, na.rm = TRUE) 

  for (k in 1 : K) {

    indices <- which(occ == k) + lag
    toberemoved <-  which(0 > indices | indices > dim(var)[3])

    if (length(toberemoved) > 0) {
        indices <- indices[-toberemoved]
    }
    if (eno == TRUE) {
        n_k <- Eno(var[, , indices], posdim = 3)
    } else {
        n_k <- length(indices)
    }
    if (length(indices) == 1) {
        composite[, , k] <- var[, , indices] 
        warning(paste("Composite", k, "has length 1 and pvalue is NA."))
    }  else {
        composite[, , k] <- Mean1Dim(var[, , indices], posdim = 3, narm = TRUE)
    }
    stdv_k <- apply(var[, , indices], c(1, 2), sd, na.rm = TRUE)
    
    tvalue <- (mean_tot - composite[, , k]) / 
               sqrt(stdv_tot ^ 2 / n_tot + stdv_k ^ 2 / n_k)
    dof <- (stdv_tot ^ 2 / n_tot + stdv_k ^ 2 / n_k) ^ 2 / 
           ((stdv_tot ^ 2 / n_tot) ^ 2 / (n_tot - 1) +
           (stdv_k ^ 2 / n_k) ^ 2 / (n_k - 1))
    pvalue[, , k] <- 2 * pt(-abs(tvalue), df = dof)
  }  

  if (is.null(fileout) == FALSE) { 
    output <- list(composite = composite, pvalue = pvalue)   
    save(output, file = paste(fileout, '.sav', sep = ''))
  }
  
  invisible(list(composite = composite, pvalue = pvalue)) 
}
