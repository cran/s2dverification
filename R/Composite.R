Composite <- function(var, occ, lag=0, eno=FALSE, fileout=NULL) {

  if ( dim(var)[3]!=length(occ) ) { stop("temporal dimension of var is not equal to length of occ") }

  K         <- max(occ)
  composite <- array(dim = c(dim(var)[1:2], K))
  tvalue    <- array(dim = dim(var)[1:2])
  dof       <- array(dim = dim(var)[1:2])
  pvalue    <- array(dim = c(dim(var)[1:2], K))

  if ( eno==TRUE ) { n_tot <- Eno(var, posdim=3) }
  else { n_tot <- length(occ) }

  mean_tot <- Mean1Dim(var, posdim=3, narm=TRUE)
  stdv_tot <- apply(var, c(1,2), sd, na.rm=TRUE) 

  for (k in 1:K) {
    indices <- which(occ==k)+lag

    toberemoved=which(0>indices|indices>dim(var)[3])
    if ( length(toberemoved) > 0 ) { indices=indices[-toberemoved] }

    if ( eno==TRUE ) { n_k <- Eno(var[,,indices], posdim=3) }
    else { n_k <- length(indices) }

    composite[,,k] <- Mean1Dim(var[,,indices], posdim=3, narm=TRUE)
    stdv_k         <- apply(var[,,indices], c(1,2), sd, na.rm=TRUE)
    
    tvalue[,]   <- (mean_tot - composite[,,k])/sqrt(stdv_tot^2/n_tot + stdv_k^2/n_k)
    dof[,]      <- (stdv_tot^2/n_tot + stdv_k^2/n_k)^2/((stdv_tot^2/n_tot)^2/(n_tot - 1) +  (stdv_k^2/n_k)^2/(n_k - 1))
    pvalue[,,k] <- 2*pt(-abs(tvalue[,]), df=dof[,])
  }  

  if ( is.null(fileout)==FALSE ) { 
    output <- list(composite=composite, pvalue=pvalue)   
    save(output,file=paste(fileout,'.sav',sep=''))
  }
  
  invisible(list(composite = composite, pvalue = pvalue)) 
}
