\name{Ano_CrossValid}
\alias{Ano_CrossValid}
\title{Computes Anomalies in Cross-Validation Mode}
\description{
  This function computes anomalies from experimental and observational 
  matrices output from \code{load()} by subtracting the climatologies 
  computed in a cross-validation mode and with a per-pair method.
}
\usage{
  Ano_CrossValid(var_exp, var_obs, memb = TRUE)
}
\arguments{
  \item{var_exp}{
    Model data:
      c(nmod/nexp, nmemb/nparam, nsdates, nltime) up to
      c(nmod/nexp, nmemb/nparam, nsdates, nltime, nlevel, nlat, nlon)
  }
  \item{var_obs}{
    Observational data:
      c(nobs, nmemb, nsdates, nltime) up to
      c(nobs, nmemb, nsdates, nltime, nlevel, nlat, nlon)
  }
  \item{memb}{
    memb: TRUE/FALSE (1 climatology for each member/1 climatology averaging all 
    the members). Default = TRUE.
  }
}
\value{
  \item{$ano_exp}{Matrix with same dimensions as var_exp}
  \item{$ano_obs}{Matrix with same dimensions as var_obs}
}
\examples{
  startDates <- c('19901101', '19951101', '20001101')
  sampleData <- Load('tos', c('i00k'), c('ERSST'), startDates, nleadtime = 124, 
                     leadtimemin = 1, leadtimemax = 60, output = 'areave', 
                     latmin = 30, latmax = 45, lonmin = 0, lonmax = 40)
  anomalies <- Ano_CrossValid(sampleData$mod, sampleData$obs)
  PlotAno(anomalies$ano_exp, anomalies$ano_obs, startDates, 
          toptitle = paste('anomalies'), ytitle = c('K', 'K', 'K'), 
          legends = 'ERSST', biglab = FALSE, fileout = 'tos_ano_crossvalid.eps')
}
\author{
  History:
  0.1  -  2011-12  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to CRAN
}
\keyword{datagen}
