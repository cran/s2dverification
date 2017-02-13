\name{RatioSDRMS}
\alias{RatioSDRMS}
\alias{.RatioSDRMS}
\title{
Computes the ratio between the ensemble spread and RMSE 
}
\description{
Arrays var_exp & var_obs should have dimensions between\cr
  c(nmod/nexp, nmemb/nparam, nsdates, nltime)\cr
and\cr
  c(nmod/nexp, nmemb/nparam, nsdates, nltime, nlevel, nlat, nlon)\cr
The ratio between the standard deviation of the members around the ensemble mean in var_exp and the RMSE between var_exp and var_obs is output for each experiment and each observational dataset.\cr
The p-value is provided by a one-sided Fischer test.\cr
\cr
.RatioSDRMS provides the same functionality but taking a matrix of ensemble members as input (exp).
}
\usage{
RatioSDRMS(var_exp, var_obs, pval = TRUE)

.RatioSDRMS(exp, obs, pval = TRUE)
}
\arguments{
 \item{var_exp}{
Model data:\cr
  c(nmod/nexp, nmemb/nparam, nsdates, nltime) up to\cr
  c(nmod/nexp, nmemb/nparam, nsdates, nltime, nlevel, nlat, nlon)
  }
  \item{var_obs}{
Observational data:\cr
  c(nobs, nmemb, nsdates, nltime) up to\cr
  c(nobs, nmemb, nsdates, nltime, nlevel, nlat, nlon)
  }
  \item{pval}{
Whether to compute the p-value of Ho : SD/RMSE = 1 or not.
  }
  \item{exp}{
N by M matrix of N forecasts from M ensemble members.
  }
  \item{obs}{
Vector of the corresponding observations of length N.
  }
}
\value{
RatioSDRMS: Array with dimensions c(nexp/nmod, nobs, 1 or 2, nltime) up to
  c(nexp/nmod, nobs, 1 or 2, nltime, nlevel, nlat, nlon). \cr
The 3rd dimension corresponds to the ratio (SD/RMSE) and the p.value (only present if \code{pval = TRUE}) of the one-sided Fisher test with Ho: SD/RMSE = 1.\cr
\cr
.RatioSDRMS:
  \itemize{
    \item{$ratio}{
The ratio of the ensemble spread and RMSE, 
    }
    \item{$p_val}{
Corresponds to the p values of the ratio (only present if \code{pval = TRUE}).
    }
  }
}
\examples{
# Load sample data as in Load() example:
example(Load)
rsdrms <- RatioSDRMS(sampleData$mod, sampleData$obs)
# Reorder the data in order to plot it with PlotVsLTime
rsdrms_plot <- array(dim = c(dim(rsdrms)[1:2], 4, dim(rsdrms)[4]))
rsdrms_plot[, , 2, ] <- rsdrms[, , 1, ]
rsdrms_plot[, , 4, ] <- rsdrms[, , 2, ]
PlotVsLTime(rsdrms_plot, toptitle = "Ratio ensemble spread / RMSE", ytitle = "", 
            monini = 11, limits = c(-1, 1.3), listexp = c('CMIP5 IC3'), 
            listobs = c('ERSST'), biglab = FALSE, siglev = TRUE, 
            fileout = 'tos_rsdrms.eps')

# The following example uses veriApply combined with .RatioSDRMS instead of RatioSDRMS
  \dontrun{
require(easyVerification)  
RatioSDRMS2 <- s2dverification:::.RatioSDRMS
rsdrms2 <- veriApply("RatioSDRMS2",
                     sampleData$mod,
                     # see ?veriApply for how to use the 'parallel' option
                     Mean1Dim(sampleData$obs, 2),
                     tdim = 3, ensdim = 2)
  }
}
\author{
History:\cr
0.1  -  2011-12  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau-manubens at ic3.cat})  -  Formatting to CRAN\cr
1.1  -  2017-02  (A. Hunter, \email{alasdair.hunter at bsc.es})  -  Adapted to veriApply()
}
\keyword{datagen}
