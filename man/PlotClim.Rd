\name{PlotClim}
\alias{PlotClim}
\title{
Plots Climatologies
}
\description{
Plots climatologies as a function of the forecast time for any index output from \code{Clim()} and organized in matrix with dimensions:\cr
  c(nmod/nexp, nmemb/nparam, nltime) or c(nmod/nexp, nltime) for the experiment data\cr
  c(nobs, nmemb, nltime) or c(nobs, nltime) for the observational data
}
\usage{
PlotClim(exp_clim, obs_clim = NULL, toptitle = "", ytitle = "", monini = 1, 
         freq = 12, limits = NULL, listexp = c("exp1", "exp2", "exp3"), 
         listobs = c("obs1", "obs2", "obs3"), biglab = FALSE, leg = TRUE, 
         fileout = "output_plotclim.eps", sizetit = 1)
}
\arguments{
  \item{exp_clim}{
Matrix containing the experimental data with dimensions:\cr
  c(nmod/nexp, nmemb/nparam, nltime) or c(nmod/nexp, nltime) 
  }
  \item{obs_clim}{
Matrix containing the observational data (optional) with dimensions:\cr
  c(nobs, nmemb, nltime) or c(nobs, nltime) 
  }
  \item{toptitle}{
Main title, optional
  }
  \item{ytitle}{
Title of Y-axis, optional.
  }
  \item{monini}{
Starting month between 1 and 12. Default = 1.
  }
  \item{freq}{
1 = yearly, 12 = monthly, 4 = seasonal, ... Default = 12.
  }
  \item{limits}{
c(lower limit, upper limit): limits of the Y-axis, optional.
  }
  \item{listexp}{
List of experiment names, optional.
  }
  \item{listobs}{
List of observational dataset names, optional.
  }
  \item{biglab}{
TRUE/FALSE for presentation/paper plot. Default = FALSE.
  }
  \item{leg}{
TRUE/FALSE to plot the legend or not.
  }
  \item{fileout}{
Name of the output eps file.
  }
  \item{sizetit}{
Multiplicative factor to scale title size, optional.
  }
}
\examples{
# Load sample data as in Load() example:
example(Load)
clim <- Clim(sampleData$mod, sampleData$obs)
PlotClim(clim$clim_exp, clim$clim_obs, toptitle = paste('climatologies'), 
         ytitle = 'K', monini = 11, listexp = c('CMIP5 IC3'), 
         listobs = c('ERSST'), biglab = FALSE, fileout = 'tos_clim.eps')
}
\author{
History:\cr
0.1  -  2011-03  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to CRAN
}
\keyword{datagen}
