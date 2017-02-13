\name{PlotVsLTime}
\alias{PlotVsLTime}
\title{
Plots A Score Along The Forecast Time With Its Confidence Interval
}
\description{
Plots The Correlation (\code{Corr()}) or the Root Mean Square Error (\code{RMS()}) between the forecasted values and their observational counterpart or the slopes of their trends (\code{Trend()}) or the InterQuartile Range, Maximum-Mininum, Standard Deviation or Median Absolute Deviation of the Ensemble Members (\code{Spread()}), or the ratio between the Ensemble Spread and the RMSE of the Ensemble Mean (\code{RatioSDRMS()}) along the forecast time for all the input experiments on the same figure with their confidence intervals.
}
\usage{
PlotVsLTime(var, toptitle = "", ytitle = "", monini = 1, freq = 12, 
            nticks = NULL, limits = NULL, listexp = c("exp1", "exp2", "exp3"), 
            listobs = c("obs1", "obs2", "obs3"), biglab = FALSE, hlines = NULL,
            leg = TRUE, siglev = FALSE, sizetit = 1, show_conf = TRUE,
            fileout = "output_plotvsltime.eps", 
            width = 8, height = 5, size_units = 'in', res = 100, ...)
}
\arguments{
  \item{var}{
Matrix containing any Prediction Score with dimensions:\cr
     (nexp/nmod, 3/4 ,nltime)\cr
  or (nexp/nmod, nobs, 3/4 ,nltime)
  }
  \item{toptitle}{
Main title, optional.
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
  \item{nticks}{
Number of ticks and labels on the x-axis, optional.
  }
  \item{limits}{
c(lower limit, upper limit): limits of the Y-axis, optional.
  }
  \item{listexp}{
List of experiment names, optional.
  }
  \item{listobs}{
List of observation names, optional.
  }
  \item{biglab}{
TRUE/FALSE for presentation/paper plot. Default = FALSE.
  }
  \item{hlines}{
c(a,b, ..) Add horizontal black lines at Y-positions a,b, ...\cr
Default = NULL.
  }
  \item{leg}{
TRUE/FALSE if legend should be added or not to the plot. Default = TRUE.
  }
  \item{siglev}{
TRUE/FALSE if significance level should replace confidence interval.\cr
Default = FALSE.
  }
  \item{sizetit}{
Multiplicative factor to change title size, optional.
  }
  \item{show_conf}{
TRUE/FALSE to show/not confidence intervals for input variables.
  }
  \item{fileout}{
Name of output file. Extensions allowed: eps/ps, jpeg, png, pdf, bmp 
and tiff.\cr
Default = 'output_plotvsltime.eps'
  }
  \item{width}{
File width, in the units specified in the parameter size_units (inches by default). Takes 8 by default.
  }
  \item{height}{
File height, in the units specified in the parameter size_units (inches by default). Takes 5 by default.
  }
  \item{size_units}{
Units of the size of the device (file or window) to plot in. Inches ('in') by default. See ?Devices and the creator function of the corresponding device.
  }
  \item{res}{ 
Resolution of the device (file or window) to plot in. See ?Devices and the creator function of the corresponding device.
  }
  \item{...}{
    Arguments to be passed to the method. Only accepts the following
    graphical parameters:

    adj ann ask bg bty cex.sub cin col.axis col.lab col.main col.sub cra crt csi cxy err family fg fig font font.axis font.lab font.main font.sub lheight ljoin lmitre mar mex mfcol mfrow mfg mkh oma omd omi page pch plt smo srt tck tcl usr xaxp xaxs xaxt xlog xpd yaxp yaxs yaxt ylbias ylog
    
    For more information about the parameters see `par`
  }
}
\details{
Examples of input:\cr
  Model and observed output from \code{Load()} then \code{Clim()} then \code{Ano()} then \code{Smoothing()}:\cr
    (nmod, nmemb, nsdate, nltime) and (nobs, nmemb, nsdate, nltime)\cr
  then averaged over the members\cr
    \code{Mean1Dim(var_exp/var_obs, posdim = 2)}:\cr
      (nmod, nsdate, nltime) and (nobs, nsdate, nltime)\cr
  then passed through\cr
    \code{Corr(exp, obs, posloop = 1, poscor = 2)} or\cr
    \code{RMS(exp, obs, posloop = 1, posRMS = 2)}:\cr
      (nmod, nobs, 3, nltime)\cr
  would plot the correlations or RMS between each exp & each obs as a function of the forecast time.
}
\examples{
# Load sample data as in Load() example:
example(Load)
clim <- Clim(sampleData$mod, sampleData$obs)
ano_exp <- Ano(sampleData$mod, clim$clim_exp)
ano_obs <- Ano(sampleData$obs, clim$clim_obs)
runmean_months <- 12
dim_to_smooth <- 4  # Smooth along lead-times
smooth_ano_exp <- Smoothing(ano_exp, runmean_months, dim_to_smooth)
smooth_ano_obs <- Smoothing(ano_obs, runmean_months, dim_to_smooth)
dim_to_mean <- 2  # Mean along members
required_complete_row <- 3  # Discard startdates for which there are NA leadtimes
leadtimes_per_startdate <- 60
corr <- Corr(Mean1Dim(smooth_ano_exp, dim_to_mean), 
             Mean1Dim(smooth_ano_obs, dim_to_mean), 
             compROW = required_complete_row, 
             limits = c(ceiling((runmean_months + 1) / 2), 
                        leadtimes_per_startdate - floor(runmean_months / 2)))
PlotVsLTime(corr, toptitle = "correlations", ytitle = "correlation", 
            monini = 11, limits = c(-1, 2), listexp = c('CMIP5 IC3'), 
            listobs = c('ERSST'), biglab = FALSE, hlines = c(-1, 0, 1), 
            fileout = 'tos_cor.eps')
}
\author{
History:\cr
0.1  -  2011-03  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
0.2  -  2013-03  (I. Andreu-Burillo, \email{isabel.andreu-burillo at ic3.cat})  -  Introduced parameter sizetit\cr
0.3  -  2013-10  (I. Andreu-Burillo, \email{isabel.andreu-burillo at ic3.cat})  -  Introduced parameter show_conf\cr
1.0  -  2013-11  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to CRAN
}
\keyword{dynamic}
