\name{PlotAno}
\alias{PlotAno}
\title{
Plot Raw Or Smoothed Anomalies
}
\description{
Plots timeseries of raw or smoothed anomalies of any variable output from \code{Load()} or \code{Ano()} or or \code{Ano_CrossValid()} or \code{Smoothing()}.
}
\usage{
PlotAno(exp_ano, obs_ano = NULL, sdates, toptitle = rep("", 15), 
        ytitle = rep("", 15), limits = NULL, legends = NULL, freq = 12, 
        biglab = FALSE, fill = TRUE, memb = TRUE, ensmean = TRUE, 
        linezero = FALSE, points = FALSE, vlines = NULL, sizetit = 1, 
        fileout = paste0("output", 1:5, "_plotano.eps"), 
        width = 8, height = 5, size_units = 'in', res = 100, ...)
}
\arguments{
  \item{exp_ano}{
Array containing the experimental data:\cr
  c(nmod/nexp, nmemb/nparam, nsdates, nltime).
  }
  \item{obs_ano}{
Optional matrix containing the observational data:\cr
  c(nobs, nmemb, nsdates, nltime)
  }
  \item{sdates}{
List of starting dates: c('YYYYMMDD','YYYYMMDD').
  }
  \item{toptitle}{
Main title for each experiment: c('',''), optional.
  }
  \item{ytitle}{
Title of Y-axis for each experiment: c('',''), optional.
  }
  \item{limits}{
c(lower limit, upper limit): limits of the Y-axis, optional.
  }
  \item{legends}{
List of observational dataset names, optional.
  }
  \item{freq}{
1 = yearly, 12 = monthly, 4 = seasonal, ... Default: 12.
  }
  \item{biglab}{
TRUE/FALSE for presentation/paper plot. Default = FALSE.
  }
  \item{fill}{
TRUE/FALSE if the spread between members should be filled. Default = TRUE.
  }
  \item{memb}{
TRUE/FALSE if all members/only the ensemble-mean should be plotted.\cr
Default = TRUE.
  }
  \item{ensmean}{
TRUE/FALSE if the ensemble-mean should be plotted. Default = TRUE.
  }
  \item{linezero}{
TRUE/FALSE if a line at y=0 should be added. Default = FALSE.
  }
  \item{points}{
TRUE/FALSE if points instead of lines should be shown. Default = FALSE.
  }
  \item{vlines}{
List of x location where to add vertical black lines, optional.
  }
  \item{sizetit}{
Multiplicative factor to scale title size, optional.
  }
  \item{fileout}{
Name of the output file for each experiment: c('',''). 
Extensions allowed: eps/ps, jpeg, png, pdf, bmp and tiff. If filenames with 
different extensions are passed, it will be considered only the first one 
and it will be extended to the rest. \cr
Default = c('output1_plotano.eps', 'output2_plotano.eps', 
            'output3_plotano.eps', 'output4_plotano.eps', 
            'output5_plotano.eps')
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
    
    adj ann ask bg bty cex.sub cin col.axis col.lab col.main col.sub cra crt csi cxy err family fg fig font font.axis font.lab font.main font.sub lend lheight ljoin lmitre mar mex mfcol mfrow mfg mkh oma omd omi page plt smo srt tck tcl usr xaxp xaxs xaxt xlog xpd yaxp yaxs yaxt ylbias ylog

    For more information about the parameters see `par`
  }
}
\examples{
# Load sample data as in Load() example:
example(Load)
clim <- Clim(sampleData$mod, sampleData$obs)
ano_exp <- Ano(sampleData$mod, clim$clim_exp)
ano_obs <- Ano(sampleData$obs, clim$clim_obs)
runmean_nb_months <- 12
dim_to_smooth <- 4  # Smooth along lead-times
smooth_ano_exp <- Smoothing(ano_exp, runmean_nb_months, dim_to_smooth)
smooth_ano_obs <- Smoothing(ano_obs, runmean_nb_months, dim_to_smooth)
PlotAno(smooth_ano_exp, smooth_ano_obs, startDates, 
        toptitle = paste('smoothed anomalies'), ytitle = c('K', 'K', 'K'), 
        legends = 'ERSST', biglab = FALSE, fileout = 'tos_ano.eps')
}
\author{
History:\cr
0.1  -  2011-03  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to CRAN
}
\keyword{dynamic}
