\name{sampleDepthData}
\docType{data}
\alias{sampleDepthData}
\title{Sample of Experimental Data for Forecast Verification In Function Of 
       Latitudes And Depths}
\description{
  This data set provides data in function of latitudes and depths for the 
  variable 'tos', i.e. sea surface temperature, from the decadal climate 
  prediction experiment run at IC3 in the context of the CMIP5 project. 
  Its name within IC3 local database is 'i00k'.
}
\usage{
  data(sampleDepthData)
}
\format{
  The data set provides with a variable named 'sampleDepthData'.

  sampleDepthData$exp is an array that contains the experimental data and the 
  dimension meanings and values are:
    c(# of experimental datasets, # of members, # of starting dates, 
      # of lead-times, # of depths, # of latitudes)
    c(1, 5, 3, 60, 7, 21)

  sampleDepthData$obs should be an array that contained the observational data 
  but in this sample is not defined (NULL).

  sampleDepthData$depths is an array with the 7 longitudes covered by the data.

  sampleDepthData$lat is an array with the 21 latitudes covered by the data.
}
\keyword{datasets}
