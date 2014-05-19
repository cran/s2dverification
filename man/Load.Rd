\name{Load}
\alias{Load}
\title{Loads Experimental And Observational Data}
\description{
  This function loads experimental data together with the observational 
  data at the corresponding dates from the specified experimental and 
  observational datasets. It organizes those data into two matrices with 
  dimensions:
   c(nexp/nmod, nmembers, nsdates, nleadtimes) for the experiments
   c(nobs/nreanalyses, nmembers, nsdates, nleadtimes) for the observations
  up to
   c(nexp/nmod, nmembers, nsdates, nleadtimes, ndepth, nlat, nlon) 
   c(nobs/nreanalyses, nmembers, nsdates, nleadtimes, ndepth, nlat, nlon )
  Only a specified variable and set of starting dates is loaded.
  Once the two matrices are filled by calling this function, other functions
  of the 2dverification package that requires this data structure can be 
  executed (e.g: \code{Clim()} to compute climatologies,  \code{Ano()} to 
  compute anomalies, ...).

  Load() is currently repository-specific, that meaning you should overwrite it 
  so as to make it work properly with the directory tree and naming conventions 
  in your data repository.

  The function Load() must return 4 a list which contains four components: 
  $mod, $obs, $lat and $lon.
  $mod is the matrix containing the experimental data.
  $obs is the matrix containing the observational data.
  $lat and $lon are the latitudes and longitudes of the model grid (0 if the
  loaded variable is a global mean).

  The experimental data matrix will contain the values of a given variable for 
  all the experiments the user wish to compare for a set of starting dates, 
  lead-times, longitudes and latitudes interpolated to the common grid and with 
  the masks applied.
  The observational data matrix will contain the observed values of the same
  variable gathered from the specified observational datasets scanning the data 
  repository. The observational data is chosen so as to date-correspond the 
  experimental data.

  The data matrices must be built with the dimensions in the following order 
  and with the following lengths:
  1- The number of experimental datasets determined by the user (we need to 
  store data of all the experiments in the same array) for the experimental 
  matrix or the number of observational datasets available for validation in 
  the observational matrix.
  2- The greatest number of members across all experiments.
  3- The number of starting dates determined by the user (we need to store data 
  of each prediction of the model in each starting date).
  4- The greatest number of lead-times.
  5- The number of latitudes of the zone we want to consider.
  6- The number of longitudes of the zone we want to consider.

  Dimensions 5 and 6 are optional and their presence will depend on whether the 
  variable loaded is 2-dimensional or an area average. 
  In the case of an area average the dimensions of the matrix will be only the 
  first 4.

  For example, at IC3, we have implemented Load() to select the experimental 
  and observational datasets from which to load data, the variable to load, the 
  range of starting dates, the range of lead-times, the frequency and step in 
  which the data is loaded from the datasets, a common grid, a method of 
  interpolation, a mask to disable certain experimental values, a mask to 
  disable certain observational values and an output format:
    1) Time series of area-averaged variables over the specified domain. 
    2) Time series of meridional averages as a function of longitudes. 
    3) Time series of zonal averages as a function of latitudes.
    4) Time series of 2d fields.

  Furthermore, if the user loads a 2-dimensional variable, he/she can ask 
  Load() to load it as area averages, as longitudinal averages in function of 
  latitudes or as in latitudinal averages in function of longitudes. In these 
  cases the 5th and/or 6th dimension will also disappear.

  The observational data matrix is built very similarly to the experimental 
  matrix. The only differences are that the first dimension turns to be the
  number of observational datasets and the length of the 2nd dimension turns to 
  be the number of members of the observational datasets. Later, for each 
  lead-time we will store observational data into it if possible.

  Any value that can't be filled (either in the experimental matrix or the 
  observational matrix) should figure as an NA value.

}
\usage{
  Load(var, exp, obs = NULL, sdates, lonmin = 0, lonmax = 360, latmin = -90, 
  latmax = 90, nleadtime = NULL, nmember = NULL, leadtimemin = 1, 
  leadtimemax = NULL, storefreq = "monthly", sampleperiod = 1, 
  output = "areave", method = "conservative", grid = NULL, maskmod = list(NULL, 
  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
  NULL), maskobs = list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
  NULL, NULL, NULL, NULL, NULL, NULL))
}
\arguments{
  \item{var}{
    One of the following: 'tas', 'prlr', 'tos', 'g500', 'g200', 'ta50', 'psl', 
	'hflsd', 'hfssd', 'rls', 'rss', 'rsds', 'uas', 'vas', 'siaN', 'sieN', 
	'sivN', 'siaS', 'sieS', 'sivS', 'moc_40N55N_1-2km', 'moc_30N40N_1-2km', 
	'max_moc_38N50N_500m-2km', 'max_moc_40N', 'heatc', '0-315_heatc', 
	'373-657_heatc', '800-5350_heatc', 'mxl_heatc', 'NAtl_10N65N_heatc', 
	'NAtl_10N65N_0-315_heatc', 'NAtl_10N65N_373-657_heatc', 
	'NAtl_10N65N_800-5350_heatc', 'mxl_NAtl_10N65N_heatc', 'TAtl_30S30N_heatc', 
    'TAtl_30S30N_0-315_heatc', 'TAtl_30S30N_373-657_heatc',
    'TAtl_30S30N_800-5350_heatc', 'mxl_TAtl_30S30N_heatc', 'NPac_10N70N_heatc', 
	'NPac_10N70N_0-315_heatc', 'NPac_10N70N_373-657_heatc', 
	'NPac_10N70N_800-5350_heatc', 'mxl_NPac_10N70N_heatc', 'TPac_30S30N_heatc',
    'TPac_30S30N_0-315_heatc', 'TPac_30S30N_373-657_heatc', 
	'TPac_30S30N_800-5350_heatc', 'mxl_TPac_30S30N_heatc', 'Arc_65N90N_heatc', 
	'Arc_65N90N_0-315_heatc', 'Arc_65N90N_373-657_heatc', 
	'Arc_65N90N_800-5350_heatc', 'mxl_Arc_65N90N_heatc', 'Ant_90S60S_heatc',
    'Ant_90S60S_0-315_heatc', 'Ant_90S60S_373-657_heatc', 
	'Ant_90S60S_800-5350_heatc', 'mxl_Ant_90S60S_heatc', 'TInd_30S30N_heatc', 
	'TInd_30S30N_0-315_heatc', 'TInd_30S30N_373-657_heatc', 
	'TInd_30S30N_800-5350_heatc', 'mxl_TInd_30S30N_heatc'.
  }
  \item{exp}{
    IMPORTANT: Place first the experiment with the largest number of members 
	and, if possible, with the largest number of leadtimes. If not possible, it 
	is mandatory to specify the argument nmember and nleadtime.
    c('EnsEcmwfDec', 'EnsUkmoDec', 'EnsCerfacsDec', 'EnsIfmDec', 'EnsEcmwfSeas', 
	  'EnsCmccSeas', 'EnsIfmSeas', 'EnsMetfrSeas', 'EnsUkmoSeas', 
	  'DePreSysAsimDec', 'DePreSysNoAsimDec', 'DePreSysAsimSeas', 
	  'ECMWF_S3Seas', 'ECMWF_S4_sea', 'ECMWF_S4_ann', 'hadcm3dec', 'miroc4dec', 
	  'miroc5dec', 'mri-cgcm3dec', 'cancm4dec1', 'cancm4dec2', 'cnrm-cm5dec',
      'knmidec', 'smhidec', 'mpimdec', 'gfdldec', 'cmcc-cmdec', 'ipsldec', 
	  'bccdec', 'gfdlhis', 'ipslhis', 'cmcc-cmhis', 'bcchis', 'i00k', 'b013', 
	  'b014', 'yve2', \dots)
  }
  \item{obs}{
    c('ERA40', 'NCEP', 'ERAint', '20thCv2', 'GHCN', 'GHCNERSSTGISS', 'ERSST', 'HadISST', 
	  'GPCP', 'GPCC', 'CRU', 'HadSLP', 'NSIDC', 'PIOMAS', 'UCL', 'DS94', 
	  'OAFlux', 'DFS4.3', 'NCDCglo', 'NCDCland', 'NCDCoc', 'GISSglo', 
	  'GISSland', 'GISSoc', 'HadCRUT3glo', 'HadCRUT4', 
	  'HadSST2oc', 'CRUTEM3land')
  }
  \item{sdates}{
    c('YYYYMMDD', 'YYYYMMDD', ... , 'YYYYMMDD') 
  }
  \item{lonmin}{
    >= 0, default: 0
  }
  \item{lonmax}{
    <= 360, default: 360. If lonmin>lonmax, data accross Greenwich will be
    loaded.
  }
  \item{latmin}{
    >= -90, default: -90
  }
  \item{latmax}{
    <= 90, default: 90
  }
  \item{nleadtime}{
    Optional argument needed only if the first experiment in the list exp 
	does not have the largest number of leadtimes.
    Default: number of leadtimes of the first experiment.
  }
  \item{nmember}{
    Some experiments have more members in starting dates other than the first. 
    If it is the case in the first experiment specified, fill nmember with the 
    largest number of members.
  }
  \item{leadtimemin}{
    Load only the leadtimes from leadtimemin. 
    Default: 1.
  }
  \item{leadtimemax}{
    Load only the leadtimes before leadtimemax.
    Default: nleadtime.
  }
  \item{storefreq}{
    Frequency at which the data to be loaded are stored in the repository. Can 
	take values 'monthly' or 'daily'. 
    Default: 'monthly'.
  }
  \item{sampleperiod}{
    To load only a subset between leadtimemin and leadtimemax with a specified 
	period of subsampling. 
    Default: 1 (all leadtimes are loaded).
  }
  \item{output}{
    'areave' / 'lon' / 'lat' / 'lonlat'.
    1) Time series of area-averaged variables over the specified domain.
    2) Time series of meridional averages as a function of longitudes.
    3) Time series of zonal averages as a function of latitudes.
    4) Time series of 2d fields.
    Default: 'areave'
  }
  \item{method}{
    'bilinear' / 'bicubic' / 'conservative' / 'distance-weighted'
    Method of interpolation for 'lon' / 'lat' / 'lonlat' output options.
    Default: 'conservative'.
  }
  \item{grid}{
    To choose the output grid.
    Possible options: rNXxNY or tTRgrid, ex: r96x72, t106grid.
    Default: model grid.
    This argument needs to be filled if various experiments are loaded with 
    various grids.

    How grids work:

    We may load one or more experiments and corresponding data from one 
    or more observations.
    The most usual is that they are on different grids.

    We can set Load() to output the data in 4 different modes:

    If output = 'areave' all the data is area-averaged.
    If output = 'lon', the data are averaged along the latitude 
    dimension and output as a function of the longitudes.
    If output = 'lat', the data are averaged along the longitude 
    dimension and output as a function of the latitudes.
    If output = 'lonlat', the data are output as a function of the 
    longitudes and latitudes.

    We refer to the last 3 modes as maps.

    If we want to get data output in maps but the experiments and
    observations are on different grids, we need to specify a common 
    gridding (with the  parameter grid) onto which the loaded data 
    will be interpolated.
    If we don't do, Load() will crash.

    If we want to get area-averages, the experiments and observations 
    can be on different grids and there is no need to specify a common 
    grid: Load() will calculate their mean over the original grids 
    unless a common grid is specified.

    Whenever we specify a common gridding, all data (observational and
    experimental) is interpolated onto it.
  }
  \item{maskmod}{
    list(mask[lon,lat]) = 1/0: kept/removed grid cell over the entire model 
	domains.
    Warning: list() compulsory even if 1 model!!!
    Default: 1 everywhere.

    How masks work:

    When we want to disable certain values we can specify masks. The 
    disabled values will be replaced with NA values.

    In the case that you wish to load data as area-averages:
    You can specify a mask for each experiment and observation.
    Each area average will be performed on the original experiment 
    or observation grid and after applying the respective mask.

    In the case that you wish to load maps:
    You will have defined a mandatory common grid which all the data 
    will be interpolated onto.
    If you want to specify a mask, you will have to provide it 
    already interpolated onto the common grid (you may use 'cdo' for 
    this purpose).
    It is not usual to apply different masks on experimental 
    datasets on the same grid. All the experiment masks are expected 
    to be the same.

    Defining masks for the observational data will be useless: the 
    same mask applied to the experimental data will be applied to
    the observational data because the observations are interpolated 
    on the common grid for which the mask has already been defined 
    when loading experimental data.
  }
  \item{maskobs}{
    list(mask[lon,lat]) = 1/0: kept/removed grid cell over the entire observed 
	domains, only necessary for 'areave' output option.
    Warning: list() compulsory even if 1 dataset !!!
    Default: 1 everywhere.
  }
}
\details{
  The two output matrices have between 2 and 6 dimensions:
    1) Number of experimental/observational datasets.
    2) Number of members.
    3) Number of startdates.
    4) Number of leadtimes.
    5) Number of latitudes (optional).
    6) Number of longitudes (optional).
  but the two matrices have the same number of dimensions and only the first
  two dimensions can have different lengths depending on the input arguments.
    
  For a detailed explanation of the process, read the documentation attached to 
  the package or check the comments in the code.
}
\value{
  \item{$mod}{
    Model outputs.
    If output = 'areave', matrix with dimensions
      c(nmod/nexp, nmemb/nparam, nsdates, nltime)
    If output = 'lat', matrix with dimensions
      c(nmod/nexp, nmemb/nparam, nsdates, nltime, nlat)
    If output = 'lon', matrix with dimensions
      c(nmod/nexp, nmemb/nparam, nsdates, nltime, nlon)
    If output = 'lonlat', matrix with dimensions
      c(nmod/nexp, nmemb/nparam, nsdates, nltime, nlat, nlon)
  }
  \item{$obs}{
    Observations. Matrix with same dimensions as '$mod' except along the first
      two.
    If output = 'areave', matrix with dimensions
      c(nobs, nmemb, nsdates, nltime)
    If output = 'lat', matrix with dimensions
      c(nobs, nmemb, nsdates, nltime, nlat)
    If output = 'lon', matrix with dimensions
      c(nobs, nmemb nsdates, nltime, nlon)
    If output = 'lonlat', matrix with dimensions
      c(nobs, nmemb, nsdates, nltime, nlat, nlon)

  }
  \item{$lat}{
    Latitudes of the output grid (default: model grid of the first experiment).
  }
  \item{$lon}{
    Longitudes of the output grid (default: model grid of the first experiment).
  }
}
\author{
  History:
  0.1  -  2011-03  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to CRAN
}
\examples{
  startDates <- c('19901101', '19951101', '20001101')
  sampleData <- Load('tos', c('i00k'), c('ERSST'), startDates, nleadtime = 124, 
                     leadtimemin = 1, leadtimemax = 60, output = 'areave', 
                     latmin = 30, latmax = 45, lonmin = 0, lonmax = 40)
  PlotAno(sampleData$mod, sampleData$obs, startDates, 
          toptitle = "Mediterannean SST", ytitle = 'K', 
          fileout = 'tos_load_data.eps')
}
\keyword{datagen}
