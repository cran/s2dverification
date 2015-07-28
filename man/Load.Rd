\name{Load}
\alias{Load}
\title{Loads Experimental And Observational Data}
\description{
This function loads monthly or daily data from a set of specified experimental datasets together with data that date-corresponds from a set of specified observational datasets. See parameters 'storefreq', 'sampleperiod', 'exp' and 'obs'.\cr\cr
Load() arranges the data into two arrays with a similar format with the following dimensions:
  \enumerate{
    \item{The number of experimental datasets determined by the user through the argument 'exp' for the experimental data array or the number of observational datasets available for validation for the observational array determined as well by the user through the argument 'obs'.}
    \item{The greatest number of members across all experiments in the experimental data array or across all observational datasets in the observational data array.}
    \item{The number of starting dates determined by the user through the 'sdates' argument (we need to store data of each prediction of the model in each starting date).}
    \item{The greatest number of lead-times.}
    \item{The number of latitudes of the zone we want to consider.}
    \item{The number of longitudes of the zone we want to consider.}
  }
Dimensions 5 and 6 are optional and their presence depends on the type of the specified variable (global mean or 2-dimensional) and on the selected output type (area averaged time series, latitude averaged time series, longitude averaged time series or 2-dimensional time series).\cr
In the case of loading an area average the dimensions of the arrays will be only the first 4.\cr\cr

Only a specified variable and a set of starting dates is loaded from each experiment. See parameters 'var' and 'sdates'.\cr
Afterwards, observational data that matches every starting date and lead-time of every experimental dataset is fetched in the file system (so, if two predictions at two different start dates overlap, some observational values will be loaded more than once).\cr
If no data is found in the file system for an experimental or observational array point it is filled with an NA value.\cr\cr

If the specified output is 2-dimensional or latitude- or longitude-averaged time series all the data is interpolated into a common grid. If the specified output type is area averaged time series the data is averaged on the individual grid of each dataset but can also be averaged after interpolating into a common grid. See parameters 'grid' and 'method'.\cr
Once the two arrays are filled by calling this function, other functions in the s2dverification package that receive as inputs data formatted in this data structure can be executed (e.g: \code{Clim()} to compute climatologies, \code{Ano()} to compute anomalies, ...).\cr\cr

Load() has many additional parameters to disable values and trim dimensions of selected variable, even masks can be applied to 2-dimensional variables. See parameters 'nmember', 'nmemberobs', 'nleadtime', 'leadtimemin', 'leadtimemax', 'sampleperiod', 'lonmin', 'lonmax', 'latmin', 'latmax', 'maskmod', 'maskobs', 'varmin', 'varmax'.\cr\cr

The parameters 'exp' and 'obs' are lists of dataset identifiers. To fetch the data in the repository associated to each dataset, Load() reads a configuration file which associates each pair of (dataset ID, variable name) to a corresponding path pattern in which the dataset and variable data is stored. These patterns can contain wildcards and variables tags that will be replaced automatically by Load() with the specified starting dates, member numbers, variable name, etc. Furthermore, each pair can be associated not only to a path but also to other values such as grid type, maximum or minimum allowed values, etc., that make Load() capable of loading data stored in multiple formats and naming conventions. See '?ConfigFileOpen' and parameters 'configfile', 'suffixexp' and 'suffixobs' in this help page for more information.\cr
Only NetCDF files are supported. OPeNDAP URLs to NetCDF files are also supported.\cr\cr

All in all, Load() can load 2-dimensional or global mean variables in any of the following formats:
  \itemize{
    \item{experiments:
      \itemize{
        \item{file per ensemble per starting date (YYYY, MM and DD somewhere in the path)}
        \item{file per member per starting date (YYYY, MM, DD and MemberNumber somewhere in the path, exps with different numbers of members supported)}
      }
(YYYY, MM and DD specify the starting dates of the predictions)
    }
    \item{observations:
      \itemize{
        \item{file per ensemble per month (YYYY and MM somewhere in the path)}
        \item{file per member per month (YYYY, MM and MemberNumber somewhere in the path, obs with different numbers of members supported)}
        \item{file per dataset (No constraints in the path but the time axes in the file have to be properly defined)}
      }
(YYYY and MM correspond to the actual month data in the file)
    }
  }
In all the formats the data can be stored in a daily or monthly frequency, or a multiple of these (see parameters 'storefreq' and 'sampleperiod').\cr
All the data files must contain the target variable with the same name and, that variable, must be defined over time, members, latitude and longitude dimensions in this order, time being the record dimension, members being only required if in 'file per member per starting date' format, and longitude and latitude being only required if is a two-dimensional variable.\cr
In the case of a two-dimensional variable, the variables longitude and latitude must be defined inside the data file too and must have the same names as the dimension for longitudes and latitudes respectively.\cr
The names of these dimensions (and longitude and latitude variables) can be configured in the configuration file (see ?ConfigFileOpen) via the variables DEFAULT_DIM_NAME_LONGITUDES, DEFAULT_DIM_NAME_LATITUDES and DEFAULT_DIM_NAME_MEMBERS. When generating a new configuration file (see ?ConfigFileCreate) these take as default 'longitude', 'latitude' and 'ensemble'.\cr
All the data files are expected to have numeric values representable with 32 bits. Be aware when choosing the fill values or infinite values in the datasets to load.\cr\cr
The Load() function returns a named list whith four components: 'mod', 'obs', 'lat' and 'lon'.
  \itemize{
    \item{'mod' is the array that contains the experimental data.}
    \item{'obs' is the array that contains the observational data.}
    \item{'lat' and 'lon' are the latitudes and longitudes of the grid into which the data is interpolated (0 if the loaded variable is a global mean or the output is an area average).}
  }
}
\usage{
Load(var, exp = NULL, obs = NULL, sdates, nmember = NULL, nmemberobs = NULL, 
     nleadtime = NULL, leadtimemin = 1, leadtimemax = NULL, 
     storefreq = 'monthly', sampleperiod = 1, lonmin = 0, lonmax = 360, 
     latmin = -90, latmax = 90, output = 'areave', method = 'conservative', 
     grid = NULL, maskmod = vector("list", 15), maskobs = vector("list", 15), 
     configfile = NULL, suffixexp = NULL, suffixobs = NULL, varmin = NULL, 
     varmax = NULL, silent = FALSE, nprocs = NULL, dimnames = NULL)
}
\arguments{
  \item{var}{
Name of the variable to load.\cr
Must have matching entries in the list of global mean or in the list of 2-dimensional variables in the configuration file.\cr
Also must have matching entries together with each dataset specified in 'exp' and 'obs' in any of the experimental dataset tables and in any of the observational dataset tables.\cr
This parameter is mandatory.\cr
Ex: 'tas'
  }
  \item{exp}{
Vector of experimental dataset IDs. Each must have matching entries together with the variable name specified in 'var' in any of the experimental datasets tables in the configuration file.\cr
IMPORTANT: Place first the experiment with the largest number of members and, if possible, with the largest number of leadtimes. If not possible, the arguments 'nmember' and/or 'nleadtime' should be filled to not miss any member or leadtime.\cr
If 'exp' is not specified or set to NULL, observational data is loaded for each start-date for a period of time as long as the time period between the first specified start date and the current date.\cr
Ex: c('ecmwf', 'i00k')
  }
  \item{obs}{
Vector of observational dataset IDs. Each must have matching entries together with the variable name specified in 'var' in any of the observational datasets tables in the configuration file.\cr
IMPORTANT: Place first the observation with the largest number of members. If not possible, the argument 'nmemberobs' should be filled to not miss any member.\cr
If 'obs' is not specified or set to NULL, no observational data is loaded.\cr
Ex: c('ERSST', 'ERAint')
  }
  \item{sdates}{
Vector of starting dates of the experimental runs to be loaded following the pattern 'YYYYMMDD'.\cr
This argument is mandatory.\cr
Ex: c('19601101', '19651101', '19701101')
  }
  \item{nmember}{
Vector with the numbers of members to load from the specified experimental datasets in 'exp'.\cr
If not specified, the number of members of the first experimental dataset is detected and replied to all the experimental datasets.\cr
If a single value is specified it is replied to all the experimental datasets.\cr
Data for each member is fetched in the file system. If not found is filled with NA values.\cr
An NA value in the 'nmember' list is interpreted as "fetch as many members as the greatest number of members among experimental datasets".\cr
Note: It is recommended to specify the number of members of the first experimental dataset if it is stored in file per member format because there are known issues in the automatic detection of members if the path to the dataset in the configuration file contains Shell Globbing wildcards such as '*'.\cr
Ex: c(4, 9)
  }
  \item{nmemberobs}{
Vector with the numbers of members to load from the specified observational datasets in 'obs'.\cr
If not specified, the number of members of the first observational dataset is detected and replied to all the observational datasets.\cr
If a single value is specified it is replied to all the observational datasets.\cr
Data for each member is fetched in the file system. If not found is filled with NA values.\cr
An NA value in the 'nmemberobs' list is interpreted as "fetch as many members as the greatest number of members among observational datasets".\cr
Note: It is recommended to specify the number of members of the first observational dataset if it is stored in file per member format because there are known issues in the automatic detection of members if the path to the dataset in the configuration file contains Shell Globbing wildcards such as '*'.\cr
Ex: c(1, 5)
  }
  \item{nleadtime}{
Largest number of lead-times among experimental datasets.\cr
Takes by default the number of lead-times of the first experimental dataset in 'exp'.\cr
If 'exp' is NULL this argument won't have any effect (see Load() description).
  }
  \item{leadtimemin}{
Only lead-times higher or equal to 'leadtimemin' are loaded. Takes by default value 1.
  }
  \item{leadtimemax}{
Only lead-times lower or equal to 'leadtimemax' are loaded. Takes by default the same value as 'nleadtime'.
  }
  \item{storefreq}{
Frequency at which the data to be loaded is stored in the file system. Can take values 'monthly' or 'daily'.\cr
By default it takes 'monthly'.\cr
Note: Data stored in other frequencies with a period which is divisible by a month can be loaded with a proper use of 'storefreq' and 'sampleperiod' parameters. It can also be loaded if the period is divisible by a day and the observational datasets are stored in a file per dataset format or 'obs' is empty.
  }
  \item{sampleperiod}{
To load only a subset between 'leadtimemin' and 'leadtimemax' with the period of subsampling 'sampleperiod'.\cr
Takes by default value 1 (all lead-times are loaded).\cr
See 'storefreq' for more information.
  }
  \item{lonmin}{
If a 2-dimensional variable is loaded, values at longitudes lower than 'lonmin' aren't loaded.\cr
Must take a value in the range [0, 360] (if negative longitudes are found in the data files these are translated to this range).\cr
It is set to 0 if not specified.\cr
If 'lonmin' > 'lonmax', data across Greenwich is loaded.
  }
  \item{lonmax}{
If a 2-dimensional variable is loaded, values at longitudes higher than 'lonmax' aren't loaded.\cr
Must take a value in the range [0, 360] (if negative longitudes are found in the data files these are translated to this range).\cr
It is set to 360 if not specified.\cr
If 'lonmin' > 'lonmax', data across Greenwich is loaded.
  }
  \item{latmin}{
If a 2-dimensional variable is loaded, values at latitudes lower than 'latmin' aren't loaded.\cr
Must take a value in the range [-90, 90].\cr
It is set to -90 if not specified.
  }
  \item{latmax}{
If a 2-dimensional variable is loaded, values at latitudes higher than 'latmax' aren't loaded.\cr
Must take a value in the range [-90, 90].\cr
It is set to 90 if not specified.
  }
  \item{output}{
This parameter determines the format in which the data is arranged in the output arrays.\cr
Can take values 'areave', 'lon', 'lat', 'lonlat'.\cr
    \itemize{
      \item{'areave': Time series of area-averaged variables over the specified domain.}
      \item{'lon': Time series of meridional averages as a function of longitudes.}
      \item{'lat': Time series of zonal averages as a function of latitudes.}
      \item{'lonlat': Time series of 2d fields.}
    }
Takes by default the value 'areave'. If the variable specified in 'var' is a global mean, this parameter is forced to 'areave'.\cr
All the loaded data is interpolated into the grid of the first experimental dataset except if 'areave' is selected. In that case the area averages are computed on each dataset original grid. A common grid different than the first experiment's can be specified through the parameter 'grid'. If 'grid' is specified when selecting 'areave' output type, all the loaded data is interpolated into the specified grid before calculating the area averages.
  }
  \item{method}{
This parameter determines the interpolation method to be used when regridding data (see 'output'). Can take values 'bilinear', 'bicubic', 'conservative', 'distance-weighted'.\cr
Takes by default the value 'conservative'.
  }
  \item{grid}{
A common grid can be specified through the parameter 'grid' when loading 2-dimensional data. Data is interpolated into this grid whichever 'output' type is specified. If the selected output type is 'areave' and a 'grid' is specified, the area averages are calculated after interpolating to the specified grid.\cr
If not specified and the selected output type is 'lon', 'lat' or 'lonlat', this parameter takes as default value the grid of the first experimental dataset, which is read from the configuration file.\cr
The grid must be supported by 'cdo' tools: rNXxNY or tTRgrid.\cr
Ex: 'r96x72' 
  }
  \item{maskmod}{
List of masks to be applied to the data of each experimental dataset respectively, if a 2-dimensional variable is specified in 'var'. Each mask is a matrix with dimensions c(longitudes, latitudes) with the same size as the common grid or with the same size of the grid of the corresponding experimental dataset if 'areave' output type is specified and no common 'grid' is specified. A value of 1 at a point of the mask keeps the original value at that point whereas a value of 0 disables it (replaces by an NA value).\cr
By default all values are kept (all ones).\cr
If you are loading maps (i.e. 'lonlat', 'lon' or 'lat' output types) all the data will be interpolated onto the common 'grid'. If you want to specify a mask, you will have to provide it already interpolated onto the common grid (you may use 'cdo' for this purpose). It is not usual to apply different masks on experimental datasets on the same grid, so all the experiment masks are expected to be the same. When loading maps, any masks defined for the observational data will be ignored to make sure the same mask is applied to the experimental and observational data.
Warning: list() compulsory even if loading 1 experimental dataset only!\cr
Ex: list(array(1, dim = c(num_lons, num_lats)))
  }
  \item{maskobs}{
List of masks to be applied to the data of each observational dataset respectively, if a 2-dimensional variable is specified in 'var'. Each mask is a matrix with dimensions c(longitudes, latitudes) with the same size as the common grid or with the same size of the associated observational dataset if 'areave' output type is specified and no common 'grid' is specified. A value of 1 at a point of the mask keeps the original value at that point whereas a value of 0 disables it (replaces by an NA value).\cr
By default all values are kept (all ones).\cr
If you are loading maps (i.e. 'lonlat', 'lon' or 'lat' output types) all the data will be interpolated onto the common 'grid'. If you want to specify a mask, you will have to provide it already interpolated onto the common grid (you may use 'cdo' for this purpose). It is not usual to apply different masks on experimental datasets on the same grid, so all the experiment masks are expected to be the same. When loading maps, any masks defined for the observational data will be ignored to make sure the same mask is applied to the experimental and observational data.
Warning: list() compulsory even if loading 1 observational dataset only!\cr
Ex: list(array(1, dim = c(num_lons, num_lats)))
  }
  \item{configfile}{
Path to the s2dverification configuration file from which information on location in file system (and other) of datasets is retreived.\cr
By default the configuration file used at IC3 is used (it is included in the package).\cr
Check the IC3's configuration file or a template of configuration file in the folder 'inst/config' in the package.\cr
Check further information on the configuration file mechanism in '?LoadConfigurationFile'.
  }
  \item{suffixexp}{
Vector of suffixes used to build the dataset path of the experimental datasets respectively. Some times paths to datasets can vary depending on slight differences in model adjustments, which leads to an addition of a suffix in the common path pattern for the rest of experimental datasets. In such situation 'suffixexp' can help.\cr
Each pair of experimental dataset ID and variable name can have a suffix associated in the configuration file.\cr
If only one suffix is specified all experimental datasets will take the same suffix.\cr
An NA value in the list of suffixes is interpreted as "take the suffix specified in the configuration file".\cr
Ex: c('_f6h', '_3hourly')
  }
  \item{suffixobs}{
Vector of suffixes used to build the dataset path of the observational datasets respectively. Some times paths to datasets can vary depending on slight differences in postprocessings or reanalises adjustments, which leads to an addition of a suffix in the common path pattern for the rest of observational datasets. In such situation 'suffixobs' can help.\cr
Each pair of observational dataset ID and variable name can have a suffix associated in the configuration file.\cr
If only one suffix is specified all observational datasets will take the same suffix.\cr
An NA value in the list of suffixes is interpreted as "take the suffix specified in the configuration file".\cr
Ex: c('_f6h', '_3hourly')
  }
  \item{varmin}{
Loaded experimental and observational data values smaller than 'varmin' will be disabled (replaced by NA values).\cr
Takes by default the value specified in the configuration file.
  }
  \item{varmax}{
Loaded experimental and observational data values greater than 'varmax' will be disabled (replaced by NA values).\cr
Takes by default the value specified in the configuration file.
  }
  \item{silent}{
Parameter to show (FALSE) or hide (TRUE) information messages.\cr
Warnings will be displayed even if 'silent' is set to TRUE.\cr
Takes by default the value 'FALSE'.
  }
  \item{nprocs}{
Number of parallel processes created to perform the fetch and computation of data.\cr
These processes will use shared memory in the processor in which Load() is launched.\cr
By default the number of logical cores in the machine will be detected and as many processes as logical cores there are will be created.\cr
A value of 1 won't create parallel processes.\cr
When running in multiple processes, if an error occurs in any of the processes, a crash message appears in the R session of the original process but no detail is given about the error. A value of 1 will display all error messages in the original and only R session.\cr
Note: the parallel process create other blocking processes each time they need to compute an interpolation via 'cdo'.
  }
  \item{dimnames}{
Named list where the name of each element is a generic name of the expected dimensions inside the NetCDF files. These generic names are 'longitudes', 'latitudes' and 'members'. 'time' is not needed because it's detected automatically by discard.\cr
The value associated to each name is the actual dimension name in the NetCDF file.\cr
The variables in the file that contain the longitudes and latitudes of the data (if the data is a 2-dimensional variable) must have the same name as the longitude and latitude dimensions.\cr
By default, these names are taken from the mandatory variables in the configuration file (DEFAULT_DIM_NAME_LONGIUDE,\cr
DEFAULT_DIM_NAME_LATITUDE and \cr
DEFAULT_DIM_NAME_MEMBERS). If any of those is defined in the 'dimnames' parameter, it takes priority and overwrites the value in the configuration file.
Ex.: list(longitudes = 'x', latitudes = 'y')
In that example, the dimension 'members' will take the \cr
DEFAULT_DIM_NAME_MEMBERS specified in the configuration file.
  }
}
\details{
The two output matrices have between 2 and 6 dimensions:\cr
  \enumerate{
    \item{Number of experimental/observational datasets.}
    \item{Number of members.}
    \item{Number of startdates.}
    \item{Number of leadtimes.}
    \item{Number of latitudes (optional).}
    \item{Number of longitudes (optional).}
  }
but the two matrices have the same number of dimensions and only the first two dimensions can have different lengths depending on the input arguments.
    
For a detailed explanation of the process, read the documentation attached to the package or check the comments in the code.
}
\value{
  \item{$mod}{
Model outputs.\cr
If output = 'areave', matrix with dimensions\cr
  c(nmod/nexp, nmemb/nparam, nsdates, nltime)\cr
If output = 'lat', matrix with dimensions\cr
  c(nmod/nexp, nmemb/nparam, nsdates, nltime, nlat)\cr
If output = 'lon', matrix with dimensions\cr
  c(nmod/nexp, nmemb/nparam, nsdates, nltime, nlon)\cr
If output = 'lonlat', matrix with dimensions\cr
  c(nmod/nexp, nmemb/nparam, nsdates, nltime, nlat, nlon)\cr
  }
  \item{$obs}{
Observations. Matrix with same dimensions as '$mod' except along the first two.\cr
If output = 'areave', matrix with dimensions\cr
  c(nobs, nmemb, nsdates, nltime)\cr
If output = 'lat', matrix with dimensions\cr
  c(nobs, nmemb, nsdates, nltime, nlat)\cr
If output = 'lon', matrix with dimensions\cr
  c(nobs, nmemb nsdates, nltime, nlon)\cr
If output = 'lonlat', matrix with dimensions\cr
  c(nobs, nmemb, nsdates, nltime, nlat, nlon)\cr
  }
  \item{$lat}{
Latitudes of the output grid (default: model grid of the first experiment).\cr
If 'areave' is selected or a global mean variable is specified, takes value 0.
  }
  \item{$lon}{
Longitudes of the output grid (default: model grid of the first experiment).\cr
If 'areave' is selected or a global mean variable is specified, takes value 0.
  }
}
\author{
History:\cr
0.1  -  2011-03  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to CRAN\cr
1.2  -  2015-02  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Generalisation + parallelisation\cr
1.3  -  2015-07  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Improvements related to configuration file mechanism\cr
}
\examples{
# Let's assume we want to perform verification with data of a variable
# called 'tos' from a model called 'model' and observed data coming from 
# an observational dataset called 'observation'.
#
# The model was run in the context of an experiment named 'experiment'. 
# It simulated from 1st November in 1985, 1990, 1995, 2000 and 2005 for a 
# period of 5 years time from each starting date. 5 different sets of 
# initial conditions were used so an ensemble of 5 members was generated 
# for each starting date.
# The model generated values for the variables 'tos' and 'tas' in a 
# 3-hourly frequency but, after some initial post-processing, it was 
# averaged over every month.
# The resulting monthly average series were stored in a file for each 
# starting date for each variable with the data of the 5 ensemble members.
# The resulting directory tree was the following:
#   model
#    |--> experiment
#          |--> monthly_mean
#                |--> tos_3hourly
#                |     |--> tos_19851101.nc
#                |     |--> tos_19901101.nc
#                |               .
#                |               .
#                |     |--> tos_20051101.nc 
#                |--> tas_3hourly
#                      |--> tas_19851101.nc
#                      |--> tas_19901101.nc
#                                .
#                                .
#                      |--> tas_20051101.nc
# 
# The observation recorded values of 'tos' and 'tas' at each day of the 
# month over that period but was also averaged over months and stored in 
# a file per month. The directory tree was the following:
#   observation
#    |--> monthly_mean
#          |--> tos
#          |     |--> tos_198511.nc
#          |     |--> tos_198512.nc
#          |     |--> tos_198601.nc
#          |               .
#          |               .
#          |     |--> tos_201010.nc
#          |--> tas
#                |--> tas_198511.nc
#                |--> tas_198512.nc
#                |--> tas_198601.nc
#                          .
#                          .
#                |--> tas_201010.nc
#
# The model data is stored in a file-per-startdate fashion and the
# observational data is stored in a file-per-month, and both are stored in 
# a monthly frequency. The file format is NetCDF.
# Hence all the data is supported by Load() (see details and other supported 
# conventions in ?Load) but first we need to configure it properly.
#
# These data files are included in the package (in the 'sample_data' folder),
# only for the variable 'tos'. They have been interpolated to a very low 
# resolution grid so as to make it on CRAN.
# The original grid names (following CDO conventions) for experimental and 
# observational data were 't106grid' and 'r180x89' respectively. The final
# resolutions are 'r20x10' and 'r16x8' respectively. 
# The experimental data comes from the decadal climate prediction experiment 
# run at IC3 in the context of the CMIP5 project. Its name within IC3 local 
# database is 'i00k'. 
# The observational dataset used for verification is the 'ERSST' 
# observational dataset.
#
# The configuration file 'sample.conf' that we will create in the example 
# has the proper entries to load these (see ?LoadConfigFile for details on 
# writing a configuration file). 
#
# The code is not run because it dispatches system calls to 'cdo' and 'nco'
# which is not allowed as per CRAN policies. You can run it in your system 
# though. Instead, the code in 'dontshow' is run, which loads the equivalent
# data already processed in R.
  \dontrun{
configfile <- paste0(tempdir(), '/sample.conf')
ConfigFileCreate(configfile, confirm = FALSE)
c <- ConfigFileOpen(configfile)
c <- ConfigEditDefinition(c, 'DEFAULT_GRID', 'r20x10', confirm = FALSE)
c <- ConfigEditDefinition(c, 'DEFAULT_VAR_MIN', '-1e19', confirm = FALSE)
c <- ConfigEditDefinition(c, 'DEFAULT_VAR_MAX', '1e19', confirm = FALSE)
c <- ConfigAddVar(c, '2d', 'tos')
data_path <- system.file('sample_data', package = 's2dverification')
exp_data_path <- paste0(data_path, '/model/$EXP_NAME$')
obs_data_path <- paste0(data_path, '/$OBS_NAME$')
c <- ConfigAddEntry(c, 'experiments', 'file-per-startdate',
     dataset_name = 'experiment', var_name = 'tos', main_path = exp_data_path,
     file_path = '$STORE_FREQ$_mean/$VAR_NAME$_3hourly/$VAR_NAME$_$START_DATE$.nc')
c <- ConfigAddEntry(c, 'observations', 'file-per-month',
     dataset_name = 'observation', var_name = 'tos', main_path = obs_data_path,
     file_path = '$STORE_FREQ$_mean/$VAR_NAME$/$VAR_NAME$_$YEAR$$MONTH$.nc')
ConfigFileSave(c, configfile, confirm = FALSE)

# Now we are ready to use Load().
startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
sampleData <- Load('tos', c('experiment'), c('observation'), startDates, 
                   output = 'areave', latmin = 27, latmax = 48, 
                   lonmin = -12, lonmax = 40, configfile = configfile)
  }
  \dontshow{
startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
sampleData <- s2dverification:::.LoadSampleData('tos', c('experiment'), 
                                                c('observation'), startDates,
                                                output = 'areave', 
                                                latmin = 27, latmax = 48, 
                                                lonmin = -12, lonmax = 40) 
  } 
}
\keyword{datagen}
