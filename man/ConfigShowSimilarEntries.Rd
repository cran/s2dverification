\name{ConfigShowSimilarEntries}
\alias{ConfigShowSimilarEntries}
\alias{ConfigShowSimilarVars}
\title{
Find Similar Entries In Lists Of Variables And Tables Of Datasets
}
\description{
These functions help in finding similar entries in lists of supported variables or in tables of supported datasets by comparing all entries with some given information.\cr
This is useful when dealing with complex configuration files and not sure if already support certain variables or datasets.\cr
At least one field must be provided in ConfigShowSimilarEntries(). Other fields can be unspecified and won't be taken into account. If more than one field is provided, sameness is avreaged over all provided fields and entries are sorted from higher average to lower.
}
\usage{
ConfigShowSimilarEntries(configuration, dataset_name = NULL, var_name = NULL, 
                         main_path = NULL, file_path = NULL, nc_var_name = NULL,
                         grid = NULL, suffix = NULL, varmin = NULL, 
                         varmax = NULL, n_results = 10)
ConfigShowSimilarVars(configuration, var_name, n_results = 5)
}
\arguments{
  \item{configuration}{
Configuration object obtained either from ConfigFileCreate() or ConfigFileOpen().
  }
  \item{dataset_name}{
Optional dataset name to look for similars of.
  }
  \item{var_name}{
Optional variable name to look for similars of.
  }
  \item{main_path}{
Optional main path to look for similars of.
  }
  \item{file_path}{
Optional file path to look for similars of.
  }
  \item{nc_var_name}{
Optional variable name inside NetCDF file to look for similars of.
  }
  \item{grid}{
Optional grid to look for similars of.
  }
  \item{suffix}{
Optional suffix to look for similars of.
  }
  \item{varmin}{
Optional variable minimum to look for similars of.
  }
  \item{varmax}{
Optional variable maximum to look for similars of.
  }
  \item{n_results}{
Top 'n_results' alike results will be shown only.\cr
Defaults to 10 in ConfigShowSimilarEntries() and to 5 in ConfigShowSimilarVars().
  }
}
\details{
Sameness is calculated with string distances as specified by Simon White in [1].
}
\value{
These functions return information about the found matches.
}
\seealso{ConfigAddVar, ConfigApplyMatchingEntries, ConfigEditDefinition, ConfigEditEntry, ConfigFileOpen, ConfigMatchingVars, ConfigShowSimilarEntries, ConfigShowTable}
\references{
[1] Simon White, string seamness: http://www.catalysoft.com/articles/StrikeAMatch.html
}
\author{
History:\cr
0.1  -  2015-05  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  First version
}
\examples{
# Create an empty configuration file
config_file <- paste0(tempdir(), "/example.conf")
ConfigFileCreate(config_file, confirm = FALSE)
# Open it into a configuration object
configuration <- ConfigFileOpen(config_file)
# Add an entry at the bottom of 4th level of file-per-startdate experiments 
# table which will associate the experiment "ExampleExperiment2" and variable 
# "ExampleVariable" to some information about its location.
configuration <- ConfigAddEntry(configuration, "experiments", 
                 "file-per-startdate", "last", "ExampleExperiment2", 
                 "ExampleVariable", "/path/to/ExampleExperiment2", 
                 "ExampleVariable/ExampleVariable_$START_DATE$.nc")
# Edit entry to generalize for any variable. Changing variable needs .
configuration <- ConfigEditEntry(configuration, "experiments", 
                 "file-per-startdate", 1, var_name = ".*", 
                 file_path = "$VAR_NAME$/$VAR_NAME$_$START_DATE$.nc")
# Now put a pattern to the list of supported variables. In this example will 
# be a two-dimensional variable.
configuration <- ConfigAddVar(configuration, "2d", "Example.*Variable")
# Look for similar entries
ConfigShowSimilarEntries(configuration, dataset_name = "Exper", 
                         var_name = "Vari")
ConfigShowSimilarVars(configuration, "Vari")
}
\keyword{datagen}
