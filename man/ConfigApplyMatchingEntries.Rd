% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConfigApplyMatchingEntries.R
\name{ConfigApplyMatchingEntries}
\alias{ConfigApplyMatchingEntries}
\title{Apply Matching Entries To Dataset Name And Variable Name To Find Related Info}
\usage{
ConfigApplyMatchingEntries(
  configuration,
  var,
  exp = NULL,
  obs = NULL,
  show_entries = FALSE,
  show_result = TRUE
)
}
\arguments{
\item{configuration}{Configuration object obtained from ConfigFileOpen() 
or ConfigFileCreate().}

\item{var}{Name of the variable to load. Will be interpreted as a string, 
regular expressions do not apply here.
Examples: 'tas' or 'tasmax_q90'.}

\item{exp}{Set of experimental dataset identifiers. Will be interpreted as 
a strings, regular expressions do not apply here. Can be NULL (not to 
check in experimental dataset tables), and takes by default NULL.
Examples: c('EnsEcmwfSeas', 'EnsUkmoSeas'), c('i00k').}

\item{obs}{Set of observational dataset identifiers. Will be interpreted as 
a strings, regular expressions do not apply here. Can be NULL (not to 
check in observational dataset tables), and takes by default NULL.
Examples: c('GLORYS', 'ERAint'), c('NCEP').}

\item{show_entries}{Flag to stipulate whether to show the found matching 
entries for all datasets and variable name.}

\item{show_result}{Flag to stipulate whether to show the result of applying 
all the matching entries (dataset main path, file path, ...).}
}
\value{
A list with the information resulting of applying the matching 
 entries is returned.
}
\description{
Given a pair of dataset name and variable name, this function determines 
applies all the matching entries found in the corresponding configuration 
table to work out the dataset main path, file path, actual name of variable 
inside NetCDF files, ...
}
\examples{
# Create an empty configuration file
config_file <- paste0(tempdir(), "/example.conf")
s2dverification:::ConfigFileCreate(config_file, confirm = FALSE)
# Open it into a configuration object
configuration <- ConfigFileOpen(config_file)
# Add an entry at the bottom of 4th level of file-per-startdate experiments 
# table which will associate the experiment "ExampleExperiment2" and variable 
# "ExampleVariable" to some information about its location.
configuration <- ConfigAddEntry(configuration, "experiments", 
                "last", "ExampleExperiment2", "ExampleVariable", 
                "/path/to/ExampleExperiment2/", 
                "ExampleVariable/ExampleVariable_$START_DATE$.nc")
# Edit entry to generalize for any variable. Changing variable needs .
configuration <- ConfigEditEntry(configuration, "experiments", 1, 
                var_name = ".*", 
                file_path = "$VAR_NAME$/$VAR_NAME$_$START_DATE$.nc")
# Now apply matching entries for variable and experiment name and show the 
# result
match_info <- ConfigApplyMatchingEntries(configuration, 'tas', 
             exp = c('ExampleExperiment2'), show_result = TRUE)
}
\seealso{
ConfigApplyMatchingEntries, ConfigEditDefinition, 
 ConfigEditEntry, ConfigFileOpen, ConfigShowSimilarEntries, 
 ConfigShowTable
}
\author{
History:\cr
 0.1 - 2015-05 (N. Manubens) - First version\cr
 1.0 - 2015-11 (N. Manubens) - Removed grid column and storage types
}
\keyword{datagen}
