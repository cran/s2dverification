\name{ConfigAddVar}
\alias{ConfigAddVar}
\alias{ConfigRemoveVar}
\title{
Functions To Add Or Remove Supported Variables
}
\description{
These functions add or remove supported variables to or from the lists of supported 2d or global-mean variables.
}
\usage{
ConfigAddVar(configuration, var_type, var_entry)
ConfigRemoveVar(configuration, var_type, var_entry)
}
\arguments{
  \item{configuration}{
Configuration object obtained from ConfigFileOpen() or ConfigFileCreate() that contains the lists of variables that will be extended.
  }
  \item{var_type}{
Name of the list of variables to modify. Can take values '2d' if willing to modify the two-dimensional variable list, or 'mean' if willing to modify global mean variable list.
  }
  \item{var_entry}{
Name of the variable to add or remove (can be a regular expression to match sets of variables. See ?ConfigFileOpen).
  }
}
\value{
A modified configuration object is returned.
}
\seealso{ConfigAddVar, ConfigApplyMatchingEntries, ConfigEditDefinition, ConfigEditEntry, ConfigFileOpen, ConfigMatchingVars, ConfigShowSimilarEntries, ConfigShowTable}
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
# table which will associate the experiment "ExampleExperiment" and variable 
# "ExampleVariable" to some information about its location.
configuration <- ConfigAddEntry(configuration, "experiments", 
                 "file-per-startdate", "last", "ExampleExperiment", 
                 "ExampleVariable", "/path/to/ExampleExperiment", 
                 "ExampleVariable/ExampleVariable_$START_DATE$.nc")
# Add another entry
configuration <- ConfigAddEntry(configuration, "experiments", 
                 "file-per-startdate", "last", "ExampleExperiment2", 
                 "ExampleVariable", "/path/to/ExampleExperiment2", 
                 "ExampleVariable/ExampleVariable_$START_DATE$.nc")
# Edit second entry to generalize for any variable. Changing variable needs .
configuration <- ConfigEditEntry(configuration, "experiments", 
                 "file-per-startdate", 2, var_name = ".*", 
                 file_path = "$VAR_NAME$/$VAR_NAME$_$START_DATE$.nc")
# Remove first entry
configuration <- ConfigRemoveEntry(configuration, "experiments", 
                 "file-per-startdate", "ExampleExperiment", "ExampleVariable")
# Show results
ConfigShowTable(configuration, "experiments", "file-per-ensemble")
# Now put variable in list of supported variables. In this example will be a 
# two-dimensional variable.
configuration <- ConfigAddVar(configuration, "2d", "ExampleMistake")
# Look in list of variables and see there's nothing similar to what we look 
# for. The previous introduction was mistaken
configuration <- ConfigShowSimilarVars(configuration, "Variable")
# Remove variable and add the desired one
configuration <- ConfigRemoveVar(configuration, "2d", "ExampleMistake")
configuration <- ConfigAddVar(configuration, "2d", "ExampleVariable")
# See results
ConfigShowVars(configuration)
# Save the configuration
ConfigFileSave(configuration, config_file, confirm = FALSE)
}
\keyword{datagen}
