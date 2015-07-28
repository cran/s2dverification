\name{ConfigMatchingVars}
\alias{ConfigMatchingVars}
\title{
Find Entries In List Of Supported Variables That Match A Variable
}
\description{
Given a specific variable name, this function checks and shows if it matches with any of the patterns in the lists of supported variables.
}
\usage{
ConfigMatchingVars(configuration, var, show = TRUE)
}
\arguments{
  \item{configuration}{
Configuration object obtained from ConfigFileCreate() or ConfigFileOpen().
  }
  \item{var}{
Variable name which to check if is supported or not.
  }
  \item{show}{
Boolean to stipulate whether to show the matching entries or not.\cr
Takes TRUE by default.
  }
}
\value{
A boolean is returned. It tells whether the variable is supported as a 2-dimensional variable (TRUE) or as a global mean variable (FALSE).
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
# Check if some variables are supported after adding this rule
ConfigMatchingVars(configuration, "Example10Variable")
ConfigMatchingVars(configuration, "ExampleSomeVariable")
ConfigMatchingVars(configuration, "ExampleVariable")
}
\keyword{datagen}
