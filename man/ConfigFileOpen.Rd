% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConfigFileOpen.R
\name{ConfigFileOpen}
\alias{ConfigFileOpen}
\alias{ConfigFileCreate}
\alias{ConfigFileSave}
\title{Functions To Create Open And Save Configuration File}
\usage{
ConfigFileOpen(file_path, silent = FALSE, stop = FALSE)

ConfigFileCreate(file_path, confirm = TRUE)

ConfigFileSave(configuration, file_path, confirm = TRUE)
}
\arguments{
\item{file_path}{Path to the configuration file to create/open/save.}

\item{silent}{Flag to activate or deactivate verbose mode.
Defaults to FALSE (verbose mode on).}

\item{stop}{TRUE/FALSE whether to raise an error if not all the mandatory 
default variables are defined in the configuration file.}

\item{confirm}{Flag to stipulate whether to ask for confirmation when 
saving a configuration file that already exists.\cr
Defaults to TRUE (confirmation asked).}

\item{configuration}{Configuration object to save in a file.}
}
\value{
ConfigFileOpen() returns a configuration object with all the information for 
 the configuration file mechanism to work.\cr
ConfigFileSave() returns TRUE if the file has been saved and FALSE otherwise.\cr
ConfigFileCreate() returns nothing.
}
\description{
These functions help in creating, opening and saving configuration files.
}
\details{
ConfigFileOpen() loads all the data contained in the configuration file 
specified as parameter 'file_path'.
Returns a configuration object with the variables needed for the 
configuration file mechanism to work.
This function is called from inside the Load() function to load the 
configuration file specified in 'configfile'.\cr\cr
ConfigFileCreate() creates an empty configuration file and saves it to 
the specified path. It may be opened later with ConfigFileOpen() to be edited.
Some default values are set when creating a file with this function, you 
can check these with ConfigShowDefinitions().\cr\cr
ConfigFileSave() saves a configuration object into a file, which may then 
be used from Load().\cr\cr
Two examples of configuration files can be found inside the 'inst/config/' 
folder in the package:
 \itemize{
   \item{BSC.conf: configuration file used at BSC-CNS. Contains location 
   data on several datasets and variables.}
   \item{template.conf: very simple configuration file intended to be used as 
   pattern when starting from scratch.}
 }
How the configuration file works:\cr
~~~~~~~~~~~~~~~~~~~~~~~~~~~~\cr
It contains one list and two tables.\cr
Each of these have a header that starts with '!!'. These are key lines and 
should not be removed or reordered.\cr
Lines starting with '#' and blank lines will be ignored.
The list should contains variable definitions and default value definitions.\cr
The first table contains information about experiments.\cr
The third table contains information about observations.\cr
Each table entry is a list of comma-separated elements.\cr
The two first are part of a key that is associated to a value formed by the 
other elements.\cr
The key elements are a dataset identifier and a variable name.\cr
The value elements are the dataset main path, dataset file path, the 
variable name inside the .nc file, a default suffix (explained below) and a 
minimum and maximum vaues beyond which loaded data is deactivated.\cr
Given a dataset name and a variable name, a full path is obtained 
concatenating the main path and the file path.\cr
Also the nc variable name, the suffixes and the limit values are obtained.\cr
Any of the elements in the keys can contain regular expressions[1] that will 
cause matching for sets of dataset names or variable names.\cr
The dataset path and file path can contain shell globbing expressions[2] 
that will cause matching for sets of paths when fetching the file in the 
full path.\cr
The full path can point to an OPeNDAP URL.\cr
Any of the elements in the value can contain variables that will be replaced 
to an associated string.\cr
Variables can be defined only in the list at the top of the file. \cr
The pattern of a variable definition is\cr
VARIABLE_NAME = VARIABLE_VALUE\cr
and can be accessed from within the table values or from within the variable 
values as\cr
 $VARIABLE_NAME$\cr
For example:\cr
 FILE_NAME = tos.nc\cr
 !!table of experiments\cr
 ecmwf, tos, /path/to/dataset/, $FILE_NAME$\cr
There are some reserved variables that will offer information about the 
store frequency, the current startdate Load() is fetching, etc:\cr
 $VAR_NAME$, $START_DATE$, $STORE_FREQ$, $MEMBER_NUMBER$\cr
 for experiments only: $EXP_NAME$\cr
 for observations only: $OBS_NAME$, $YEAR$, $MONTH$, $DAY$\cr
Additionally, from an element in an entry value you can access the other 
elements of the entry as:\cr
 $EXP_MAIN_PATH$, $EXP_FILE_PATH$, \cr$VAR_NAME$, $SUFFIX$, $VAR_MIN$, $VAR_MAX$\cr
\cr
The variable $SUFFIX$ is useful because it can be used to take part in the 
main or file path. For example: '/path/to$SUFFIX$/dataset/'.\cr
It will be replaced by the value in the column that corresponds to the 
suffix unless the user specifies a different suffix via the parameter 
'suffixexp' or 'suffixobs'.\cr
This way the user is able to load two variables with the same name in the 
same dataset but with slight modifications, with a suffix anywhere in the 
path to the data that advices of this slight modification.\cr\cr
The entries in a table will be grouped in 4 levels of specificity:
 \enumerate{
   \item{
General entries:\cr
- the key dataset name and variable name are both a regular expression 
matching any sequence of characters (.*) that will cause matching for any 
pair of dataset and variable names\cr
  Example:  .*, .*, /dataset/main/path/, file/path, nc_var_name, suffix, 
var_min, var_max
   }
   \item{
Dataset entries:\cr
- the key variable name matches any sequence of characters\cr
  Example:  ecmwf, .*, /dataset/main/path/, file/path, nc_var_name, 
 suffix, var_min, var_max
   }
   \item{
Variable entries:\cr
- the key dataset name matches any sequence of characters\cr
  Example:  .*, tos, /dataset/main/path/, file/path, nc_var_name, 
 suffix, var_min, var_max
   }
   \item{
 Specific entries:\cr
- both key values are specified\cr
  Example:  ecmwf, tos, /dataset/main/path/, file/path, nc_var_name, 
 suffix, var_min, var_max
   }
 }
Given a pair of dataset name and variable name for which we want to know the 
full path, all the rules that match will be applied from more general to 
more specific.\cr
If there is more than one entry per group that match a given key pair, 
these will be applied in the order of appearance in the configuration file 
(top to bottom).\cr\cr
An asterisk (*) in any value element will be interpreted as 'leave it as is 
or take the default value if yet not defined'.\cr
The default values are defined in the following reserved variables:\cr
 $DEFAULT_EXP_MAIN_PATH$, $DEFAULT_EXP_FILE_PATH$, $DEFAULT_NC_VAR_NAME$, 
$DEFAULT_OBS_MAIN_PATH$, $DEFAULT_OBS_FILE_PATH$, $DEFAULT_SUFFIX$, 
$DEFAULT_VAR_MIN$, $DEFAULT_VAR_MAX$, \cr
$DEFAULT_DIM_NAME_LATITUDES$, $DEFAULT_DIM_NAME_LONGITUDES$, \cr
$DEFAULT_DIM_NAME_MEMBERS$\cr\cr
Trailing asterisks in an entry are not mandatory. For example\cr
 ecmwf, .*, /dataset/main/path/, *, *, *, *, *\cr
will have the same effect as\cr
 ecmwf, .*, /dataset/main/path/ \cr\cr
A double quote only (") in any key or value element will be interpreted as 
'fill in with the same value as the entry above'.
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
# Finally save the configuration file.
ConfigFileSave(configuration, config_file, confirm = FALSE)

}
\references{
[1] \url{https://stat.ethz.ch/R-manual/R-devel/library/base/html/regex.html}\cr
[2] \url{https://tldp.org/LDP/abs/html/globbingref.html}
}
\seealso{
ConfigApplyMatchingEntries, ConfigEditDefinition, 
 ConfigEditEntry, ConfigFileOpen, ConfigShowSimilarEntries, ConfigShowTable
}
\author{
History:
 0.1 - 2015-05 (N. Manubens) - First version
 1.0 - 2015-11 (N. Manubens) - Removed grid column and storage formats
}
\keyword{datagen}
