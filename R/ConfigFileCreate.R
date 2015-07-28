ConfigFileCreate <- function(file_path, confirm = TRUE) {
  success <- ConfigFileSave(list(definitions = list(
    DEFAULT_EXP_MAIN_PATH = "$EXP_NAME$", DEFAULT_EXP_FILE_PATH = "$STORE_FREQ$/$VAR_NAME$_$START_DATE$.nc",
    DEFAULT_GRID = "t106grid", DEFAULT_NC_VAR_NAME = "$VAR_NAME$", 
    DEFAULT_SUFFIX = "", DEFAULT_VAR_MIN = "", 
    DEFAULT_VAR_MAX = "", EXP_FULL_PATH = "$EXP_MAIN_PATH$/$EXP_FILE_PATH$", 
    DEFAULT_OBS_MAIN_PATH = "$OBS_NAME$", 
    DEFAULT_OBS_FILE_PATH = "$STORE_FREQ$/$VAR_NAME$_$YEAR$$MONTH$.nc", 
    OBS_FULL_PATH = "$OBS_MAIN_PATH$/$OBS_FILE_PATH$",
    DEFAULT_DIM_NAME_LONGITUDES = "longitude", DEFAULT_DIM_NAME_LATITUDES = "latitude", 
    DEFAULT_DIM_NAME_MEMBERS = "ensemble")), file_path, confirm = confirm)
  if (success) {
    cat("WARNING: You have just created a configuration file. You can edit the defaults according to your needs with the functions ConfigFileOpen(), ConfigEditDefinition() and ConfigFileSave() or edit the file manually as specified in ?ConfigFileOpen.\n")
  }
}
