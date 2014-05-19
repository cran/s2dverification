LeapYear <- function(year) {
  # This function tells whether a year is leap year or not.
  #
  # Args:
  #   year: The year to tell whether is leap year or not.
  #
  # Returns:
  #   Boolean telling whether the year is a leap year or not.
  #
  # History:
  #   1.0  #  2011-03  (V. Guemas, vguemas@ic3.cat)  #  Original code

  leap <- FALSE
  if (year %% 4 == 0) {
    leap <- TRUE
    if (year %% 100 == 0) {
      leap <- FALSE
      if (year %% 400 == 0) {
        leap <- TRUE
      }
    } 
  }
  
  #
  # Output
  #
  leap
}
