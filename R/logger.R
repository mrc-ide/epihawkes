#' Logger function
#'
#' Prints messages to screen at different levels of logging
#' @param msg Message to be printed.
#' @param log_level Level of the logger at which this message should be printed.
#' @param print_level The level the code is running at
#' @return Message printed to console.
#' @export
print_message <- function(msg, log_level, print_level){
  if (log_level <= print_level)
    print(sprintf(msg))
}
