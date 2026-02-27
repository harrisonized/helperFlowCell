## Functions
## format_number
## title_to_snake_case
## substr_right
## txt_strip


#' Format Number
#' 
#' Truncate to two digits for plotting
#' 
#' 1500000000 -> "1.5B"
#' 50000000   -> "50M"
#' 999999     -> "1000K"
#' 5432.897   -> "5.43K"
#' 123.456    -> "123.46"
#' 0.5        -> "0.5"
#' -2500000   -> "-2.5M"
#' 
format_number <- function(x, sigfigs = 2) {
  abs_x <- abs(x)
  sign  <- ifelse(x < 0, "-", "")

  scaled <- ifelse(abs_x >= 1e9, abs_x / 1e9,
            ifelse(abs_x >= 1e6, abs_x / 1e6,
            ifelse(abs_x >= 1e3, abs_x / 1e3,
                                 abs_x)))

  suffix <- ifelse(abs_x >= 1e9, "B",
            ifelse(abs_x >= 1e6, "M",
            ifelse(abs_x >= 1e3, "K",
                                 "")))

  # Apply sigfigs after scaling so significant digits are counted correctly
  scaled <- signif(scaled, sigfigs)

  # Use enough decimal places to show all sig figs
  digits <- pmax(0, sigfigs - ceiling(log10(pmax(abs(scaled), 1))))
  formatted <- mapply(function(s, d) formatC(s, format = "f", digits = d), scaled, digits)
  formatted <- sub("\\.?0+$", "", formatted)

  return(paste0(sign, formatted, suffix))
}


#' Standardize Space-separated Titles
#' 
#' @description Converts "Column Title" to column_title
#' 
#' @examples
#' title_to_snake_case('Column Title')
#' 
title_to_snake_case <- function(text) {
    return(tolower(
        paste(
            unlist(strsplit(text, '[ ]')), collapse='_')
        )
    )
}


#' Extract last n characters
#' 
#' @description Extracts the last n characters from a string
#' 
#' @examples
#' substr_right('hi there', 5)
#' 
#' @references
#' \href{https://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r}{StackOverflow post}
#'
substr_right <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


#' Remove leading and traling characters
#' 
#' @description Removes special characters from beginning and end of a string
#' 
#' @examples
#' txt_strip(" hi there ", chars=" ")
#' 
txt_strip <- function(x, chars=' ') {
    chars <- unique(strsplit(chars, '')[[1]])

    for (char in chars) {

        # special characters
        if (char %in% c('(', ')')) {
            char <- paste0('\\', char)
        }

        x <- gsub(paste0('.*^', char, '+'), '', x)
        x <- gsub(paste0('+', char, '$'), '', x)
    }
    return(x)
}
