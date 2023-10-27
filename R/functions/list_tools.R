## Functions
## items_in_a_not_b
## find_first_match_index
## replace_specific_items
## multiple_replacement


#' https://stackoverflow.com/questions/10298662/find-elements-not-in-smaller-character-vector-list-but-in-big-list
#' 
#' @export
items_in_a_not_b <- function(a, b) {
    return((new <- a[which(!a %in% b)]))
}


#' return the first match given a pattern
#'
#' @export
find_first_match_index <- function(pattern, items) {
    return (grep(pattern, items)[[1]])
}


#' Replace specific items
#' Use this to rename columns
#'
#' @export
replace_specific_items <- function(items, replacer) {
    replace_ids <- which(items %in% intersect(names(replacer), items))
    for (idx in replace_ids) {
        items[idx] <- replacer[items[idx]]
    }
    return(items)
}


#' Convenience function to perform multiple replacements on a list
#' 
#' Example:
#' replacement_dict <- c(
#'     '[A-Za-z]' = '',
#'     '[0-9]+' = ''
#' )
#' 
#' @export
multiple_replacement <- function(items, replace_dict, func='gsub') {

    if (func == 'gsub') {
        replace_func = gsub
    } else if (func == 'sub') {
        replace_func = sub
    } else {
        return(items)
    }

    for (pattern in names(replace_dict)) {
        replacement = replace_dict[[pattern]]
        items <- unlist( lapply(items, function(x) replace_func(pattern, replacement, x)) )
    }

    return (items)
}
