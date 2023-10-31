import::from('stringi', 'stri_replace_all_regex', .character_only=TRUE)

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


#' Convenience function to perform multiple replacements on a list or dataframe column
#' 
#' Example:
#' replacements <- c(
#'     '[A-Za-z]' = '',
#'     '[0-9]+' = ''
#' )
#' 
#' @export
multiple_replacement <- function(items, replace_dict) {

    # for (pattern in names(replace_dict)) {
    #     replacement = replace_dict[[pattern]]
    #     items <- sapply(items, function(x) gsub(pattern, replacement, x))
    # }

    patterns <- names(replace_dict)
    replacements <- sapply(unname(replace_dict), function(x) gsub('\\\\', '$', x))
    
    items <- sapply(items,
        function(x) stri_replace_all_regex(
            x,
            pattern = patterns,
            replacement = replacements,
            vectorize_all = FALSE)
    )

    return (items)
}
