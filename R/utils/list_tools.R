import::here(stringi, 'stri_replace_all_regex')

## Functions
## items_in_a_not_b
## filter_list_for_match
## find_first_match_index
## replace_specific_items
## multiple_replacement


#' https://stackoverflow.com/questions/10298662/find-elements-not-in-smaller-character-vector-list-but-in-big-list
#' 
#' @export
items_in_a_not_b <- function(a, b) {
    return((new <- a[which(!a %in% b)]))
}


#' return elements of a list matching a particular substring
#'
#' @examples
#' filter_list_for_match(c("gene_id_pat", "gene_id_mat", "count"), "pat")
#' 
#' @export
filter_list_for_match <- function(items, patterns) {
    # filter
    for (i in 1:length(patterns)){
        items <- lapply(items, grep, pattern=patterns[[i]], value=TRUE)
    }
    return (unlist(items[!sapply(items, identical, character(0))]))  # remove character(0)
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
#' replace_dict <- c(
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
