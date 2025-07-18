import::here(stringi, 'stri_replace_all_regex')

## Functions

## dict_zip
## fill_missing_keys
## filter_list_for_match
## find_first_match_index
## interleave
## items_in_a_not_b
## list2matrix
## matrix2list
## move_list_items_to_front
## multiple_replacement
## replace_specific_items


#' Convert a matrix into a list of lists
#'
#' Input:
#'      [,1] [,2] [,3] [,4] [,5] [,6]
#' [1,]    1    1    1    2    2    3
#' [2,]    2    3    4    3    4    4
#' 
#' Output: 
#' list(c(1, 2), c(1, 3), c(1, 4), c(2, 4), c(2, 4), c(3, 4))
#' 
collect_matrix_cols <- function(mat) {
    return(split(mat, col(mat)))
}


#' Dictionary
#' 
#' @description Simple dictionary implementation using R environment
#'
#' @export
dict_zip <- function(keys, values) {
    if (!is.list(values)) {
        values <- as.list(values)
    }
    named_list <- setNames(values, keys)
    env <- list2env(named_list)
    return(env)
}


#' Fill Missing Keys
#' 
#' Add missing keys to a named list (ie. dict)
#' 
fill_missing_keys <- function(x, keys, val=NA) {

    missing <- setdiff(keys, names(x))
    x[missing] <- val  # fill
    x <- x[keys]  # sort
    return(x)
}


#' Find matches based on substring
#' 
#' @description Return elements of a list matching a particular substring
#' 
#' @param items list or vector
#' @param patterns a string or collection of strings
#' @return Returns a list or vector with any matching items
#' 
#' @examples
#' filter_list_for_match(c("a_suffix", "b_suffix", "c"), "suffix")
#' 
filter_list_for_match <- function(items, patterns) {
    # filter
    for (i in 1:length(patterns)){
        items <- lapply(items, grep, pattern=patterns[[i]], value=TRUE)
    }
    return (unlist(items[!sapply(items, identical, character(0))]))  # remove character(0)
}


#' Returns the first match given a pattern
#' 
#' @description Thin wrapper around [grep()]
#' 
#' @param items list or vector
#' @param pattern a single regex pattern, can also be an exact string
#' @return Returns the index of the first match
#' 
#' @examples
#' find_first_match_index('_suffix', c('a', 'b_suffix', 'c'))
#' 
find_first_match_index <- function(pattern, items) {
    return (grep(pattern, items)[[1]])
}


#' Interleave
#' 
#' @description
#' 
#' @param a list or vector
#' @param b list or vector
#' @return Returns the interleaved list
#' 
#' @references
#' \href{https://stackoverflow.com/questions/16443260/interleave-lists-in-r}{Stack Overflow}
#' 
#' @examples
#' interleave(c('A', 'B', 'C'), c('a', 'b', 'c'))
#' 
interleave <- function(a, b) {
    idx <- order(c(seq_along(a), seq_along(b)))
    return(unlist(c(a,b))[idx])
}


#' Return items unique to vector a
#' 
#' @description
#' Return all the items found in a and not b. Useful for filtering dataframe columns.
#' 
#' @param a list or vector
#' @param b list or vector
#' @return Returns the filtered list or vector (same type as a)
#' 
#' @references
#' \href{https://stackoverflow.com/questions/10298662/find-elements-not-in-smaller-character-vector-list-but-in-big-list}{Stack Overflow}
#' 
#' @examples
#' items_in_a_not_b(c('a', 'b', 'c', '1', '2'), c('1', '2'))
#' 
items_in_a_not_b <- function(a, b) {
    return((new <- a[which(!a %in% b)]))
}


#' Rearrange a list into a matrix
#' 
#' @description Thin wrapper around matrix, inserts NA
#' 
#' @param z list or vector
#' @param ncol
#' @param byrow
#' @param fill
#' @return 
#' 
#' @examples
#' list2matrix(c('a', 'b', 'c', 'd', 'e'), ncol=3 fill=NA)
#' 
list2matrix <- function(z, ncol=2, byrow=TRUE, fill=NA) {
    remain <- length(z) %% ncol
    if (remain > 0) {
        z <- c(z, rep(fill, ncol-remain))
    }
    mtx <- matrix(z, ncol = ncol, byrow = byrow)
    return(mtx)
}


#' Rearrange a matrix to a named list
#' 
matrix2list <- function(mat) {
    vec <- as.vector(mat)
    names(vec) <- with(
        expand.grid(rownames(mat), colnames(mat)), 
        paste(Var1, Var2, sep = "-")
    )
    return(unlist(vec))
}


#' Rearrange items according to an ordered list
#' 
#' @description Front-loads the ordered list, then adds the remaining items
#' 
#' @param items list or vector
#' @param ordered
#' @return 
#' 
#' @examples
#' move_list_items_to_front(c('a', 'b', 'c', 'd'), c('b', 'c'))
#' 
move_list_items_to_front <- function(items, ordered) {
    return(c(ordered[ordered %in% items],
             items_in_a_not_b(items, ordered)))
}


#' Replaces each item in a list or vector with all replacements
#' 
#' @description
#' Convenience function to perform multiple replacements on a list or dataframe column.
#' Unlike [replace_specific_items()], `multiple_replacement()` can recognize patterns.
#' Note: this is equivalent to [plyr::mapvalues()]
#' 
#' @param items list or vector
#' @param replacements a named list of replacements. uses names to match and values to replace.
#' @return Returns a vector with replaced items
#' 
#' @examples
#' replacements <- c('prefix_' = '', '_suffix' = '')
#' items <- c('prefix_a_suffix', 'prefix_b_suffix')
#' multiple_replacement(items, replacements)
#' 
#' @seealso [replace_specific_items()], [stringi::stri_replace_all_regex()]
#' 
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


#' Replaces each item in a list or vector with all replacements
#' 
#' @description Use this to rename columns
#' 
#' @param items list or vector
#' @param replacements a named list of replacements. uses names to match and values to replace.
#' @return Returns a vector with replaced items
#' 
#' @examples
#' replace_specific_items(c('a', 'b', 'c'), c('a'="A", 'c'="C"))
#' 
#' @seealso [multiple_replacement()]
#' 
replace_specific_items <- function(items, replacer) {
    replace_ids <- which(items %in% intersect(names(replacer), items))
    for (idx in replace_ids) {
        items[idx] <- replacer[items[idx]]
    }
    return(items)
}
