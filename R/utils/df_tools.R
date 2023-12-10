import::here(file.path(wd, 'R', 'utils', 'list_tools.R'),
    'items_in_a_not_b', 'replace_specific_items', .character_only=TRUE)

## Functions
## rename_columns
## fillna
## reset_index
## append_dataframe
## dataframe_row_from_named_list
## stranspose


#' rename specific dataframe columns
#' 
#' @export
rename_columns <- function(df, columns, inplace=FALSE) {
    df_name <- deparse(substitute(df))
    colnames(df) <- replace_specific_items(colnames(df), columns)
    if (inplace) {
        # see: https://stackoverflow.com/questions/3969852/update-data-frame-via-function-doesnt-work
        assign(df_name, df, envir=.GlobalEnv)
    } else {
        return(df)
    }
}


#' fill a specific column with na
#' 
#' @export
fillna <- function(df, cols, val=0, inplace=FALSE) {
    df_name <- deparse(substitute(df))
    for (col in cols) {
        df[is.na(df[, col]), col] <- val
    }
    if (inplace) {
        assign(df_name, df, envir=.GlobalEnv)
    } else {
        return(df)
    }
}


#' See: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
#' 
#' @export
reset_index <- function(df, index_name='index', drop=FALSE) {
    df <- cbind(index = rownames(df), df)
    rownames(df) <- 1:nrow(df)
    colnames(df)[colnames(df) == "index"] = index_name
    if (drop == TRUE) {
        df <- df[, items_in_a_not_b(colnames(df), 'index')]
    }
    return (df)
}


#' append df2 to df1
#' Avoids bind_rows errors: Can't combine `..1$1` <character> and `..2$1` <double>.
#'
#' @export
append_dataframe <- function(df1, df2, infront=FALSE, reset_index=TRUE) {

    missing_cols <- items_in_a_not_b(colnames(df1), colnames(df2))
    for (col in missing_cols) {
        df2[[col]] <- NA
    }
    if (infront) {
        df <- rbind(df2[, intersect(colnames(df1), colnames(df2))], df1)
    } else {
        df <- rbind(df1, df2[, intersect(colnames(df1), colnames(df2))])
    }
    if (reset_index) {
        df <- reset_index(df, drop=TRUE)
    }

    return(df)
}


#' convert a named list into a row in a dataframe
#'
#' @export
dataframe_row_from_named_list <- function(items) {
    df <- data.frame('1'=unname(items))
    df <- as.data.frame(t(df))  # transpose
    rownames(df) <- 1
    colnames(df) <- names(items)
    return(df)
}


#' convenience function to also set the column name with the transpose
#'
#' @export
stranspose <- function(df, colname=NULL) {
    tdf <- as.data.frame(t(df))
    if (!is.null(colname)) {
        colnames(tdf) <- tdf[colname,]
        tdf <- tdf[items_in_a_not_b(rownames(tdf), colname), ]
    }
    return(tdf)
}
