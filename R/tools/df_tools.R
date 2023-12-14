import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', 'replace_specific_items', .character_only=TRUE)

## Functions
## append_dataframe
## dataframe_row_from_named_list
## fillna
## rename_columns
## reset_index
## stranspose


#' Appends df1 with df2
#' 
#' @description
#' Use this to avoids errors like this one:
#' 
#' Can't combine `..1$1` <character> and `..2$1` <double>.
#' 
#' If you combine dataframes with different types, the data gets coerced to `character`.
#'
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


#' Convert a named list into a row in a dataframe
#'
#' @description
#' Use this to create a custom dataframe with a single row.
#' Helpful for appending a metadata row to an existing dataframe.
#'
#' @examples
#' dataframe_row_from_named_list(c('col1'=1, 'col2'=2, 'col3'=3))
#' 
dataframe_row_from_named_list <- function(items) {
    df <- data.frame('1'=unname(items))
    df <- as.data.frame(t(df))  # transpose
    rownames(df) <- 1
    colnames(df) <- names(items)
    return(df)
}


#' Fill specific column with NA
#' 
#' @description Mirrors Pandas' \href{https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.fillna.html}{fillna}
#' 
#' @param df a dataframe
#' @param cols a list of columns to replace
#' @param val the value to fill with
#' @param inplace TRUE allows you to avoid re-assigning the variable
#' @return Returns a dataframe.
#' 
#' @examples
#' mtcars['new_col'] <- NA
#' head(fillna(mtcars, c('new_col'), 1))
#' 
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


#' Rename dataframe columns
#' 
#' @description This is a thin wrapper around replace_specific_items that acts on dataframe columns
#' 
#' @param df a dataframe
#' @param columns a named list of replacements. uses names to match and values to replace
#' @param inplace TRUE allows you to avoid re-assigning the variable
#' @return Returns a dataframe.
#' 
#' @examples
#' head(rename_columns(mtcars, c('mpg'="MPG", 'disp'="DISP")))
#' 
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


#' Reset index
#' 
#' @description
#' Moves the values in the index to a column. Resets the index to the default integer index.
#' Mirrors Pandas' \href{https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.reset_index.html}{reset_index}.
#' 
#' @param df a dataframe
#' @param index_name select a new name for the index column
#' @param drop if TRUE, does not copy the index values to the new column
#' @return Returns a dataframe with index renamed and new integer index 
#' 
#' @examples
#' reset_index(mtcars, index_name='model')
#' 
#' @references
#' \href{https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column}{StackOverflow post}
#' 
reset_index <- function(df, index_name='index', drop=FALSE) {
    df <- cbind(index = rownames(df), df)
    rownames(df) <- 1:nrow(df)
    colnames(df)[colnames(df) == "index"] = index_name
    if (drop == TRUE) {
        df <- df[, items_in_a_not_b(colnames(df), 'index')]
    }
    return (df)
}


#' Special transpose
#' 
#' @description Allows you to choose which column becomes the new column name
#'
#' @param df dataframe
#' @param colname name of the column that will become the new column name
#' @return Returns the transposed dataframe
#' 
stranspose <- function(df, colname=NULL) {
    tdf <- as.data.frame(t(df))
    if (!is.null(colname)) {
        colnames(tdf) <- tdf[colname,]
        tdf <- tdf[items_in_a_not_b(rownames(tdf), colname), ]
    }
    return(tdf)
}
