import::here(plyr, 'rbind.fill')
import::here(readxl, 'read_excel')
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'filter_list_for_match', .character_only=TRUE)


## Functions
## append_many_csv
## list_files
## join_many_csv
## read_excel_or_csv


#' Aggregate csv files by appending them rowwise
#' 
#' @description
#' Read all the csv files from a directory and append them into a single dataframe
#' 
#' @export
append_many_csv <- function(dir_path, sep=',', row_names=NULL, return_list=FALSE) {
    filenames <- list.files(dir_path, full.names=TRUE)

    dfs <- new.env()
    for (file in filenames) {
        df <- read.csv(
            file, sep=sep, row.names=row_names,
            header=TRUE, check.names=FALSE
        )

        # remove nulls
        df <- Filter(function(x) !all(is.na(x)), df)  # columns
        df <- df[rowSums(is.na(df)) != ncol(df), ]  # rows

        # repair colnames
        missing_colnames <- which(nchar(colnames(df))==0)
        num_missing <- length(missing_colnames)
        if (num_missing > 0) {
            colnames(df)[which(nchar(colnames(df))==0)] <- paste0("X", 1:num_missing)
        }

        # add filename
        cols <- colnames(df)
        df[['filename']] <- basename(file)
        df <- df[, c('filename', cols)]

        dfs[[file]] <- df
    }
    
    dfs <- as.list(dfs)

    if (return_list) {
        return(dfs)
    } else {
        combined <- do.call(rbind.fill, dfs)
        return(combined)  
    }
}


#' List all files with a specific extension
#' 
#' @description
#' This is a thin wrapper around [list.files()].
#' 
#' @references
#' \href{https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching}{StackOverflow post}
#' 
list_files <- function(dir_path, ext=NULL, recursive = TRUE) {
    all_files = list.files(dir_path, recursive = recursive, full.name=TRUE)

    if (!is.null(ext)) {
        # See: https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
        return (all_files[tools::file_ext(all_files)==ext])
    } else {
        return (all_files)
    }
}


#' Aggregate csv files by joining them column-wise
#' 
#' @description
#' Read all the csv or tsv files from a directory and left join them into a single dataframe
#' 
#' @references
#' \href{https://stackoverflow.com/questions/5319839/read-multiple-csv-files-into-separate-all_reads-frames}{StackOverflow post}
#' 
join_many_csv <- function(
    paths, index_cols, value_cols=NULL,
    all.x=FALSE, all.y=FALSE, recursive=TRUE
) {

    # distinguish if paths is a directory or a list of files
    if (length(paths)==1) {
        if (dir.exists(paths)) {
             paths <- list_files(paths, recursive=recursive)
        }
    } else if (length(paths[file.exists(paths)]) == 0) {
        stop(paste("no files found!"))
    }

    # split into csv and tsv files
    csv_paths = filter_list_for_match(paths, 'csv')
    csv_list <- lapply(csv_paths, read.csv, sep=',')
    tsv_paths = filter_list_for_match(paths, 'tsv')
    tsv_list <- lapply(tsv_paths, read.csv, sep='\t') 
    df_list = c(csv_list, tsv_list)

    filenames = c(tools::file_path_sans_ext(basename(paths)))
    # Warning: column names ‘count.x’, ‘count.y’ are duplicated in the result
    # See: https://stackoverflow.com/questions/38603668/suppress-any-emission-of-a-particular-warning-message
    withCallingHandlers({
        all_data <- Reduce(
            function(...) merge(..., by=index_cols, all.x=all.x, all.y=all.y),
            lapply(df_list, "[", c(index_cols, value_cols))
        )
    }, warning = function(w) {
        # print(conditionMessage(w))
        if (startsWith(conditionMessage(w), "column names")) {
            invokeRestart( "muffleWarning" )
        }
    })
    
    # rename columns
    colnames(all_data) = c(
        index_cols,  # index_cols
        as.list(outer(value_cols, filenames, paste, sep='-'))  # suffix value_cols with filename
    )
    return(all_data)
}


#' switch case to read excel or csv based on the extension
#' 
#' @export
read_excel_or_csv <- function(filepath) {
    ext=tools::file_ext(filepath)
    if (ext == 'xlsx') {
        df <- read_excel(filepath, .name_repair = "unique_quiet")
    } else if (ext == 'csv') {
        df <- read.csv(filepath, header=TRUE, check.names=FALSE)
    } else {
        stop('Please enter a xlsx or csv file.')
    }
    return(df)
}
