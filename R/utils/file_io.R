import::here(readxl, 'read_excel')
import::here(file.path(wd, 'R', 'utils', 'list_tools.R'),
    'filter_list_for_match', .character_only=TRUE)

## Functions
## list_files
## read_csv_from_text
## read_excel_or_csv
## join_many_csv


#' list all files in all subdirectories with a given extension
#' 
#' @export
list_files <- function(dir_path, ext=NULL, recursive = TRUE) {
    all_files = list.files(dir_path, recursive = recursive, full.name=TRUE)

    if (!is.null(ext)) {
        # See: https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
        return (all_files[tools::file_ext(all_files)==ext])
    } else {
        return (all_files)
    }
}


#' Reads and parses csv embedded in .txt files
#'
#' This is useful for data exported from plate readers
#' 
#' @examples
#' # for 96-well plates:
#' df <- read_csv_from_text(
#'   file_path,
#'   skiprows=3, nrows=8,
#'   skipcols=2, ncols=12,
#'   index=LETTERS[1:8],
#'   columns=seq(1, 12)
#' )
#' @export
read_csv_from_text <- function(
    file_path,
    encoding='UTF-16', sep='\t',
    skiprows=0, nrows=NULL,
    skipcols=0, ncols=NULL,
    index=NULL,
    columns=NULL,
    numeric=FALSE
) {

    con = file(file_path, encoding=encoding)
    rawData <- readLines(con)
    close(con)
    
    # autodetermine ranges if not specified
    if(is.null(nrows)) {
        nrows <- length(rawData)
    }
    if(is.null(ncols)) {
        rowArr = unlist(strsplit(rawData[1+skiprows], split='\t'))
        ncols = length(rowArr)-skipcols
    }
    
    # instantiate empty dataframe and append row-by-row
    df <- data.frame(matrix(ncol=ncols, nrow=0))
    for (row in rawData[(1+skiprows):(nrows+skiprows)]) {
        rowArr <- unlist(strsplit(row, split='\t'))
        df[nrow(df) + 1,] = rowArr[(1+skipcols):(ncols+skipcols)]
    }
    
    # rename columns
    colnames(df) <- columns
    rownames(df) <- index

    if(numeric) {
        df[] <- lapply(df, function(x) as.numeric(as.character(x)))
    }

    return(df)
}


#' switch case to read excel or csv based on the extension
#' 
#' @export
read_excel_or_csv <- function(filepath) {
    ext=tools::file_ext(filepath)
    if (ext == 'xlsx') {
        df <- read_excel(filepath)
    } else if (ext == 'csv') {
        df <- read.csv(filepath, header=TRUE, check.names=FALSE)
    } else {
        stop('Please enter a xlsx or csv file.')
    }
    return(df)
}


#' Read all the csv files from a directory and left join them into a single dataframe
#' See: https://stackoverflow.com/questions/5319839/read-multiple-csv-files-into-separate-all_data-frames
#' The paths argument can be a single directory or a list of individual files
#' index_cols=c('gene_id', 'gene_name', 'chromosome')
#' value_cols=c('count')
#' 
#' @export
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
