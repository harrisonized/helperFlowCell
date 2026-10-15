import::here(tidyr, 'pivot_longer')
import::here(plyr, 'mapvalues')
import::here(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'preprocess_flowjo_export', 'preprocess_flowjo_export_v2', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'reset_index', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', .character_only=TRUE)

## Functions
## import_flowjo_export
## import_flowjo_export_v2
## import_flow_metadata


#' Import Flowjo Export
#' 
#' Read data exported from flowjo, use this to import counts, gmfi, or sdev
#' If metric=='count', the import will include the Single Cells and Live Cells gates
#'
import_flowjo_export <- function(
    dirpath,
    metric_name='num_cells',
    last_initial_gate=NA
) {
    raw_tables <- append_many_csv(dirpath, recursive=TRUE,
        na_strings=c('n/a'), include_filepath=FALSE, return_list=TRUE)
    if (length(raw_tables)==0) {
        return(NULL)
    }

    dfs <- lapply(raw_tables, function(x) preprocess_flowjo_export(
        x, metric_name=metric_name, last_initial_gate=last_initial_gate
    ))
    df <- do.call(rbind, dfs)
    df <- reset_index(df, drop=TRUE)

    return(df)
}


#' Import Flowjo Export
#' 
#' Read data exported from flowjo, use this to import counts, gmfi, or sdev
#' If metric=='count', the import will include the Single Cells and Live Cells gates
#' 
#' This is part of a refactoring effort to generalize data import so that files with
#' different initial gates can be merged.
#'
import_flowjo_export_v2 <- function(
    dirpath,
    metric_name='num_cells',
    divisor="Live Cells",
    recursive=FALSE,
    na_strings=c('n/a')
) {

    filepaths <- list.files(dirpath, recursive=recursive, full.names=TRUE)
    if (length(filepaths)==0) {
        return(NULL)
    }

    dfs <- new.env()
    for (filepath in filepaths) {
        raw_table <- read.csv(
            filepath, sep=',', row.names=NULL,
            header=TRUE, check.names=FALSE,
            na.strings=na_strings
        )

        # remove entirely null columns and rows
        raw_table <- Filter(function(x) !all(is.na(x)), raw_table)  # columns
        raw_table <- raw_table[rowSums(is.na(raw_table)) != ncol(raw_table), ]  # rows

        # repair column names
        missing_colnames <- which(nchar(colnames(raw_table))==0)
        num_missing <- length(missing_colnames)
        if (num_missing > 0) {
            colnames(raw_table)[which(nchar(colnames(raw_table))==0)] <- paste0("X", 1:num_missing)
        }

        df <- preprocess_flowjo_export_v2(raw_table, metric_name=metric_name, divisor=divisor)

        dfs[[filepath]] <- df
    }

    dfs <- as.list(dfs)
    if (length(dfs) > 1) {
        df <- do.call(rbind, dfs[filepaths])
    }
    df <- reset_index(df, drop=TRUE)

    return(df)
}


#' Import Flow Metadata
#' 
import_flow_metadata <- function(dirpath) {

    flow_metadata <- append_many_csv(
        dirpath, recursive=TRUE, include_filepath=FALSE)

    # note: this is required
    if (is.null(flow_metadata)) {
        msg <- paste("No metadata found. Please check", dirpath, '...')
        stop(msg)
    }

    # single value column with only 'F' causes R to interpret this as FALSE
    if ('sex' %in% colnames(flow_metadata)) {
        flow_metadata[['sex']] <- mapvalues(
            flow_metadata[['sex']],
            from = c(FALSE, TRUE),
            to = c('F', 'T'),
            warn_missing=FALSE
        )
    }

    for (col in c('viable_cells_conc', 'total_vol', 'total_viable_cells')) {
        if (col %in% colnames(flow_metadata)) {
            flow_metadata[[col]] <- as.numeric(flow_metadata[[col]])
        }
    }

    return(flow_metadata)
}
