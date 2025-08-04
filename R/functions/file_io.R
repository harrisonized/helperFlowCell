import::here(tidyr, 'pivot_longer')
import::here(plyr, 'mapvalues')
import::here(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'preprocess_flowjo_export', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'reset_index', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', .character_only=TRUE)

## Functions
## import_flowjo_export
## import_flow_metadata


#' Import Flowjo Export
#' 
#' Read data exported from flowjo, use this to import counts, gmfi, or sdev
#' If metric=='count', the import will include the Single Cells and Live Cells gates
#'
import_flowjo_export <- function(
    dirpath,
    metric_name='num_cells',
    include_initial_gates=FALSE
) {
    raw_tables <- append_many_csv(dirpath, recursive=TRUE,
        na_strings=c('n/a'), include_filepath=FALSE, return_list=TRUE)
    if (length(raw_tables)==0) {
        return(NULL)
    }

    dfs <- lapply(raw_tables, function(x) preprocess_flowjo_export(
        x, metric_name=metric_name, include_initial_gates=include_initial_gates
    ))
    df <- do.call(rbind, dfs)
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
        if ( typeof(flow_metadata[['sex']])=='logical' ) {
            flow_metadata[['sex']] <- mapvalues(
                flow_metadata[['sex']],
                from = c(FALSE, TRUE),
                to = c('F', 'T'),
                warn_missing=FALSE
            )
        }
    }

    return(flow_metadata)
}
