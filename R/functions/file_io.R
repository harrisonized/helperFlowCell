import::here(tidyr, 'pivot_longer')
import::here(plyr, 'mapvalues')
import::here(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'preprocess_flowjo_export', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', .character_only=TRUE)
import::here(file.path(wd, 'R', 'config', 'flow.R'),
    'id_cols', 'initial_gates', 'cell_type_spell_check', 'cell_type_ignore',
    .character_only=TRUE)

## Functions
## import_flowjo_export
## import_flow_metadata


#' Import Flowjo Export
#' 
#' Read data exported from flowjo, use this to import counts, gmfi, or sdev
#' TODO: change this to pivot before append
#'
import_flowjo_export <- function(
    dirpath,
    metric_name='num_cells',  # num_cells, gmfi, or rsdev
    import_error=TRUE
) {

    raw_table <- append_many_csv(dirpath, recursive=TRUE, na_strings=c('n/a'))
    if (import_error) {
        if (is.null(raw_table)) {
                msg <- paste("No data found. Please check", dirpath, '...')
                stop(msg)
            }
    }
    raw_table <- preprocess_flowjo_export(raw_table)

    # Reshape so each row is a gate in each sample
    df <- pivot_longer(raw_table,
        names_to = "gate", values_to = metric_name,
        cols=items_in_a_not_b(colnames(raw_table), c(id_cols, 'Count', initial_gates)),
        values_drop_na = TRUE
    )
    df[['cell_type']] <- unlist(lapply(strsplit(df[['gate']], '/'), function(x) x[length(x)]))
    df[['cell_type']] <- multiple_replacement(df[['cell_type']], cell_type_spell_check)
    
    # filter ignored cell types
    for (cell_type in c(cell_type_ignore, 'mNeonGreen+')) {
        df <- df[!str_detect(df[['cell_type']], cell_type), ]
    }

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
