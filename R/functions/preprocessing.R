import::here(magrittr, '%>%')
import::here(dplyr, 'group_by', 'slice_min', 'ungroup', 'select', 'all_of')
import::here(stringr, 'str_detect')
import::here(plyr, 'mapvalues')
import::here(XML, 'xmlOutputDOM')
import::here(file.path(wd, 'R', 'config', 'replacements.R'),
    'fluorophore_replacements', 'antibody_replacements', 'instr_cfg_colreps', 'ab_inv_colreps',
    .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'append_dataframe', 'rename_columns', 'reset_index', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'text_tools.R'),
    'title_to_snake_case', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'multiple_replacement', 'find_first_match_index', .character_only=TRUE)
import::here(file.path(wd, 'R', 'config', 'flow.R'),
    'id_cols', 'initial_gates', 'cell_type_spell_check', 'cell_type_ignore',
    .character_only=TRUE)

## Functions
## create_survival_table
## preprocess_flowjo_export
## preprocess_antibody_inventory
## preprocess_instrument_config
## sort_groups_by_metric
## spillover_to_xml
## parse_flowjo_metadata


#' Create Survival Table
#' 
#' Reshape for survfit
#' Use mouse_id as the id column
#' 
create_survival_table <- function(start_info, end_info, end_date) {

    df <- merge(start_info, end_info,
       by=c("mouse_id", "group"),
       all.x=TRUE, all.y=TRUE, suffixes=c('_start', '_end'))

    # calculate day from date
    df[['day']] <- as.numeric(difftime(
        as.Date(df[['date_end']], format = '%m/%d/%y'),
        as.Date(df[['date_start']], format = '%m/%d/%y'),
        units = paste0('day', 's')
    ))
    df[['week']] <- df[['day']]/7
    df[['is_dead']] <- as.numeric( df[['date_end']] < end_date )

    # extend line to experiment start for survfit
    dummy_start <- data.frame(
        group = unique(df[["group"]]),
        date_start = min(df[['date_start']]),
        date_end = min(df[['date_start']]),
        day = 0,
        week = 0,
        is_dead = 0
    )

    dummy_end <- df %>%
        group_by(.data[['group']]) %>%
        slice_min(order_by = date_end, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        select(all_of(c('group', 'date_end', 'day', 'week')))
    dummy_end[['date_start']] = min(df[['date_start']])
    dummy_end[['day']] = dummy_end[['day']] - 0.001
    dummy_end[['week']] = dummy_end[['week']] - 0.001
    dummy_end[['is_dead']] = 0

    df <- append_dataframe(df, dummy_start)
    df <- append_dataframe(df, dummy_end)

    return(df)
}


#' Preprocess Each Flowjo Export
#' 
preprocess_flowjo_export <- function(raw_table,
    metric_name='num_cells',  # num_cells, gmfi, or rsdev
    include_initial_gates=FALSE  # turn on for analyze_counts
) {

    raw_table <- rename_columns(raw_table,
        c('X1'='fcs_name',
          'Count'='Ungated',
          'Mean (Comp-Alexa Fluor 488-A)'='Ungated',  # find a better solution
          'Mode (Comp-Alexa Fluor 488-A)'='Ungated',
          'Geometric Mean (Comp-Alexa Fluor 488-A)'='Ungated'
        )
    )
    colnames(raw_table) <- unname(sapply(colnames(raw_table), function(x) strsplit(x, ' \\| ')[[1]][1]))
    raw_table <- raw_table[!(raw_table[['fcs_name']] %in% c('Mean', 'SD')), ]  # drop summary statistics
    raw_table <- raw_table[!str_detect(raw_table[['fcs_name']], 'unstained'), ]  # drop unstained cells
    raw_table <- reset_index(raw_table, drop=TRUE)

    # Reshape so each row is a gate in each sample
    df <- pivot_longer(raw_table,
        names_to = "gate", values_to = metric_name,
        cols=items_in_a_not_b(colnames(raw_table), c(id_cols, 'Count', if (include_initial_gates) {initial_gates} else {NULL} )),
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


#' Preprocess Antibody Inventory
#'
preprocess_antibody_inventory <- function(df) {

    df <- df[, 1:(find_first_match_index('\\.{3}\\d{2}', colnames(df))-1)]  # filter extra columns
    colnames(df) <- unlist(lapply(colnames(df), title_to_snake_case))  # column names
    colnames(df) <- unlist(lapply(colnames(df), function(x) gsub('[.]', '', x)))  # column names
    df <- rename_columns(df, ab_inv_colreps)
    
    # remove 'c2 a0', the "no-break space"
    # see: https://stackoverflow.com/questions/68982813/r-string-encoding-best-practice
    for (col in colnames(df)) {
        df[(!is.na(df[[col]]) & (df[[col]] == enc2utf8("\u00a0"))), col] <- NA
        df[[col]] = unlist(
            lapply(df[[col]], function(x) gsub(paste0('.*^', enc2utf8("\u00a0"), '+'), '', x))
        )    
    }

    # standardize names
    if ('alternative_name' %in% colnames(df)) {
        df[['alternative_name']] = multiple_replacement(
            df[['alternative_name']], antibody_replacements
        )
    }
    df[['antibody']] = multiple_replacement(df[['antibody']], antibody_replacements)
    df[['fluorophore']] = multiple_replacement(df[['fluorophore']], fluorophore_replacements)
    
    # fill in missing fixable_dyes
    fixable_dyes <- c("DAPI", "Zombie UV", "Zombie Aqua")
    mask <- (df[['antibody']] %in% fixable_dyes) & (is.na(df[['fluorophore']]))
    df[mask, 'fluorophore'] <- df[mask, 'antibody']

    return(df)
}


#' Preprocess_instrument_config
#' 
preprocess_instrument_config <- function(df) {

    # format and rename columns
    colnames(df) <- unlist(lapply(colnames(df), title_to_snake_case))
    df <- rename_columns(df, instr_cfg_colreps)

    # standardize fluorophore names
    df[['fluorophore']] = multiple_replacement(
        df[['fluorophore']], fluorophore_replacements
    )

    return(df)

    # order list
    df <- df[order(
        df[['bandpass_filter']], df[['excitation']],
        decreasing = c(FALSE, TRUE),
        method = "radix"), ]
}


#' Sort Groups by Metric
#' 
#' Sort rows by groups in decreasing order of the median percent cells
#' 
sort_groups_by_metric <- function(
    df,
    x='cell_type',  # group
    y='pct_cells',  # metric
    groups=c('groupby')  # enter a collection of groups
) {

    # find the median of each group
    x_axis_order <- df %>%
        group_by(.data[[x]]) %>%
        summarize(median = median(.data[[y]], na.rm=TRUE)) %>%
        arrange(-.data[['median']]) %>%
        select(.data[[x]])

    # sort each group by decreasing median
    sorted <- df %>% arrange(
        match(.data[[x]], unlist(x_axis_order)),
        !!! rlang::syms(groups)
    )

    return(sorted)
}


#' Constructs a mtx file for the compensation matrix for Flowjo import 
#' Make sure to include the diagonal 1 values
#' 
spillover_to_xml <- function(spillover_table, channels, name='test') {

    xml <- xmlOutputDOM("gating:gatingML")

    # <transforms:spilloverMatrix>
    xml$addTag("transforms:spilloverMatrix",
        attrs=list(spectral="0",
            weightOptAlgorithmType="OLS",
            prefix="Comp-",
            name=name,
            editable="1",
            color="#00ccff",
            version="FlowJo-10.10.0",
            status="FINALIZED",
            'transforms:id'="",
            suffix=""
        ),
        close=FALSE
    )

        # <data-type:parameters>
        xml$addTag("data-type:parameters", close=FALSE)
        for (channel in channels) {
            channel_name <- strsplit(channel, ' :: ')[[1]][1]
            xml$addTag("data-type:parameter",
                attrs = list(
                    'data-type:name'=channel_name,
                    userProvidedCompInfix=paste0('Comp-', channel_name))
            )
        }
        xml$closeTag()

        # build table
        for (channel in channels) {

            # subset by channel
            spillover_subtable <- spillover_table[
                (spillover_table[['- % Fluorochrome']]==channel),
                c('Fluorochrome', 'Spectral Overlap')]
            spillover_subtable[['Fluorochrome']] <- unname(sapply(
                spillover_subtable[['Fluorochrome']],
                function(x) strsplit(x, ' :: ')[[1]][1]
            ))
            spillover_subtable <- rename_columns(spillover_subtable, c(
                'Fluorochrome'='data-type:parameter',
                'Spectral Overlap'='transforms:value'
            ))

            # <transforms:spillover >
            channel_name <- strsplit(channel, ' :: ')[[1]][1]
            xml$addTag(
                "transforms:spillover",
                attrs = list(
                    'data-type:parameter'=channel_name,
                    userProvidedCompInfix=paste0('Comp-', channel_name)),
                close=FALSE
            )

            # <transforms:coefficient >
            for (row in row.names(spillover_subtable)) {
                xml$addTag("transforms:coefficient", attrs = spillover_subtable[row,])
            }

            xml$closeTag()
        }

    xml$closeTag()
    return(xml)
}


# ----------------------------------------------------------------------
# To be deprecated


#' Parse metadata saved in the fcs name
#' 
#' @description This function is specific to the naming convention of the experiment. To be deprecated.
#'
parse_flowjo_metadata <- function(
    df,
    cols=c('organ', 'mouse_id', 'treatment_group', 'strain')
) {

    df[['metadata']] <- as.character(
        lapply(strsplit(df[['fcs_name']], '_'),
        function(x) paste( x[4:length(x)-1], collapse='_' )) 
    )

    # single value column with only 'F' causes R to interpret this as FALSE
    if ('sex' %in% colnames(df)) {
        if ( typeof(df[['sex']])=='logical' ) {
            df[['sex']] <- mapvalues(
                tmp[['sex']], from = c(FALSE, TRUE), to = c('F', 'T')
            )
        }        
    }

    for (i in 1:length(cols)) {
        col <- cols[[i]]
        tryCatch({
            df[[col]] <- as.character(
                lapply(strsplit(df[['metadata']], '-'), function(x) x[[i]])
            )
        },
        error = function(condition) {
            log_print(paste("Column", col, "not in metadata"))
        })
    }

    df[['strain']] <- unlist(as.character(
        lapply(df[['mouse_id']], function(x) strsplit(x, '_')[[1]][[1]])
    ))

    return(df)
}
