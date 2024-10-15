import::here(stringr, 'str_detect')
import::here(wrapr, 'orderv')
import::here(XML, 'xmlOutputDOM')
import::here(file.path(wd, 'R', 'config', 'replacements.R'),
    'fluorophore_replacements', 'antibody_replacements', 'instr_cfg_colreps', 'ab_inv_colreps',
    .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'rename_columns', 'reset_index', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'text_tools.R'),
    'title_to_snake_case', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'multiple_replacement', 'find_first_match_index', .character_only=TRUE)

## Functions
## compensation_csv_to_mtx
## parse_flowjo_metadata
## preprocess_flowjo_export
## preprocess_antibody_inventory
## preprocess_instrument_config


#' Constructs a mtx file for the compensation matrix for Flowjo import
#' 
compensation_csv_to_mtx <- function(cytometer_settings, channels, name='test') {

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

            channel_name <- strsplit(channel, ' :: ')[[1]][1]
            spillover_table <- cytometer_settings[
                (cytometer_settings[['- % Fluorochrome']]==channel),
                c('Fluorochrome', 'Spectral Overlap')]
            spillover_table[['Fluorochrome']] <- unname(sapply(
                spillover_table[['Fluorochrome']],
                function(x) strsplit(x, ' :: ')[[1]][1]
            ))
            spillover_table <- rename_columns(spillover_table, c(
                'Fluorochrome'='data-type:parameter',
                'Spectral Overlap'='transforms:value'
            ))

            # <transforms:spillover >
            xml$addTag(
                "transforms:spillover",
                attrs = list(
                    'data-type:parameter'=channel_name,
                    userProvidedCompInfix=paste0('Comp-', channel_name)),
                close=FALSE
            )

            # <transforms:coefficient >
            for (row in row.names(spillover_table)) {
                xml$addTag("transforms:coefficient", attrs = spillover_table[row,])
            }

            xml$closeTag()
        }

    xml$closeTag()
    return(xml)
}


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


#' Preprocess Flowjo Export
#' 
preprocess_flowjo_export <- function(df) {

    df <- rename_columns(df, c('X1'='fcs_name'))
    colnames(df) <- gsub('[[:space:]]\\|[[:space:]]Count', '', colnames(df))  # remove suffix
    df <- df[!(df[['fcs_name']] %in% c('Mean', 'SD')), ]  # drop summary statistics
    df <- df[!str_detect(df[['fcs_name']], 'unstained'), ]  # drop unstained cells
    df <- reset_index(df, drop=TRUE)
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

    # order list
    df <- df[orderv(df[, c('excitation', 'bandpass_filter')], decreasing=TRUE), ]
}
