import::here(wrapr, 'orderv')
import::here(file.path(wd, 'R', 'config', 'replacements.R'),
    'fluorophore_replacements', 'antibody_replacements', 'instr_cfg_colreps', 'ab_inv_colreps',
    .character_only=TRUE)
import::here(file.path(wd, 'R', 'utils', 'df_tools.R'),
    'rename_columns', .character_only=TRUE)
import::here(file.path(wd, 'R', 'utils', 'text_tools.R'),
    'title_to_snake_case', .character_only=TRUE)
import::here(file.path(wd, 'R', 'utils', 'list_tools.R'),
    'multiple_replacement', 'find_first_match_index', .character_only=TRUE)

## Functions
## preprocess_instrument_config
## preprocess_antibody_inventory


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
