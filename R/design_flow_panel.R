## Build the design matrix

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
library('readxl')
library('tidyr')
suppressMessages(library('dplyr'))
suppressMessages(library('tibble'))
library('optparse')
library('logr')
source(file.path(wd, 'R', 'preprocessing.R'))
source(file.path(wd, 'R', 'functions', 'file_io.R'))  # read_excel_or_csv
source(file.path(wd, 'R', 'functions', 'df_tools.R'))  # rename_columns
source(file.path(wd, 'R', 'functions', 'list_tools.R'))  # items_in_a_not_b


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-file"), default='data/panel.txt',
                metavar='data/panel.txt', type="character",
                help="specify the antibodies to use in your flow panel"),

    make_option(c("-a", "--antibody-inventory"), default='ref/antibody_inventory.xlsx',
                metavar='ref/antibody_inventory.xlsx', type="character",
                help="antibody inventory table"),

    make_option(c("-c", "--instrument-config"), default='ref/instrument_config.xlsx',
                metavar='ref/instrument_config.xlsx', type="character",
                help="instrument configuration file"),

    make_option(c("-o", "--output-dir"), default="data/output",
                metavar="data/output", type="character",
                help="set the output directory for the data"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]

output_dir <- file.path(wd, opt[['output-dir']])
if (!troubleshooting) {
    troubleshooting_dir = file.path(output_dir, 'troubleshooting')
    if (!dir.exists(file.path(troubleshooting_dir))) {
        dir.create(file.path(troubleshooting_dir), recursive=TRUE)
    }
}

# Start Log
start_time = Sys.time()
log <- log_open(paste0("design_flow_panel-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Instrument Config

instr_cfg <- read_excel_or_csv(file.path(wd, opt[['instrument-config']]))
instr_cfg <- preprocess_instrument_config(instr_cfg)
instr_cfg_long <- separate_rows(instr_cfg, 'fluorophore', sep=', ')

# save
if (!troubleshooting) {
    filepath = file.path(troubleshooting_dir, 
        paste0('_', tools::file_path_sans_ext(basename(opt[['instrument-config']])), '.csv')
    )
    write.table(instr_cfg, file = filepath, row.names = FALSE, sep=',')
}


# ----------------------------------------------------------------------
# Antibody Inventory

ab_inv <- read_excel_or_csv(file.path(wd, opt[['antibody-inventory']]))
ab_inv <- preprocess_antibody_inventory(ab_inv)

# add fluorophores in this list to your instrument-config file
unavailable_fluorophores = sort(items_in_a_not_b(
    unique( ab_inv[['fluorophore']] ), unique( instr_cfg_long[['fluorophore']] )
))

# save
if (!troubleshooting) {
    filepath = file.path(troubleshooting_dir, 'antibodies.txt')
    write.table(sort(unique(ab_inv[['antibody']])), filepath,
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    filepath = file.path(troubleshooting_dir, 'unavailable_fluorophores.txt')
    write.table(unavailable_fluorophores, filepath,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
}


# ----------------------------------------------------------------------
# Build design_matrix

# Read data
panel_list <- read_csv_from_text(file.path(wd, opt[['input-file']]),
    columns='antibody', encoding='UTF-8')
ab_shortlist <- merge(panel_list, ab_inv, by='antibody',
    all.x=FALSE, all.y=FALSE, na_matches = 'never', sort=FALSE)
ab_shortlist <- ab_shortlist[, colnames(ab_inv)]

antibodies_not_found = sort(items_in_a_not_b(
    unique( panel_list[['antibody']] ), unique( ab_inv[['antibody']] )
))

# save
if (!troubleshooting) {
    filepath = file.path(output_dir, 
        paste0(tools::file_path_sans_ext(basename(opt[['antibody-inventory']])),
               '-shortlist', '.csv')
    )
    write.table(ab_shortlist , file = filepath, row.names = FALSE, sep=',')

    if (length(antibodies_not_found) > 0) {
        filepath = file.path(troubleshooting_dir, 'antibodies_not_found.txt')
        write.table(antibodies_not_found, filepath,
            row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
}

# compute base table with antibody + instrument channel information in each row
ab_counts <- ab_shortlist[, c('antibody', 'company', 'catalog_no', 'fluorophore')] %>% 
    merge(instr_cfg_long[, c('fluorophore', 'channel_id')],
          by='fluorophore', suffixes=c('', '_'),
          all.x=TRUE, all.y=FALSE, na_matches = 'never') %>% 
    fillna(c('channel_id'), val='other') %>%
    # group_by(antibody, channel_id) %>%
    # group_by(.data[['antibody']], .data[['channel_id']]) %>%
    group_by(!!!syms(c('antibody', 'channel_id'))) %>%
    summarise(
        num_fluorophores=n_distinct(fluorophore),
        fluorophores = toString(unique(fluorophore)),
        most_common_fluorophore=names(table(fluorophore))[which.max(table(fluorophore))],
        .groups='keep'
    )

# pivot to table where rows are antibodies and columns are instrument channels
design_matrix <- pivot_wider(
        ab_counts,
        id_cols = 'antibody', names_from = 'channel_id',
        values_from = 'num_fluorophores', values_fill = 0
    ) %>%
    ungroup()


# ----------------------------------------------------------------------
# Sort

# Note: rowSums should be computed before colSums, so that the
# 'num_channels_per_ab' sums to the number of antibodies

# Sort rows
num_channels_per_ab <- select_if(design_matrix, is.numeric) %>%
    mutate_at(colnames(select_if(design_matrix, is.numeric)),
              function(x) ifelse(x >= 1, 1, 0)) %>%
    rowSums()
design_matrix[['num_channels_per_ab']] <- num_channels_per_ab
design_matrix <- design_matrix[order(design_matrix[['num_channels_per_ab']]), ]

# Sort columns
num_abs_per_channel <- select_if(design_matrix, is.numeric) %>%
    mutate_at(colnames(select_if(design_matrix, is.numeric)),
              function(x) ifelse(x >= 1, 1, 0)) %>%
    colSums()
design_matrix <- append_dataframe(
    design_matrix, dataframe_row_from_named_list(num_abs_per_channel),
    reset_index=FALSE
)
design_matrix <- design_matrix[ ,
    order(unlist(design_matrix[nrow(design_matrix), ]))
]

# Move channel_cols to the middle
id_cols <- c('antibody', 'num_channels_per_ab')
channel_cols <- items_in_a_not_b(colnames(design_matrix), c(id_cols, 'other'))
design_matrix <- design_matrix[ , c(id_cols, channel_cols, 'other')]


# ----------------------------------------------------------------------
# Add metadata

# Add comma-separated fluorophore list to each row
row_metadata <- ab_counts %>%
    group_by(antibody) %>%
    summarise(fluorophores = toString(unique(most_common_fluorophore)))
design_matrix = merge(design_matrix, row_metadata,
    by='antibody', suffixes=c('', '_'),
    all.x=TRUE, all.y=FALSE, na_matches = 'never', sort=FALSE
)

# Add available fluorophores and instrument info to each column
col_metadata <- ab_counts %>%
    group_by(channel_id) %>%
    summarise(all_fluorophores = toString(unique(most_common_fluorophore))) %>%
    merge(instr_cfg[, c('channel_id', 'laser', 'pmt')],
          by='channel_id', suffixes=c('', '_'),
          all.x=TRUE, all.y=FALSE, na_matches = 'never', sort=FALSE)
design_matrix <- design_matrix %>% append_dataframe(
    stranspose(col_metadata, colname='channel_id'),
    infront=TRUE, reset_index=FALSE
)


# ----------------------------------------------------------------------
# Postprocessing

# Rename columns from channel_id (1, 2, 3, etc.) to representative_fluorophore (AF488, FITC, PE, etc.)
fluorophore_for_channel <- ab_counts %>%
    group_by(channel_id) %>%
    summarise(representative_fluorophore=names(
        table(most_common_fluorophore)[which.max(table(most_common_fluorophore))]
    )) %>%
    deframe()
design_matrix <- rename_columns(design_matrix, fluorophore_for_channel)

# Rename last row
row.names(design_matrix)[nrow(design_matrix)] <- 'count'
design_matrix <- reset_index(design_matrix)  # export index information without frameshift

# save
if (!troubleshooting) {
    filepath = file.path(output_dir, 'design_matrix.csv')
    write.table(design_matrix, file = filepath, row.names = FALSE, sep=',')
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
