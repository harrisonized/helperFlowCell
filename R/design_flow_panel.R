## Assist in designing flow panels

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
library('readxl')
library('tidyr')
suppressMessages(library('dplyr'))
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

    make_option(c("-a", "--antibody-inventory"), default='data/antibody_inventory.xlsx',
                metavar='data/antibody_inventory.xlsx', type="character",
                help="antibody inventory table"),

    make_option(c("-c", "--instrument-config"), default='data/instrument_config.xlsx',
                metavar='data/instrument_config.xlsx', type="character",
                help="instrument configuration file"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]

if (!troubleshooting) {
    # note: except for the last file, all files are output in the troubleshooting folder
    data_dir = file.path(wd, dirname(opt[['input-file']]), 'troubleshooting')
    if (!dir.exists(file.path(data_dir, 'lists'))) {
        dir.create(file.path(data_dir, 'lists'), recursive=TRUE)
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
available_fluorophores <- instr_cfg_long[['fluorophore']]

# save
if (!troubleshooting) {
    filepath = file.path(data_dir, '_instrument_config.csv')
    write.table(instr_cfg, file = filepath, row.names = FALSE, sep=',')

    filepath = file.path(data_dir, 'lists', 'available_fluorophores.txt')
    write.table(sort(unique(available_fluorophores)), filepath,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
}


# ----------------------------------------------------------------------
# Antibody Inventory

ab_inv <- read_excel_or_csv(file.path(wd, opt[['antibody-inventory']]))
ab_inv <- preprocess_antibody_inventory(ab_inv)

unavailable_fluorophores = sort(items_in_a_not_b(
    unique( ab_inv[['fluorophore']] ),
    unique( available_fluorophores )
))

# save
if (!troubleshooting) {

    filepath = file.path(data_dir, '_antibody_inventory.csv')
    write.table(ab_inv, file = filepath, row.names = FALSE, sep=',')

    filepath = file.path(data_dir, 'lists', 'antibodies.txt')
    write.table(sort(unique(ab_inv[['antibody']])), filepath,
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    filepath = file.path(data_dir, 'lists', 'alternative_names.txt')
    write.table(sort(unique(ab_inv[['alternative_name']])), filepath,
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    filepath = file.path(data_dir, 'lists', 'all_fluorophores.txt')
    write.table(sort(unique(ab_inv[['fluorophore']])), filepath,
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    filepath = file.path(data_dir, 'lists', 'unavailable_fluorophores.txt')
    write.table(sort(unique(unavailable_fluorophores)), filepath,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
}


# ----------------------------------------------------------------------
# Build table

# Read data
panel_list <- read_csv_from_text(file.path(wd, opt[['input-file']]), encoding='UTF-8')
panel <- ab_inv[
        (ab_inv[['antibody']] %in% panel_list[[1]]),
        c('antibody', 'company', 'catalog_no', 'fluorophore')
    ] %>% 
    merge(instr_cfg_long[, c('fluorophore', 'channel_id')],
          by='fluorophore', suffixes=c('', '_'),
          all.x=TRUE, all.y=FALSE, na_matches = 'never') %>% 
    fillna(c('channel_id'), val='other')  # Biotin

# compute metadata for the panel
ab_counts <- panel %>%
    group_by(antibody, channel_id) %>%
    summarise(
        num_fluorophores=n_distinct(fluorophore),
        fluorophores = toString(unique(fluorophore)),
        most_common_fluorophore=names(table(fluorophore))[which.max(table(fluorophore))],
        .groups='keep'
    )

# pivot to table where rows are antibodies and columns are instrument channels
df <- pivot_wider(
        ab_counts,
        id_cols = 'antibody', names_from = 'channel_id',
        values_from = 'num_fluorophores', values_fill = 0
    ) %>% ungroup()


# ----------------------------------------------------------------------
# Sort

# Note: rowSums should be computed before colSums, so that the
# 'num_channels_per_ab' sums to the number of antibodies

# Sort rows
num_channels_per_ab <- df %>%
    select_if(is.numeric) %>%
    mutate_at(colnames(select_if(df, is.numeric)),
              function(x) ifelse(x >= 1, 1, 0)) %>%
    rowSums()
df[['num_channels_per_ab']] <- num_channels_per_ab
df <- df[order(df[['num_channels_per_ab']]), ]

# Sort columns
num_abs_per_channel <- df %>%
    select_if(is.numeric) %>%
    mutate_at(colnames(select_if(df, is.numeric)),
              function(x) ifelse(x >= 1, 1, 0)) %>%
    colSums()
df <- append_dataframe(df, dataframe_row_from_named_list(num_abs_per_channel),
                       reset_index=FALSE)
df <- df[ , order(unlist(df['count', ]), decreasing=TRUE)]  # sort columns

# Move channel_cols to the middle
non_channel_cols <- c('antibody', 'num_channels_per_ab', 'other')
channel_cols <- items_in_a_not_b(colnames(df), non_channel_cols)
df <- df[ , c('antibody', 'num_channels_per_ab', channel_cols, 'other')]


# ----------------------------------------------------------------------
# Add metadata

# for each antibody, aggregate fluorophores into comma-separated list
abs_vs_fluorophores <- ab_counts %>%
    group_by(antibody) %>%
    summarise(fluorophores = toString(unique(most_common_fluorophore)))

# add comma-separated fluorophore list to each row
df = merge(df, abs_vs_fluorophores,
           by='antibody', suffixes=c('', '_'),
           all.x=TRUE, all.y=FALSE, na_matches = 'never', sort=FALSE)

# compute instrument metadata
instr_info <- ab_counts %>%
    group_by(channel_id) %>%
    summarise(fluorophores = toString(unique(most_common_fluorophore)),
              representative_fluorophore=names(
                  table(most_common_fluorophore)[which.max(table(most_common_fluorophore))]
              )) %>%
    merge(instr_cfg[, c('channel_id', 'laser', 'pmt')],
          by='channel_id', suffixes=c('', '_'),
          all.x=TRUE, all.y=FALSE, na_matches = 'never', sort=FALSE)

# add available fluorophores and instrument info to each column
selected_instr_info <- instr_info[, c('channel_id', 'fluorophores', 'laser', 'pmt')]
df <- df %>% append_dataframe(
    stranspose(selected_instr_info, colname='channel_id'),
    infront=TRUE, reset_index=FALSE
)


# ----------------------------------------------------------------------
# Postprocessing

# Rename columns from channel_id (1, 2, 3, etc.) to representative_fluorophore
fluorophores_for_channels = instr_info[['representative_fluorophore']]
names(fluorophores_for_channels) = instr_info[['channel_id']]
df <- rename_columns(df, fluorophores_for_channels[channel_cols])

# Rename last row
row.names(df)[nrow(df)] <- 'count'
df <- reset_index(df)  # export index information without frameshift

# save
if (!troubleshooting) {
    filepath = file.path(dirname(data_dir), 'antibodies_vs_channels.csv')
    write.table(df, file = filepath, row.names = FALSE, sep=',')
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
