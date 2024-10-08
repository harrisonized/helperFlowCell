## Build the design matrix

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
suppressPackageStartupMessages(library('dplyr'))
library("openxlsx")
library('optparse')
suppressPackageStartupMessages(library('logr'))
import::from(magrittr, '%>%')
import::from(tidyr, 'separate_rows', 'pivot_wider')
import::from(reshape2, 'melt')
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'read_excel_or_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'fillna', 'reset_index',
    'append_dataframe', 'dataframe_row_from_named_list', 'stranspose',
    .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'preprocessing.R'),
   'preprocess_instrument_config', 'preprocess_antibody_inventory',
    .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-file"), default='data/design/panel.csv',
                metavar='data/design/panel.csv', type="character",
                help="specify the antibodies to use in your flow panel"),

    make_option(c("-a", "--antibody-inventory"), default='ref/antibody_inventory.xlsx',
                metavar='ref/antibody_inventory.xlsx', type="character",
                help="antibody inventory table"),

    make_option(c("-c", "--instrument-config"), default='ref/instrument_config.xlsx',
                metavar='ref/instrument_config.xlsx', type="character",
                help="instrument configuration file"),

    make_option(c("-o", "--output-dir"), default="data/design",
                metavar="data/design", type="character",
                help="set the output directory for the data"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]
output_dir <- file.path(wd, opt[['output-dir']])
troubleshooting_dir = file.path(output_dir, 'troubleshooting')


# Start Log
start_time = Sys.time()
log <- log_open(paste0("build_design_matrix-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Configuration files

# Instrument config
instr_cfg <- read_excel_or_csv(file.path(wd, opt[['instrument-config']]))
instr_cfg <- preprocess_instrument_config(instr_cfg)
instr_cfg_long <- separate_rows(instr_cfg, 'fluorophore', sep=', ')

# Antibody inventory
ab_inv <- read_excel_or_csv(file.path(wd, opt[['antibody-inventory']]))
ab_inv <- preprocess_antibody_inventory(ab_inv)


# ----------------------------------------------------------------------
# Build design_matrix

# Read data
panel <- read_excel_or_csv(file.path(wd, opt[['input-file']]))
ab_shortlist <- merge(panel['antibody'], ab_inv, by='antibody',
    all.x=FALSE, all.y=FALSE, na_matches = 'never', sort=FALSE)
ab_shortlist <- ab_shortlist[, colnames(ab_inv)]

antibodies_not_found = sort(items_in_a_not_b(
    unique( panel[['antibody']] ), unique( ab_inv[['antibody']] )
))

# save
if (!troubleshooting) {
    if (!dir.exists(file.path(output_dir))) {
        dir.create(file.path(output_dir), recursive=TRUE)
    }
    filepath = file.path(output_dir, 
        paste0(tools::file_path_sans_ext(basename(opt[['antibody-inventory']])),
               '-shortlist', '.csv')
    )
    write.table(ab_shortlist, file = filepath, na = "", row.names = FALSE, sep=',')

    if (length(antibodies_not_found) > 0) {
        if (!dir.exists(file.path(troubleshooting_dir))) {
            dir.create(file.path(troubleshooting_dir), recursive=TRUE)
        }
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
        fluorophores = toString(unique(fluorophore)),  # fluorophores not sorted
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

# order by excitation, then emission
channel_order <- instr_cfg[
    order(instr_cfg[['excitation']], instr_cfg[['bandpass_filter']]),
    'channel_id'
][[1]]
design_matrix <- design_matrix[ ,
    c(intersect(channel_order, colnames(design_matrix)),
      items_in_a_not_b(colnames(design_matrix), channel_order)
    )
]

# Move channel_cols to the middle
id_cols <- c('antibody', 'num_channels_per_ab')
channel_cols <- items_in_a_not_b(colnames(design_matrix), c(id_cols, 'other'))
design_matrix <- design_matrix[ , c(id_cols, channel_cols, 'other')]


# ----------------------------------------------------------------------
# Postprocessing

# Add comma-separated fluorophore list to each row
row_metadata <- ab_counts %>%
    group_by(antibody) %>%
    summarise(fluorophores = toString(unique(most_common_fluorophore)))
design_matrix[['antibody_']] <- design_matrix[['antibody']]
design_matrix = merge(design_matrix, row_metadata,
    by='antibody', suffixes=c('', '_'),
    all.x=TRUE, all.y=FALSE, na_matches = 'never', sort=FALSE
)

# Add available fluorophores and instrument info to each column
col_metadata <- ab_counts %>%
    group_by(channel_id) %>%
    summarise(all_fluorophores = toString(unique(most_common_fluorophore))) %>%
    merge(instr_cfg[, c('channel_id', 'bandpass_filter', 'excitation', 'laser')],
          by='channel_id', suffixes=c('', '_'),
          all.x=TRUE, all.y=FALSE, na_matches = 'never', sort=FALSE)
design_matrix <- design_matrix %>% append_dataframe(
    stranspose(col_metadata, colname='channel_id'),
    infront=TRUE, reset_index=FALSE
)

# Rename last row
row.names(design_matrix)[nrow(design_matrix)] <- 'count'
design_matrix <- reset_index(design_matrix)  # export index information without frameshift
# design_matrix[design_matrix==0] <- NA


# ----------------------------------------------------------------------
# Build Excel file


# table properties
row_start <- 5  # 4 column headers
row_end <- nrow(design_matrix)-1  # last row is count
col_start <- 4  # 3 row headers
col_end <- ncol(design_matrix)-2  # last columns are 'antibody_' and 'fluorophore'
nrows <- row_end-row_start+1
ncols <- col_end-col_start+1


wb <- createWorkbook()
addWorksheet(wb, "Sheet 1", gridLines = TRUE)


# add colors
submatrix_long <- melt(
    reset_index(design_matrix[row_start:row_end, col_start:col_end], 'rowname'),
    id.vars=c('rowname'), variable.name='colname'
)
submatrix_long[['row']] <- rep(1:nrows, times=ncols, each=1)
submatrix_long[['col']] <- rep(1:ncols, times=1, each=nrows)
pos_cells <- submatrix_long[(submatrix_long['value'] > 0), c('row', 'col')]


# design_matrix[(design_matrix ==0)] <- NA
styling_object <- createStyle(fgFill = "#C3C3C3")  # light gray
addStyle(
    wb,
    sheet = 1,
    style = styling_object,
    rows = pos_cells[['row']]+row_start,
    cols = pos_cells[['col']]+col_start-1
)


# set cell dimensions
setColWidths(wb, sheet = 1,
    cols = seq(col_start, col_end, 1),
    widths = 5.83,  # 40 pixels
    ignoreMergedCells = TRUE
)
setRowHeights(
    wb, sheet = 1,
    rows = seq(row_start+1, row_end+2, 1),
    heights = 20.00  # 20 pixels
)

writeData(wb, sheet = 1, design_matrix, rowNames = FALSE)


# save
if (!troubleshooting) {
    filepath = file.path(output_dir, 'design_matrix.xlsx')
    saveWorkbook(wb, filepath, overwrite  = TRUE)
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()


# ----------------------------------------------------------------------
# Troubleshooting only

if (FALSE) {

    if (!dir.exists(file.path(troubleshooting_dir))) {
        dir.create(file.path(troubleshooting_dir), recursive=TRUE)
    }

    # troubleshoot fluorophore naming
    filepath = file.path(troubleshooting_dir, 
        paste0('_', tools::file_path_sans_ext(basename(opt[['instrument-config']])), '.csv')
    )
    write.table(instr_cfg, file = filepath, na = "", row.names = FALSE, sep=',')

    # all fluorophores
    all_fluorophores = sort(unique( ab_inv[['fluorophore']] ))
    filepath = file.path(troubleshooting_dir, 'all_fluorophores.txt')
    write.table(all_fluorophores, filepath,
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    # troubleshoot fluorophore naming
    unavailable_fluorophores = sort(items_in_a_not_b(
        unique( ab_inv[['fluorophore']] ), unique( instr_cfg_long[['fluorophore']] )
    ))
    filepath = file.path(troubleshooting_dir, 'unavailable_fluorophores.txt')
    write.table(unavailable_fluorophores, filepath,
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    # all antibodies
    all_antibodies = sort(unique( ab_inv[['antibody']] ))
    filepath = file.path(troubleshooting_dir, 'all_antibodies.txt')
    write.table(all_antibodies, filepath, na = "",
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    # troubleshoot antibody naming
    filepath = file.path(troubleshooting_dir, 'antibodies.txt')
    write.table(sort(unique(ab_inv[['antibody']])), filepath,
                row.names = FALSE, col.names = FALSE, quote = FALSE)

}
