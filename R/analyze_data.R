## Build the design matrix

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
suppressPackageStartupMessages(library('dplyr'))
library('optparse')
library('logr')
import::from(magrittr, '%>%')
import::from(stringr, 'str_detect')
import::from(tidyr, 'pivot_longer')

import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'rename_columns', 'reset_index', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_dots', .character_only=TRUE)

import::from(file.path(wd, 'R', 'config', 'populations.R'),
    'bm_populations', 'pb_populations', 'pc_populations', 'spleen_populations',
    .character_only=TRUE)

# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/flow-tables',
                metavar='data/flow-tables', type="character",
                help="set the input directory, all csv files will be read in"),

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
troubleshooting_dir = file.path(output_dir, 'troubleshooting')


# Start Log
start_time = Sys.time()
log <- log_open(paste0("analyze_data-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


populations_for_organ <- list2env(list(
    'bm'= bm_populations,
    'pb'= pb_populations,
    'pc'= pc_populations,
    'spleen'= spleen_populations
))


# ----------------------------------------------------------------------
# Preprocessing

# Read counts data
df <- append_many_csv(file.path(wd, opt[['input-dir']]))
df <- rename_columns(df, c('X1'='fcs_name'))
colnames(df) <- gsub('[[:space:]]\\|[[:space:]]Count', '', colnames(df))  # remove suffix
df <- df[!(df[['fcs_name']] %in% c('Mean', 'SD')), ]  # drop summary statistics
df <- df[!str_detect(df[['fcs_name']], 'unstained'), ]  # drop unstained cells

# extract metadata from the BD FACSDiva specimen name
df[['metadata']] <- as.character(
    lapply(strsplit(df[['fcs_name']], '_'),
    function(x) paste( x[4:length(x)-1], collapse='_' )) 
)
df[['organ']] <- as.character(lapply(strsplit(df[['metadata']], '-'), function(x) x[[2]]))
df[['mouse_id']] <- as.character(lapply(strsplit(df[['metadata']], '-'), function(x) x[[3]]))
df[['treatment_group']] <- as.character(lapply(strsplit(df[['metadata']], '-'), function(x) x[[4]]))
df[['strain']] <- as.character(strsplit(strsplit(df[['mouse_id']], '-')[[1]], '_')[[1]])
df <- reset_index(df, drop=TRUE)


# ----------------------------------------------------------------------
# Plot percentages of each cell type


for (organ in c('bm', 'pb', 'pc', 'spleen')) {
    populations <- populations_for_organ[[organ]]
    cell_types <- unlist(lapply(strsplit(populations, '/'), function(x) x[length(x)]))

    # preprocess data
    tmp <- df[
        (df['organ']==organ),
        c('mouse_id', 'organ', 'strain', 'treatment_group', populations)
    ]
    colnames(tmp) <- lapply(strsplit(colnames(tmp), '/'), function(x) x[length(x)])
    tmp <- rename_columns(tmp, c('Live Cells'='live_cells'))
    tmp <- pivot_longer(tmp,
        names_to = "cell_type", values_to = "num_cells",
        cols=cell_types[2:length(cell_types)])
    tmp[['pct_cells']] <- tmp[['num_cells']] / tmp[['live_cells']] * 100

    # plot
    fig <- plot_dots(
        tmp,
        x='cell_type', y='pct_cells', color='treatment_group',
        title=organ
    )
    print(fig)
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
