wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('logr'))
library('optparse')
import::from(progress, 'progress_bar')
import::from(jsonlite, 'toJSON')
import::from(stringr, 'str_detect')
import::from(tidyr, 'pivot_longer')
import::from(ggplot2, 'ggsave')

import::from(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'find_initial_gates', 'preprocess_flowjo_export', 'sort_groups_by_metric',
    .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'file_io.R'),
    'import_flowjo_export', 'import_flow_metadata', .character_only=TRUE)

import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', 'multiple_replacement', 'move_list_items_to_front',
    .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'math.R'),
    'apply_unpaired_t_test', 'apply_multiple_comparisons', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'save_fig', 'plot_violin', 'plot_multiple_comparisons', .character_only=TRUE)
import::from(file.path(wd, 'R', 'config', 'flow.R'),
    'mouse_db_ignore', .character_only=TRUE)
import::from(file.path(wd, 'R', 'config', 'user_input.R'),
    'custom_group_order', .character_only=TRUE)

reticulate::py_config()  # required on linux to access reticulate


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/weights',
                metavar='data/weights', type="character",
                help="input directory, all csv files will be read in"),

    make_option(c("-o", "--output-dir"), default="hfc-output",
                metavar="hfc-output", type="character",
                help="output directory for the data"),

    make_option(c("-m", "--metadata-dir"), default="data/flow-metadata",
                metavar="data/flow-metadata", type="character",
                help="metadata directory containing csvs mapping fcs to mouse ids"),

    make_option(c("-r", "--ref-dir"), default='data/mice',
                metavar='data/mice', type="character",
                help="directory of files containing all the mouse data"),

    make_option(c("-g", "--group-by"), default='genotype,treatment',
                metavar='genotype,treatment', type="character",
                help="enter a column or comma-separated list of columns, no spaces"),

    make_option(c("-s", "--stat"), default='fishers_lsd',
                metavar='fishers_lsd', type="character",
                help="Choose 'fishers_lsd', 't_test', 'tukey', or 'bonferroni'"),

    make_option(c("-n", "--show-numbers"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="in violin with significance, show all numbers instead of stars"),

    make_option(c("-l", "--height"), default=500,
                metavar="500", type="integer",
                help="height in px"),

    make_option(c("-w", "--width"), default=2000,
                metavar="2000", type="integer",
                help="width in px, max width is 200000"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]
output_dir <- file.path(wd, opt[['output-dir']])
troubleshooting_dir = file.path(output_dir, 'troubleshooting')

if (!(opt[['stat']] %in% c('fishers_lsd', 't_test', 'tukey', 'bonferroni'))) {
    stop("Choose from stat= 'fishers_lsd', 't_test', 'tukey', or 'bonferroni'")
}

# args
metadata_cols <- unlist(strsplit(opt[['group-by']], ','))

# Start Log
start_time = Sys.time()
log <- log_open(paste0("analyze_counts-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Raw Data

log_print(paste(Sys.time(), 'Reading data...'))

# Read counts data exported from flowjo
df <- append_many_csv(
    file.path(wd, opt[['input-dir']]), recursive=TRUE, include_filepath=FALSE
)

if (length(metadata_cols) > 1) {
    df[['group_name']] <- apply( df[ , metadata_cols ] , 1 , paste , collapse = ", " )
} else {
    df[['group_name']] <- df[[metadata_cols]]
}


# ----------------------------------------------------------------------
# Compute statistics

log_print(paste(Sys.time(), paste0('Computing ', opt[['stat']], '...')))

# if (opt[['stat']]=='t_test') {
#     pval_tbl <- apply_unpaired_t_test(
#         df,
#         index_cols=c('organ'),
#         group_name='group_name',
#         metric='weight'
#     )
# } else {
#     pval_tbl <- apply_multiple_comparisons(
#         df,
#         index_cols=c('organ'),
#         group_name='group_name',
#         metric='weight',
#         correction=opt[['stat']]
#     )
# }


withCallingHandlers({
    fig <- plot_multiple_comparisons(
        df,
        x='group_name', y='weight',
        ymin=0,
        ylabel="Weight (g)",
        title="Spleen Weight",
        test=opt[['stat']],
        show_numbers=opt[['show-numbers']],
        custom_group_order=c('WT, DMSO', 'KO, DMSO', 'WT, R848', 'KO, R848')
    )
}, warning = function(w) {
    if ( any(grepl("containing non-finite values", w),
             grepl("outside the scale range", w)) ) {
        invokeRestart("muffleWarning")
    }
})

n_groups <- length(unique(df[['group_name']]))

# save
if (!troubleshooting) {
    
    dirpath <- file.path(wd, opt[['output-dir']], 'figures')
    if (!dir.exists(dirpath)) { dir.create(dirpath, recursive=TRUE) }
    filepath = file.path(dirpath, 'spleen_weight.svg')
    withCallingHandlers({
        ggsave(
            filepath,  
            plot=fig,
            height=5000, width=5000, dpi=500, units='px',
            scaling=(if (n_groups==2) {1} else 2)
        )
    }, warning = function(w) {
        if ( any(grepl("rows containing non-finite values", w),
                 grepl("fewer than two data points", w),
                 grepl("argument is not numeric or logical: returning NA", w)) ) {
            invokeRestart("muffleWarning")
        }
    })
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
