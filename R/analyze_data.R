## Analyze data exported from flowjo

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
suppressPackageStartupMessages(library('dplyr'))
library('optparse')
library('logr')
# import::from(magrittr, '%>%')
import::from(stringr, 'str_detect')
import::from(tidyr, 'pivot_longer')
import::here(plyr, 'rbind.fill')
import::from(ggplot2, 'ggsave')

import::from(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'preprocess_flowjo_export', .character_only=TRUE)

import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'rename_columns', 'reset_index', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_dots', .character_only=TRUE)

import::from(file.path(wd, 'R', 'config', 'populations.R'),
    'populations_for_organ', .character_only=TRUE)

# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/flow-tables',
                metavar='data/flow-tables', type="character",
                help="set the input directory, all csv files will be read in"),

    make_option(c("-o", "--output-dir"), default="data/analysis",
                metavar="data/analysis", type="character",
                help="set the output directory for the data"),

    make_option(c("-f", "--figures-dir"), default="figures/analysis",
                metavar="figures/analysis", type="character",
                help="set the output directory for the figures"),

    make_option(c("-l", "--length"), default=800,
                metavar="800", type="integer",
                help="height in px"),

    make_option(c("-w", "--width"), default=1200,
                metavar="1200", type="integer",
                help="width in px"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]
output_dir <- file.path(wd, opt[['output-dir']])
troubleshooting_dir = file.path(output_dir, 'troubleshooting')


# # Start Log
# start_time = Sys.time()
# log <- log_open(paste0("analyze_data-",
#     strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
# log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Preprocessing

# Read counts data. Note that counts are in a wide format
counts_list <- append_many_csv(file.path(wd, opt[['input-dir']]), return_list=TRUE)
for (i in 1:length(counts_list)) {
    colnames(counts_list[[i]]) <- sub("^Macrophages \\|", "Cells |", colnames(counts_list[[i]]))
    colnames(counts_list[[i]]) <- sub("^Macrophages/", "Cells/", colnames(counts_list[[i]]))
}
counts <- do.call(rbind.fill, counts_list)
counts <- preprocess_flowjo_export(counts)


# ----------------------------------------------------------------------
# Plot percentages of each cell type

for (organ in c('bm', 'pb', 'pc', 'spleen')) {
    log_print(paste(Sys.time(), 'Processing...', organ))
    populations <- populations_for_organ[[organ]]
    cell_types <- unlist(lapply(strsplit(populations, '/'), function(x) x[length(x)]))

    # select data
    counts_subset <- counts[
        (counts['organ']==organ),
        c('mouse_id', 'organ', 'strain', 'treatment_group', populations)
    ]
    colnames(counts_subset) <- lapply(
        strsplit(colnames(counts_subset), '/'), function(x) x[length(x)]
    )  # get name of cell type from last gate
    counts_subset <- aggregate(counts_subset[, cell_types],
        by=list(
            'mouse_id'=counts_subset[['mouse_id']],
            'organ'=counts_subset[['organ']],
            'strain'=counts_subset[['strain']],
            'treatment_group'=counts_subset[['treatment_group']]
        ),
        FUN=function(x) sum(x, na.rm=TRUE)
    )

    # reshape
    df <- pivot_longer(counts_subset,
        names_to = "cell_type", values_to = "num_cells",
        cols=cell_types[2:length(cell_types)])
    df <- rename_columns(df, c('Live Cells'='live_cells'))
    df[['pct_cells']] <- df[['num_cells']] / df[['live_cells']] * 100


    # plot
    fig <- plot_dots(
        df,
        x='cell_type', y='pct_cells', color='treatment_group',
        title=organ
    )

    # export
    if (!troubleshooting) {
        if (!dir.exists(file.path(wd, opt[['figures-dir']]))) {
            dir.create(file.path(wd, opt[['figures-dir']]), recursive=TRUE)
        }
        ggsave(
            file.path(wd, opt[['figures-dir']],
                paste0('pct_cells-cell_type-treatment_group-', organ, '.png')),  # filename
            plot=fig,
            height=opt[['length']], width=opt[['width']], dpi=200,
            units="px", scaling=0.5
        )
    }
}


# end_time = Sys.time()
# log_print(paste('Script ended at:', Sys.time()))
# log_print(paste("Script completed in:", difftime(end_time, start_time)))
# log_close()


# ----------------------------------------------------------------------
# Troubleshooting

if (FALSE) {

    common_cols <- c(
        "filename",
        "fcs_name",
        "Cells",
        "Cells/Single Cells",
        "Cells/Single Cells/Single Cells",
        "Cells/Single Cells/Single Cells/Live Cells",

        # figure out how to combine this with Cells
        "Macrophages",
        "Macrophages/Single Cells",
        "Macrophages/Single Cells/Single Cells",
        "Macrophages/Single Cells/Single Cells/Live Cells",

        "metadata",
        "staining",
        "organ",
        "mouse_id",
        "treatment_group",
        "strain"
    )

    counts_long <- pivot_longer(counts,
        names_to = "gating", values_to = "num_cells",
        cols=items_in_a_not_b(colnames(counts), common_cols)
    )
    counts_long <- counts_long[
        ( rowSums(is.na(counts_long[, c('Cells', 'Macrophages')]))==1 ) & 
        ( !is.na(counts_long[, c('num_cells')]) ),
    ]
}
