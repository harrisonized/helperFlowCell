## Compute cell proportions of different organs from flowjo

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
suppressPackageStartupMessages(library('dplyr'))
library('optparse')
suppressPackageStartupMessages(library('logr'))
import::from(magrittr, '%>%')
import::from(stringr, 'str_detect', 'str_extract')
import::from(tidyr, 'pivot_longer')
import::from(plyr, 'rbind.fill')
import::from(dplyr, 'group_by', 'summarize')
import::from(ggplot2, 'ggsave')
import::from(plotly, 'save_image')
import::from(htmlwidgets, 'saveWidget')  # brew install pandoc

import::from(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'preprocess_flowjo_export', .character_only=TRUE)

import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'rename_columns', 'reset_index', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'dict_zip', 'interleave', 'items_in_a_not_b', 'multiple_replacement',
    .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_dots', 'plot_violin', .character_only=TRUE)
import::from(file.path(wd, 'R', 'config', 'flow.R'),
    'id_cols', 'numerical_cols', 'ignored_cell_types', 'cell_type_replacements',
    .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/flow-data',
                metavar='data/flow-data', type="character",
                help="set the input directory, all csv files will be read in"),

    make_option(c("-o", "--output-dir"), default="data/analysis",
                metavar="data/analysis", type="character",
                help="set the output directory for the data"),

    make_option(c("-f", "--figures-dir"), default="figures/analysis",
                metavar="figures/analysis", type="character",
                help="set the output directory for the figures"),

    make_option(c("-s", "--save-html"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="save html files in addition to PNG"),

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


# ----------------------------------------------------------------------
# Preprocessing

# Read counts data exported from flowjo
log_print(paste(Sys.time(), 'Reading data...'))
counts <- append_many_csv(file.path(wd, opt[['input-dir']]), recursive=TRUE)
counts <- preprocess_flowjo_export(counts)
counts <- counts[, colnames(counts)[!str_detect(colnames(counts), 'mNeonGreen')]]  # drop mNeonGreen columns
counts <- counts[!str_detect(counts[['fcs_name']], '(U|u)nstained'), ]  # drop unstained cells

# Reshape
df <- pivot_longer(counts,
    names_to = "gate", values_to = "num_cells",
    cols=items_in_a_not_b(colnames(counts), c(id_cols, numerical_cols)),
    values_drop_na = TRUE
)
df[['organ']] <- unlist(lapply(strsplit(df[['filename']], '-'), function(x) x[1]))
df[['cell_type']] <- unlist(lapply(strsplit(df[['gate']], '/'), function(x) x[length(x)]))
df[['cell_type']] <- multiple_replacement(df[['cell_type']], cell_type_replacements)
for (ignored_cell_type in ignored_cell_types) {
    df <- df[!str_detect(df[['cell_type']], ignored_cell_type), ]
}
df[['pct_cells']] <- (df[['num_cells']] /
    df[['Cells/Single Cells/Single Cells/Live Cells']] * 100
)

# reorder columns
df <- df[, c(
    id_cols, 'organ', 'gate', 'cell_type',
    numerical_cols, 'num_cells', 'pct_cells')
]


# ----------------------------------------------------------------------
# Plot

organs <- sort(unique(df[['organ']]))
log_print(paste(Sys.time(), 'Groups found...', paste(organs, collapse = ', ')))
for (organ in sort(organs)) {
    log_print(paste(Sys.time(), 'Processing...', organ))

    fig <- plot_violin(
        df[(df[['organ']]==organ), ],
        x='cell_type', y='pct_cells', group_by=NULL,
        ylabel='Percent of Live Cells',
        ymin=0, ymax=100,
        title=organ
    )

    # export
    if (!troubleshooting) {

        if (!dir.exists( file.path(wd, opt[['figures-dir']], 'cell_proportions') )) {
            dir.create( file.path(wd, opt[['figures-dir']], 'cell_proportions'), recursive=TRUE)
        }

        # save PNG
        suppressWarnings(save_image(fig,
            file=file.path(wd, opt[['figures-dir']], 'cell_proportions',
                paste0('violin-pct_cells-', organ, '.png')),  # filename
            height=500, width=800, scale=3
        ))

        # save HTML
        if (opt[['save-html']]) {
            if (!dir.exists( file.path(wd, opt[['figures-dir']], 'cell_proportions', 'html') )) {
                dir.create( file.path(wd, opt[['figures-dir']], 'cell_proportions', 'html'), recursive=TRUE)
            }
            suppressWarnings(saveWidget(
                widget = fig,
                file=file.path(wd, opt[['figures-dir']], 'cell_proportions', 'html',
                    paste0('violin-pct_cells-', organ, '.html')),  # filename
                selfcontained = TRUE
            ))
            unlink(file.path(
                wd, opt[['figures-dir']], 'cell_proportions', 'html',
                paste0('violin-pct_cells-', organ, '_files')
            ), recursive=TRUE)            
        }
    }

}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
