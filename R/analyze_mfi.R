## Graphs MFI for each cell population

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
suppressPackageStartupMessages(library('dplyr'))
library('optparse')
suppressPackageStartupMessages(library('logr'))
import::from(magrittr, '%>%')
import::from(stringr, 'str_detect', 'regex')
import::from(tidyr, 'pivot_longer')

import::from(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'preprocess_flowjo_export', 'sort_for_graphpad', .character_only=TRUE)

import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'rename_columns', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', 'multiple_replacement',
    .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'save_fig', 'plot_scatter', 'plot_violin', .character_only=TRUE)
import::from(file.path(wd, 'R', 'config', 'flow.R'),
    'id_cols', 'initial_gates', 'cell_type_spell_check', 'cell_type_ignore',
    'mouse_db_ignore', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/flow-data',
                metavar='data/flow-data', type="character",
                help="input directory, all csv files will be read in"),

    make_option(c("-m", "--metadata-dir"), default="data/flow-metadata",
                metavar="data/flow-metadata", type="character",
                help="metadata directory containing csvs mapping fcs to mouse ids"),

    make_option(c("-o", "--output-dir"), default="output",
                metavar="output", type="character",
                help="output directory for the data"),

    make_option(c("-r", "--ref-dir"), default='data/mice',
                metavar='data/mice', type="character",
                help="directory of files containing all the mouse data"),

    make_option(c("-g", "--group-by"), default='sex,genotype,treatment',
                metavar='treatment', type="character",
                help="enter a column or comma-separated list of columns, no spaces"),

    make_option(c("-p", "--png-only"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="only save png and not HTML, useful for setting png height/width"),

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

# args
groupby <- unlist(strsplit(opt[['group-by']], ','))


# Start Log
start_time = Sys.time()
log <- log_open(paste0("quantify_mfi-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Raw Data

log_print(paste(Sys.time(), 'Reading data...'))

# Read counts data exported from flowjo
counts <- append_many_csv(file.path(wd, opt[['input-dir']]), recursive=TRUE)
if (is.null(counts)) {
    msg <- paste("No counts data found. Please check", opt[['input-dir']], '...')
    stop(msg)
}
counts <- preprocess_flowjo_export(counts)
counts <- counts[!str_detect(counts[['fcs_name']], '(U|u)nstained'), ]  # drop unstained cells

# reference files
mouse_db <- append_many_csv(file.path(wd, opt[['ref-dir']]), recursive=TRUE, include_filepath=FALSE)
flow_metadata <- append_many_csv(file.path(wd, opt[['metadata-dir']]), recursive=TRUE, include_filepath=FALSE)


# ----------------------------------------------------------------------
# Preprocessing

# Reshape so each row is a gate in each sample
df <- pivot_longer(counts,
    names_to = "gate", values_to = "mfi",
    cols=items_in_a_not_b(colnames(counts), c(id_cols, initial_gates)),
    values_drop_na = TRUE
)
df[['gate']] <- unlist(lapply(strsplit(df[['gate']], ' \\| '), function(x) x[1]))
df[(df[['gate']]=='Geometric Mean (Comp-Alexa Fluor 488-A)'), 'gate'] <- 'Ungated'
df[['cell_type']] <- unlist(lapply(strsplit(df[['gate']], '/'), function(x) x[length(x)]))
df[['cell_type']] <- multiple_replacement(df[['cell_type']], cell_type_spell_check)

# filter ignored cell types
for (cell_type in c(cell_type_ignore,
    'mNeonGreen', 'mNeonGreen', 'Ungated', regex("^Cells$"), regex("^Single Cells$"))
    ) {
        df <- df[!str_detect(df[['cell_type']], cell_type), ]
        
    }

# inner join flow metadata
df <- merge(df, flow_metadata,
    by="fcs_name", all.x=FALSE, all.y=FALSE, suffixes=c('', '_'))
df <- df[(df[['is_unstained']]==FALSE), ]  # drop unstained cells
if (length(groupby) > 1) {
    df[['groupby']] <- apply( df[ , groupby ] , 1 , paste , collapse = ", " )
} else {
    df[['groupby']] <- df[[groupby]]
}


# left join mouse data
df <- merge(df, mouse_db[, items_in_a_not_b(colnames(mouse_db), mouse_db_ignore)],
    by="mouse_id", all.x=TRUE, all.y=FALSE, suffixes=c('', '_'))
df[['weeks_old']] <- round(df[['age']]/7, 1)


# # save
# if (!troubleshooting) {
#     if (!dir.exists(file.path(wd, opt[['output-dir']], 'data'))) {
#         dir.create(file.path(wd, opt[['output-dir']], 'data'), recursive=TRUE)
#     }
#     filepath = file.path(wd, opt[['output-dir']], 'data', 'populations.csv')
#     write.table(df, file = filepath, row.names = FALSE, sep = ',' )
# }


# ----------------------------------------------------------------------
# Plot cell type frequency per organ

log_print(paste(Sys.time(), 'Quantifying cell type frequencies...'))

organs <- sort(unique(df[['organ']]))
log_print(paste(Sys.time(), 'Groups found...', paste(organs, collapse = ', ')))
for (organ in sort(organs)) {

    log_print(paste(Sys.time(), 'Processing...', organ))


    # ----------------------------------------------------------------------
    # Plot

    legend_order <- c(
        "Female, WT, DMSO",
        "Female, mNG+, DMSO",
        "Female, WT, R848",
        "Female, mNG+, R848",
        "Male, WT, DMSO",
        "Male, mNG+, DMSO",
        "Male, WT, R848",
        "Male, mNG+, R848"
    )

    fig <- plot_violin(
        df[(df[['organ']]==organ), ],
        x='cell_type', y='mfi', group_by='groupby',
        legend_order=legend_order,
        ylabel='MFI', title=organ,
        ymin=0,
        hover_data=unique(c(
            'mouse_id', 'groupby', groupby, 'sex', 'treatment', 'weeks_old', 
            'Cells', 'num_cells', 'fcs_name')),
        sort=TRUE
    )

    if (!troubleshooting) {
        save_fig(
            fig=fig,
            height=opt[['height']], width=opt[['width']],
            dirpath=file.path(wd, opt[['output-dir']], 'figures', 'mfi'),
            filename=paste('violin-pct_cells', organ, gsub(',', '_', opt[['group-by']]), sep='-'),
            save_html=!opt[['png-only']]
        )
    }


    # ----------------------------------------------------------------------
    # Data

    tmp <- sort_for_graphpad(
        df[(df[['organ']]==organ),
        c('cell_type', 'groupby', groupby, 'mouse_id',
          'organ', 'sex', 'treatment', 'weeks_old', 'fcs_name', 'mfi')],
        y='mfi',
        groups=c('groupby')
    )

    if (!troubleshooting) {
        dirpath <- file.path(wd, opt[['output-dir']], 'data', 'mfi')
        if (!dir.exists(dirpath)) {
            dir.create(dirpath, recursive=TRUE)
        }
        filepath = file.path(dirpath,
            paste0(organ, '-', gsub(',', '_', opt[['group-by']]), '.csv')
        )
        write.table(tmp, file = filepath, row.names = FALSE, sep = ',')
    }

}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
