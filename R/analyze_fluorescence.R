## Compute cell proportions from each organ using exported flowjo data

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
suppressPackageStartupMessages(library('dplyr'))
library('optparse')
suppressPackageStartupMessages(library('logr'))
import::from(magrittr, '%>%')
import::from(stringr, 'str_detect')
import::from(tidyr, 'pivot_longer')

import::from(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'preprocess_flowjo_export', .character_only=TRUE)

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

    make_option(c("-g", "--group-by"), default='treatment',
                metavar='treatment', type="character",
                help="enter a column or comma-separated list of columns, no spaces"),

    make_option(c("-p", "--png-only"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="only save png and not HTML, useful for setting png height/width"),

    make_option(c("-l", "--height"), default=500,
                metavar="500", type="integer",
                help="height in px"),

    make_option(c("-w", "--width"), default=800,
                metavar="800", type="integer",
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
group_by <- unlist(strsplit(opt[['group-by']], ','))


# Start Log
start_time = Sys.time()
log <- log_open(paste0("analyze_fluorescence-",
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
    names_to = "gate", values_to = "num_cells",
    cols=items_in_a_not_b(colnames(counts), c(id_cols, initial_gates)),
    values_drop_na = TRUE
)
df[['cell_type']] <- unlist(lapply(strsplit(df[['gate']], '/'), function(x) x[length(x)]))
df[['cell_type']] <- multiple_replacement(df[['cell_type']], cell_type_spell_check)
df[['pct_cells']] <- round(df[['num_cells']] /
    df[['Cells/Single Cells/Single Cells/Live Cells']] * 100, 2)
# filter ignored cell types
for (cell_type in c(cell_type_ignore)) {
    df <- df[!str_detect(df[['cell_type']], cell_type), ]
}


# inner join flow metadata
df <- merge(df, flow_metadata,
    by="fcs_name", all.x=FALSE, all.y=FALSE, suffixes=c('', '_'))
df <- df[(df[['is_unstained']]==FALSE), ]  # drop unstained cells
if (length(group_by) > 1) {
    df[['group_by']] <- apply( df[ , group_by ] , 1 , paste , collapse = ", " )
} else {
    df[['group_by']] <- df[[group_by]]
}


# left join mouse data
df <- merge(df, mouse_db[, items_in_a_not_b(colnames(mouse_db), mouse_db_ignore)],
    by="mouse_id", all.x=TRUE, all.y=FALSE, suffixes=c('', '_'))
df[['weeks_old']] <- round(df[['age']]/7, 1)


# unpivot mNeonGreen+ num_cells into its own column
if ( any(str_detect(colnames(counts), 'mNeonGreen')) ) {

    # split out mNeonGreen+ rows
    fp_quant <- df[(df[['cell_type']]=='mNeonGreen+'), c(id_cols, 'gate', 'num_cells')]

    # second-to-last gate
    fp_quant[['gate']] <- unlist(lapply(
        strsplit(fp_quant[['gate']], '/'), function(x) paste0(x[1:length(x)-1], collapse='/')
    ))
    fp_quant[['cell_type']] <- unlist(lapply(
        strsplit(fp_quant[['gate']], '/'), function(x) x[length(x)]
    ))
    fp_quant <- rename_columns(fp_quant, c('num_cells'='num_mneongreen_pos'))
    
    df = merge(
        df[(df[['cell_type']]!='mNeonGreen+'), ],
        fp_quant[, c(id_cols, 'gate', 'num_mneongreen_pos')],
        by=c(id_cols, 'gate'), all.x=TRUE, all.y=FALSE
    )
    df[['pct_mneongreen_pos']] <- round(df[['num_mneongreen_pos']] / df[['num_cells']] * 100, 2)
}


# ----------------------------------------------------------------------
# Plot mNeonGreen positivity per cell type

log_print(paste(Sys.time(), 'Quantifying mNeonGreen...'))

organs <- sort(unique(df[['organ']]))
for (organ in sort(organs)) {
    if ( !all(is.na(df[(df[['organ']]==organ), 'pct_mneongreen_pos'])) ) {

        log_print(paste(Sys.time(), 'Processing...', organ))

        fig <- plot_violin(
            df[(df[['organ']]==organ) & (df[['num_cells']]>10), ],
            x='cell_type', y='pct_mneongreen_pos', group_by='zygosity',
            ylabel='Percent mNeonGreen+', title=organ,
            ymin=0, ymax=100,
            xaxis_angle=-90,
            hover_data=c('mouse_id', 'sex', 'zygosity', 'treatment', 'weeks_old',
                         'Cells', 'num_cells', 'num_mneongreen_pos', 'fcs_name'),
            color_discrete_map=c(
                'heterozygous'='#2ca02c',  # green
                'wild type'='#62c1e5'  # blue
            )
        )

        # export
        if (!troubleshooting) {
            save_fig(
                fig=fig,
                dirpath=file.path(wd, opt[['output-dir']], 'figures', 'mneongreen'),
                filename=paste0('violin-pct_mneongreen_pos-', organ),
                save_html=!opt[['png-only']]
            )
        }

        fig <- plot_scatter(
            df[((df[['organ']]==organ) & 
                (df[['num_cells']]>10) &
                (df[['zygosity']]=='heterozygous') &
                (df[['cell_type']] %in% c('Ly6C-hi Monocytes', 'Ly6C-int Monocytes', 'Neutrophils'))
                ), ],
            x='weeks_old', y='pct_mneongreen_pos', group_by='cell_type',
            xlabel='Age (Weeks)', ylabel='Percent mNeonGreen+', title=organ,
            ymin=0, ymax=100,
            hover_data=c('mouse_id', 'sex', 'zygosity', 'treatment', 'weeks_old',
                         'Cells', 'num_cells', 'num_mneongreen_pos', 'fcs_name'),
            color_discrete_map=c(
                'Ly6C-hi Monocytes'='#ff7f0e',  # orange
                'Ly6C-int Monocytes'='#2ca02c',  # green
                'Neutrophils'='#62c1e5'  # blue
            )
        )

        # export
        if (!troubleshooting) {
            save_fig(
                fig=fig,
                height=1000, width=1600, scale=3,
                dirpath=file.path(wd, opt[['output-dir']], 'figures', 'mneongreen'),
                filename=paste0('scatter-pct_mneongreen_pos-', organ),
                save_html=!opt[['png-only']]
            )
        }
    }
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
