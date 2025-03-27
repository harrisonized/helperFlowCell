## Compute cell proportions from each organ using exported flowjo data

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
suppressPackageStartupMessages(library('dplyr'))
library('optparse')
suppressPackageStartupMessages(library('logr'))
import::from(magrittr, '%>%')
import::from(stringr, 'str_detect')
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
    'save_fig', 'plot_violin', .character_only=TRUE)
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

    make_option(c("-g", "--group-by"), default='zygosity',
                metavar='zygosity', type="character",
                help="Ex: zygosity, genotype, ..."),

    make_option(c("-f", "--fluorescence"), default='mNeonGreen',
                metavar='mNeonGreen', type="character",
                help="Ex: mNeonGreen, GFP, YFP, ..."),

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
fp_pos <- paste0(opt[['fluorescence']],'+')
num_fp_pos <- paste0('num_', tolower(opt[['fluorescence']]),'_pos')
pct_fp_pos <- paste0('pct_', tolower(opt[['fluorescence']]),'_pos')


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
    df[['Cells/Single Cells/Single Cells/Live Cells']] * 100, 4)
# filter ignored cell types
for (cell_type in c(cell_type_ignore)) {
    df <- df[!str_detect(df[['cell_type']], cell_type), ]
}


# inner join flow metadata
df <- merge(df, flow_metadata,
    by="fcs_name", all.x=FALSE, all.y=FALSE, suffixes=c('', '_'))
df <- df[(df[['is_unstained']]==FALSE), ]  # drop unstained cells


# left join mouse data
df <- merge(df, mouse_db[, items_in_a_not_b(colnames(mouse_db), mouse_db_ignore)],
    by="mouse_id", all.x=TRUE, all.y=FALSE, suffixes=c('', '_'))
df[['weeks_old']] <- round(df[['age']]/7, 1)


# unpivot mNeonGreen+ num_cells into its own column
if ( any(str_detect(colnames(counts), opt[['fluorescence']])) ) {

    # split out mNeonGreen+ rows
    fp_quant <- df[
        (df[['cell_type']]==fp_pos),
        c(id_cols, 'gate', 'num_cells')
    ]

    # second-to-last gate
    fp_quant[['gate']] <- unlist(lapply(
        strsplit(fp_quant[['gate']], '/'), function(x) paste0(x[1:length(x)-1], collapse='/')
    ))
    fp_quant[['cell_type']] <- unlist(lapply(
        strsplit(fp_quant[['gate']], '/'), function(x) x[length(x)]
    ))
    fp_quant <- rename_columns(fp_quant, c('num_cells'=num_fp_pos))
    
    df = merge(
        df[(df[['cell_type']]!=fp_pos), ],
        fp_quant[, c(id_cols, 'gate', num_fp_pos)],
        by=c(id_cols, 'gate'), all.x=TRUE, all.y=FALSE
    )
    df[[pct_fp_pos]] <- round(df[[num_fp_pos]] / df[['num_cells']] * 100, 4)
}


# ----------------------------------------------------------------------
# Plot mNeonGreen positivity per cell type

log_print(paste(Sys.time(), 'Quantifying mNeonGreen...'))

organs <- sort(unique(df[['organ']]))
for (organ in sort(organs)) {
    if ( !all(is.na(df[(df[['organ']]==organ), num_fp_pos])) ) {

        log_print(paste(Sys.time(), 'Processing...', organ))


        # ----------------------------------------------------------------------
        # Plot

        fig <- plot_violin(
            df[(df[['organ']]==organ), ],
            x='cell_type', y=pct_fp_pos, group_by=opt[['group-by']],
            ylabel=paste('Percent', fp_pos), title=organ,
            ymin=0, ymax=100,
            xaxis_angle=-90,
            hover_data=unique(
                c('mouse_id', 'sex', opt[['group-by']], 'treatment', 'weeks_old',
                  'Cells', 'num_cells', num_fp_pos, 'fcs_name')
            )
            # color_discrete_map=c(
            #     'heterozygous'='#2ca02c',  # green
            #     'wild type'='#62c1e5'  # blue
            # )
        )

        # export
        if (!troubleshooting) {
            save_fig(
                fig=fig,
                dirpath=file.path(wd, opt[['output-dir']], 'figures', tolower(opt[['fluorescence']])),
                filename=paste0('violin-pct_', tolower(opt[['fluorescence']]), '_pos-', organ),
                save_html=!opt[['png-only']]
            )
        }


        # ----------------------------------------------------------------------
        # Data

        tmp <- sort_for_graphpad(
            df[(df[['organ']]==organ) & (!is.na(df[[pct_fp_pos]])),
                unique(
                    c('cell_type', opt[['group-by']], 'zygosity', 'sex', 'treatment',
                      'mouse_id', pct_fp_pos, num_fp_pos, 'num_cells', 'pct_cells',
                      'organ', 'sex', 'treatment', 'weeks_old', 'fcs_name')
                )
            ],
            groups=c('zygosity', 'sex', 'treatment')
        )

        if (!troubleshooting) {
            dirpath <- file.path(wd, opt[['output-dir']], 'data', tolower(opt[['fluorescence']]))
            if (!dir.exists(dirpath)) {
                dir.create(dirpath, recursive=TRUE)
            }
            filepath = file.path(dirpath,
                paste0(organ, '-', gsub(',', '_', opt[['group-by']]), '.csv')
            )
            write.table(tmp, file = filepath, row.names = FALSE, sep = ',')
        }

    }
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
