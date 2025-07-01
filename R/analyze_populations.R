## Graphs percentage of live cells of each cell population

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('logr'))
library('optparse')
import::from(progress, 'progress_bar')
import::from(jsonlite, 'toJSON')
import::from(magrittr, '%>%')
import::from(stringr, 'str_detect')
import::from(tidyr, 'pivot_longer', 'pivot_wider')
import::from(ggplot2, 'ggsave')

import::from(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'preprocess_flowjo_export', 'sort_groups_by_metric', .character_only=TRUE)

import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'rename_columns', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'fill_missing_keys', 'items_in_a_not_b', 'multiple_replacement',
    .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'math.R'),
    'apply_unpaired_t_test', 'apply_multiple_comparisons', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'save_fig', 'plot_scatter', 'plot_violin',
    'plot_multiple_comparisons', .character_only=TRUE)
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

    make_option(c("-o", "--output-dir"), default="output",
                metavar="output", type="character",
                help="output directory for the data"),

    make_option(c("-m", "--metadata-dir"), default="data/flow-metadata",
                metavar="data/flow-metadata", type="character",
                help="metadata directory containing csvs mapping fcs to mouse ids"),

    make_option(c("-r", "--ref-dir"), default='data/mice',
                metavar='data/mice', type="character",
                help="directory of files containing all the mouse data"),

    make_option(c("-g", "--group-by"), default='treatment,genotype',
                metavar='treatment', type="character",
                help="enter a column or comma-separated list of columns, no spaces"),

    make_option(c("-s", "--stat"), default='fishers_lsd',
                metavar='fishers_lsd', type="character",
                help="Choose 'fishers_lsd', 't_test', 'tukey', or 'bonferroni'"),

    make_option(c("-p", "--plotly-overview"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="generate overview interactive violin plot"),

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
log <- log_open(paste0("analyze_populations-",
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
flow_metadata <- append_many_csv(
    file.path(wd, opt[['metadata-dir']]), recursive=TRUE, include_filepath=FALSE)
mouse_db <- append_many_csv(
    file.path(wd, opt[['ref-dir']]), recursive=TRUE, include_filepath=FALSE)


# ----------------------------------------------------------------------
# Preprocessing

log_print(paste(Sys.time(), 'Preprocessing...'))

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
for (cell_type in c(cell_type_ignore, 'mNeonGreen+')) {
    df <- df[!str_detect(df[['cell_type']], cell_type), ]
}


# inner join flow metadata
df <- merge(df, flow_metadata,
    by="fcs_name", all.x=FALSE, all.y=FALSE, suffixes=c('', '_'))
df <- df[(df[['is_unstained']]==FALSE), ]  # drop unstained cells
if (length(metadata_cols) > 1) {
    df[['group_name']] <- apply( df[ , metadata_cols ] , 1 , paste , collapse = ", " )
} else {
    df[['group_name']] <- df[[metadata_cols]]
}


# left join mouse data if available
if (!is.null(mouse_db)) {
    log_print(paste(Sys.time(), 'Merging mouse data...'))
    df <- merge(df, mouse_db[, items_in_a_not_b(colnames(mouse_db), mouse_db_ignore)],
        by="mouse_id", all.x=TRUE, all.y=FALSE, suffixes=c('', '_'))
    df[['weeks_old']] <- round(df[['age']]/7, 1)
}


# filter and sort
df <- df[(df[['num_cells']]>50), ]
df <- df[order(df[['cell_type']], df[['organ']], df[['group_name']]), ]  # sort rows
group_names <- sort(unique( df[['group_name']] ))


# ----------------------------------------------------------------------
# Compute statistics

log_print(paste(Sys.time(), paste0('Computing ', opt[['stat']], '...')))

if (opt[['stat']]=='t_test') {
    pval_tbl <- apply_unpaired_t_test(
        df,
        index_cols=c('organ', 'cell_type'),
        group_name='group_name',
        metric='pct_cells'
    )
} else {
    pval_tbl <- apply_multiple_comparisons(
        df,
        index_cols=c('organ', 'cell_type'),
        group_name='group_name',
        metric='pct_cells',
        correction=opt[['stat']]
    )
}

# save
if (!troubleshooting) {
    tmp <- pval_tbl
    for (col in group_names) {
        tmp[[col]] <- mapply(toJSON, tmp[[col]])
    }

    if (!dir.exists(file.path(wd, opt[['output-dir']], 'data'))) {
        dir.create(file.path(wd, opt[['output-dir']], 'data'), recursive=TRUE)
    }
    filepath = file.path(
        wd, opt[['output-dir']], 'data', paste0('pvals-', opt[['stat']], '.csv'))
    write.table(tmp, file = filepath, row.names = FALSE,  sep = ',' )
}


# ----------------------------------------------------------------------
# Data overview

organs <- sort(unique(df[['organ']]))
log_print(paste(Sys.time(), 'Groups found...', paste(organs, collapse = ', ')))

if (opt[['plotly-overview']]) {
    log_print(paste(Sys.time(), 'Exporting data with plots...'))
} else {
    log_print(paste(Sys.time(), 'Exporting data...'))
}

for (organ in sort(organs)) {

    log_print(paste0(Sys.time(), ' Processing ', organ, '...'))


    # ----------------------------------------------------------------------
    # Data

    tmp <- sort_groups_by_metric(
        df[(df[['organ']]==organ),
            unique(intersect(
                c('organ', 'genotype', 'treatment',
                   'group_name', metadata_cols, 'cell_type',
                   'Cells/Single Cells/Single Cells/Live Cells',
                   'num_cells', 'pct_cells', 
                   'mouse_id', 'sex', 'weeks_old', 'fcs_name'),
                colnames(df)
            ))],
        groups=c('group_name')
    )

    if (!troubleshooting) {
        dirpath <- file.path(wd, opt[['output-dir']], 'data', 'populations')
        if (!dir.exists(dirpath)) {
            dir.create(dirpath, recursive=TRUE)
        }
        filepath = file.path(dirpath,
            paste0(organ, '.csv')
        )
        write.table(tmp, file = filepath, row.names = FALSE, sep = ',')
    }


    # ----------------------------------------------------------------------
    # Plot

    if (opt[['plotly-overview']]) {

        fig <- plot_violin(
            df[(df[['organ']]==organ), ],
            x='cell_type', y='pct_cells', group_by='group_name',
            ylabel='Percent of Live Cells', title=organ,
            ymin=0, ymax=100,
            hover_data=unique(intersect(
                    c('mouse_id', 'group_name', metadata_cols,
                      'sex', 'treatment', 'weeks_old', 
                      'Cells', 'num_cells', 'fcs_name'),
                    colnames(df)
            ))
            # color_discrete_map=c(
            #     'heterozygous'='#2ca02c',  # green
            #     'wild type'='#62c1e5'  # blue
            # )
        )

        if (!troubleshooting) {
            save_fig(
                fig=fig,
                height=opt[['height']], width=opt[['width']],
                dirpath=file.path(wd, opt[['output-dir']], 'figures', 'overview'),
                filename=paste('violin-pct_cells', 
                    organ, gsub(',', '_', opt[['group-by']]), sep='-'),
                save_html=TRUE
            )
        }
    }
}


# ----------------------------------------------------------------------
# Violin with significance

log_print(paste(Sys.time(), paste(
    nrow(pval_tbl), 'populations found across',
    length(organs), 'organs') )
)

log_print(paste(Sys.time(), paste('Plotting mulitple comparisons using', opt[['stat']])))

pbar <- progress_bar$new(total = nrow(pval_tbl))
for (idx in 1:nrow(pval_tbl)) {

    # subset data
    organ <- pval_tbl[idx, 'organ'][[1]]
    cell_type <- pval_tbl[idx, 'cell_type'][[1]]
    df_subset <- df[
        (df[['organ']] == organ) & (df[['cell_type']] == cell_type),
        c('organ', 'cell_type', 'group_name', 'pct_cells')
    ]

    fig <- plot_multiple_comparisons(
        df_subset,
        x='group_name', y='pct_cells',
        ylabel='Percent of\nLive Cells',
        title=paste(toupper(organ), cell_type),
        test=opt[['stat']],
        show_numbers=opt[['show-numbers']]
    )

    # save
    if (!troubleshooting) {
        
        dirpath <- file.path(wd, opt[['output-dir']], 'figures',
            paste( opt[['stat']], gsub(',', '_', opt[['group-by']]), sep='-' ),
            organ
        )
        if (!dir.exists(dirpath)) {
            dir.create(dirpath, recursive=TRUE)
        }
        filepath = file.path(dirpath,
            paste0(opt[['stat']], '-', organ, '-', gsub(' ', '_', tolower(cell_type)), '.svg' )
        )

        withCallingHandlers({
            ggsave(
                filepath,  
                plot=fig,
                height=5000, width=5000, dpi=500, units='px', scaling=2 
            )
        }, warning = function(w) {
            if ( any(grepl("rows containing non-finite values", w),
                     grepl("fewer than two data points", w),
                     grepl("argument is not numeric or logical: returning NA", w)) ) {
                invokeRestart("muffleWarning")
            }
        })
    }

    pbar$tick()
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
