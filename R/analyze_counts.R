## Graphs percentage of live cells for each cell population
## Expects the first four gates to be Cells/Single Cells/Single Cells/Live Cells
## If abs_count is included as a column in your metadata file, this script
## will also calculate absolute counts

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


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/flow-counts',
                metavar='data/flow-counts', type="character",
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

    make_option(c("-g", "--group-by"), default='sex,treatment,zygosity',
                metavar='treatment', type="character",
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

# get last_initial_gate
last_initial_gate <- 'Live Cells'  # try CD45+

# Read counts data exported from flowjo
df <- import_flowjo_export(
    file.path(wd, opt[['input-dir']]),
    metric_name='num_cells',
    last_initial_gate=last_initial_gate  # change to CD45+
)
if (length(df)==0) {
    msg <- paste("No data found. Please check", file.path(wd, opt[['input-dir']]), '...')
    stop(msg)
}

divisor <- tail(find_initial_gates(colnames(df), last_initial_gate), 1)[[1]]
df[['pct_cells']] <- round(df[['num_cells']] /
    df[[divisor]] * 100, 4)

# reference files
flow_metadata <- import_flow_metadata(file.path(wd, opt[['metadata-dir']]))
mouse_db <- append_many_csv(
    file.path(wd, opt[['ref-dir']]),recursive=TRUE, include_filepath=FALSE)


# ----------------------------------------------------------------------
# Preprocessing

log_print(paste(Sys.time(), 'Preprocessing...'))

# inner join flow metadata
df <- merge(df, flow_metadata,
    by="fcs_name", all.x=FALSE, all.y=FALSE, suffixes=c('', '_'))
df <- df[(df[['is_unstained']]==FALSE), ]  # drop unstained cells
if (length(metadata_cols) > 1) {
    df[['group_name']] <- apply( df[ , metadata_cols ] , 1 , paste , collapse = ", " )
} else {
    df[['group_name']] <- df[[metadata_cols]]
}

# calculate
if ('total_viable_cells' %in% colnames(df)) {
    df[['abs_count']] <- round(df[['total_viable_cells']] * df[['pct_cells']]/100, 0)
}

# left join mouse data if available
if (!is.null(mouse_db)) {
    log_print(paste(Sys.time(), 'Merging mouse data...'))
    df <- merge(df, mouse_db[, items_in_a_not_b(colnames(mouse_db), mouse_db_ignore)],
        by="mouse_id", all.x=TRUE, all.y=FALSE, suffixes=c('', '_'))
    df[['weeks_old']] <- round(df[['age']]/7, 1)
}


# filter and sort
df <- df[order(df[['cell_type']], df[['organ']], df[['group_name']]), ]  # sort rows
group_names <- sort(unique( df[['group_name']] ))
df <- df[(df[['num_cells']]>10), ]

# exclude negative gates
cell_types <- unique(df[['cell_type']])
negative_gates <- cell_types[grepl("-$", cell_types)]
df <- df[((df[['cell_type']] %in% negative_gates)==FALSE), ]


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

    dirpath <- file.path(wd, opt[['output-dir']], 'data',
        gsub(',', '_', opt[['group-by']]), 'pct_cells')
    if (!dir.exists(dirpath)) {
        dir.create(dirpath, recursive=TRUE)
    }
    filepath = file.path(dirpath, paste0('pct_cells-', opt[['stat']], '.csv'))
    write.table(tmp, file = filepath, row.names = FALSE,  sep = ',' )
}


# Absolute Counts
if ('abs_count' %in% colnames(df)) {

    if (opt[['stat']]=='t_test') {
        abs_count_tbl <- apply_unpaired_t_test(
            df[!is.na(df[['abs_count']]), ],
            index_cols=c('organ', 'cell_type'),
            group_name='group_name',
            metric='abs_count'
        )
    } else {
        abs_count_tbl <- apply_multiple_comparisons(
            df[!is.na(df[['abs_count']]), ],
            index_cols=c('organ', 'cell_type'),
            group_name='group_name',
            metric='abs_count',
            correction=opt[['stat']]
        )
    }

    # save
    if (!troubleshooting) {
        tmp <- abs_count_tbl
        for (col in group_names) {
            tmp[[col]] <- mapply(toJSON, tmp[[col]])
        }

        dirpath <- file.path(wd, opt[['output-dir']], 'data',
            gsub(',', '_', opt[['group-by']]), 'abs_count')
        if (!dir.exists(dirpath)) {
            dir.create(dirpath, recursive=TRUE)
        }
        filepath = file.path(dirpath, paste0('abs_count-', opt[['stat']], '.csv'))
        write.table(tmp, file = filepath, row.names = FALSE,  sep = ',' )
    }
}


# ----------------------------------------------------------------------
# Data Overview

organs <- sort(unique(df[['organ']]))
log_print(paste(Sys.time(), 'Groups found...', paste(organs, collapse = ', ')))
log_print(paste(Sys.time(), 'Exporting data...'))

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
                   'num_cells', 'pct_cells', 'abs_count', 'viable_cells_conc',
                   'total_vol', 'total_viable_cells',
                   'mouse_id', 'sex', 'weeks_old', 'fcs_name'),
                colnames(df)
            ))],
        x='cell_type',
        y='pct_cells',
        groups=c('group_name')
    )

    if (!troubleshooting) {
        dirpath <- file.path(wd, opt[['output-dir']], 'data',
            gsub(',', '_', opt[['group-by']]), 'pct_cells')  # created above
        filepath = file.path(dirpath,
            paste0(organ, '.csv')
        )
        write.table(tmp, file = filepath, row.names = FALSE, sep = ',')
    }


    # ----------------------------------------------------------------------
    # Plot


    # ------------------------------------------------------------
    # Percent Cells

    fig <- plot_violin(
        tmp,
        x='cell_type', y='pct_cells', group_by='group_name',
        ylabel=paste('Percent of', last_initial_gate), title=organ,
        ymin=0, ymax=100,
        hover_data=unique(intersect(
                c('mouse_id', 'group_name', metadata_cols,
                  'sex', 'treatment', 'weeks_old',
                  'Cells', 'num_cells', 'fcs_name'),
                colnames(df)
        ))
    )

    if (!troubleshooting) {
        save_fig(
            fig=fig,
            height=opt[['height']], width=opt[['width']],
            dirpath=file.path(wd, opt[['output-dir']], 'figures',
                gsub(',', '_', opt[['group-by']]), 'pct_cells', 'interactive'),
            filename=paste(organ, 'pct_cells', sep='-'),
            save_html=TRUE
        )
    }

    # ------------------------------------------------------------
    # Absolute Counts

    if ('abs_count' %in% colnames(df)) {
        fig <- plot_violin(
            tmp,
            x='cell_type', y='abs_count', group_by='group_name',
            ylabel='Absolute Count', title=organ,
            ymin=0,
            hover_data=unique(intersect(
                    c('mouse_id', 'group_name', metadata_cols,
                      'sex', 'treatment', 'weeks_old', 
                      'Cells', 'num_cells', 'fcs_name'),
                    colnames(df)
            ))
        )

        if (!troubleshooting) {
            save_fig(
                fig=fig,
                height=opt[['height']], width=opt[['width']],
                dirpath=file.path(wd, opt[['output-dir']], 'figures',
                    gsub(',', '_', opt[['group-by']]), 'abs_count', 'interactive'),
                filename=paste(organ, 'abs_count', sep='-'),
                save_html=TRUE
            )
        }
    }
}


# ----------------------------------------------------------------------
# Multiple Comparisons

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
        unique(intersect(
                c('organ', 'cell_type', 'group_name', 'pct_cells', 'abs_count'),
                colnames(df)
        ))
    ]

    # ----------------------------------------------------------------------
    # Percent Cells

    final_group_order <- move_list_items_to_front(
        unique(df_subset[['group_name']]),
        custom_group_order
    )
    n_groups <- length(unique(final_group_order))
    
    withCallingHandlers({
        fig <- plot_multiple_comparisons(
            df_subset[, c("organ", "cell_type", "group_name", "pct_cells")],
            x='group_name', y='pct_cells',
            ymin=min( df_subset[['pct_cells']]*1.1, 0 ),
            ylabel='Percent of\nLive Cells',
            title=paste(toupper(organ), cell_type),
            test=opt[['stat']],
            show_numbers=opt[['show-numbers']],
            custom_group_order=final_group_order  # manual input
        )
    }, warning = function(w) {
        if ( any(grepl("containing non-finite values", w),
                 grepl("outside the scale range", w)) ) {
            invokeRestart("muffleWarning")
        }
    })

    # save
    if (!troubleshooting) {
        
        dirpath <- file.path(wd, opt[['output-dir']], 'figures',
            gsub(',', '_', opt[['group-by']]), 'pct_cells', organ,
            opt[['stat']]
        )
        if (!dir.exists(dirpath)) { dir.create(dirpath, recursive=TRUE) }
        
        filepath = file.path(dirpath,
            paste0(gsub(' ', '', tolower(organ)), '-',
                gsub(' ', '_', tolower(cell_type)), '.svg' )
        )
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

    # ----------------------------------------------------------------------
    # Absolute Counts

    if ('abs_count' %in% colnames(df)) {

        tmp <- df_subset[
            !is.na(df_subset[['abs_count']]), 
            c("organ", "cell_type", "group_name", "abs_count")
        ]

        if (nrow(tmp)>0) {
            fig <- plot_multiple_comparisons(
                tmp,
                x='group_name', y='abs_count',
                ymin=min( df_subset[['abs_count']]*1.1, 0 ),
                ylabel='Absolute Count',
                title=paste(toupper(organ), cell_type),
                test=opt[['stat']],
                show_numbers=opt[['show-numbers']],
                custom_group_order=final_group_order
            )

            # save
            if (!troubleshooting) {
                
                dirpath <- file.path(wd, opt[['output-dir']], 'figures',
                    gsub(',', '_', opt[['group-by']]), 'abs_count',
                    paste(organ, 'absolute', sep='-'),
                    opt[['stat']]
                )
                if (!dir.exists(dirpath)) { dir.create(dirpath, recursive=TRUE) }
                
                filepath = file.path(dirpath,
                    paste0(gsub(' ', '', tolower(organ)), '-',
                        gsub(' ', '_', tolower(cell_type)), '.svg' )
                )

                # need bugfix for compressed images when more than 2 groups
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
        }
    }

    pbar$tick()
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
