## Graphs percentage of Live Cells for each cell population
## 
## To plot absolute counts, include abs_count as a column in the metadata file
## If your initial gate is not found, this script will falls back on the top-level gate as the divisor

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('logr'))
library('optparse')
import::from(progress, 'progress_bar')
import::from(jsonlite, 'toJSON')
import::from(ggplot2, 'ggsave')

import::from(file.path(wd, 'R', 'functions', 'file_io.R'),
    'import_flowjo_export_v2', 'import_flow_metadata', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', 'move_list_items_to_front',
    .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'math.R'),
    'apply_unpaired_t_test', 'apply_multiple_comparisons', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'save_fig', 'plot_violin', 'plot_ribbon_chart', 'plot_multiple_comparisons',
    .character_only=TRUE)
import::from(file.path(wd, 'R', 'config', 'flow.R'),
    'mouse_db_ignore', .character_only=TRUE)
import::from(file.path(wd, 'R', 'config', 'user_input.R'),
    'min_cells', 'cell_type_to_color', 'composite_gates', 'custom_group_order',
    .character_only=TRUE)

reticulate::py_config()  # required on linux to access reticulate


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/flow-counts',
                metavar='data/flow-counts', type="character",
                help="input directory, all csv files will be read in"),

    make_option(c("-m", "--metadata-dir"), default="data/flow-metadata",
                metavar="data/flow-metadata", type="character",
                help="metadata directory containing csvs mapping fcs to mouse ids"),

    make_option(c("-o", "--output-dir"), default="hfc-output",
                metavar="hfc-output", type="character",
                help="output directory for the data"),

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

    make_option(c("-d", "--divisor"), default="Live Cells",
                metavar="Live Cells", type="character",
                help="Can also choose: 'Single Cells', 'CD45+', etc."),

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
divisor <- opt[['divisor']]


# Start Log
start_time = Sys.time()
log <- log_open(paste0("analyze_counts-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read Data

log_print(paste(Sys.time(), 'Reading data...'))

# Required files
df <- import_flowjo_export_v2(file.path(wd, opt[['input-dir']]), metric_name='num_cells', divisor=divisor)
if (is.null(df)) { stop(paste('No files found in', file.path(wd, opt[['input-dir']])))}
flow_metadata <- import_flow_metadata(file.path(wd, opt[['metadata-dir']]))
if (is.null(flow_metadata)) { stop(paste('No files found in', file.path(wd, opt[['metadata-dir']])))}

# Optional files
mouse_db <- append_many_csv(file.path(wd, opt[['ref-dir']]), recursive=TRUE, include_filepath=FALSE)


# ----------------------------------------------------------------------
# Preprocessing

log_print(paste(Sys.time(), 'Preprocessing...'))

df[['pct_cells']] <- round(df[['num_cells']] / df[["divisor"]] * 100, 4)

# inner join flow metadata
df <- merge(df, flow_metadata,
    by="fcs_name", all.x=FALSE, all.y=FALSE, suffixes=c('', '_'))
if ('is_unstained' %in% colnames(df)) {
    df <- df[(df[['is_unstained']]==FALSE), ]  # backwards compatibility
}
if ('is_excluded' %in% colnames(df)) {
    df <- df[(df[['is_excluded']]==FALSE), ]  # drop unstained cells
}

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
df <- df[(df[['num_cells']] > min_cells), ]    # filter
df <- df[rowSums(is.na(df)) != ncol(df), ]  # filter NA rows caused by the above line

# exclude negative gates
cell_types <- unique(df[['cell_type']])
negative_gates <- cell_types[grepl("-$", cell_types)]
df <- df[((df[['cell_type']] %in% negative_gates)==FALSE), ]


# ----------------------------------------------------------------------
# Overview Box Plot

organs <- sort(unique(df[['organ']]))
log_print(paste(Sys.time(), 'Groups found...', paste(organs, collapse = ', ')))

for (organ in sort(organs)) {

    log_print(paste0(Sys.time(), ' Processing ', organ, '...'))

    # save data
    if (!troubleshooting) {
        dirpath <- file.path(wd, opt[['output-dir']], 'data',
            gsub(',', '_', opt[['group-by']]), 'counts')
        filepath = file.path(dirpath, paste0('counts-', organ, '.csv'))
        if (!dir.exists(dirpath)) { dir.create(dirpath, recursive=TRUE) }
        write.table(df[(df[['organ']]==organ), ], file = filepath, row.names = FALSE, sep = ',')
    }

    # ------------------------------------------------------------
    # Percent Cells

    fig <- plot_violin(
        df[(df[['organ']]==organ), ],
        x='cell_type', y='pct_cells', group_by='group_name',
        ylabel=paste('Percent of', divisor), title=organ,
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
            height=500, width=2000,
            dirpath=file.path(wd, opt[['output-dir']], 'figures',
                gsub(',', '_', opt[['group-by']]),
                'overview', 'pct_cells'),
            filename=paste(organ, 'pct_cells', sep='-'),
            save_html=TRUE
        )
    }

    # ------------------------------------------------------------
    # Absolute Counts

    if ('abs_count' %in% colnames(df)) {
        fig <- plot_violin(
            df[(df[['organ']]==organ), ],
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
                height=500, width=2000,
                dirpath=file.path(wd, opt[['output-dir']], 'figures',
                    gsub(',', '_', opt[['group-by']]),
                    'overview', 'abs_count'),
                filename=paste(organ, 'abs_count', sep='-'),
                save_html=TRUE
            )
        }
    }
}


# ----------------------------------------------------------------------
# Ribbon Charts

log_print(paste0(Sys.time(), ' Plotting ribbon charts...'))

ribbon_df <- df %>%
    group_by(across(all_of(c("organ", "cell_type", metadata_cols)))) %>%
    summarise(
        mean_pct_cells = mean(.data[["pct_cells"]], na.rm = TRUE),
        !!!if ("abs_count" %in% names(df)) {
            list(mean_abs_count = rlang::expr(mean(.data[["abs_count"]], na.rm = TRUE)))
        } else {
            list()
        },
        .groups = "drop"
    )
ribbon_df <- as.data.frame(ribbon_df)
ribbon_df[['group_name']] <- apply( ribbon_df[ , metadata_cols ] , 1 , paste , collapse = ", " )
ribbon_df <- ribbon_df[!ribbon_df[['cell_type']] %in% composite_gates, ]

for (organ in organs) {

    # mean_pct_cells
    fig <- plot_ribbon_chart(
        df = ribbon_df[(ribbon_df[['organ']]==organ), ],
        x = "group_name", y = "mean_pct_cells", subcat = "cell_type",
        title = organ,
        bar_width = 0.7,
        normalize = TRUE,
        show_values = TRUE,
        show_axes = TRUE,
        show_grid = TRUE,
        show_xlabels = TRUE,
        color=cell_type_to_color,
        group_order = sort(unique(df[['group_name']]))
    )

    if (!troubleshooting) {
        save_fig(
            fig=fig,
            height=400, width=600, scale=4,
            dirpath=file.path(wd, opt[['output-dir']], 'figures',
                gsub(',', '_', opt[['group-by']]),
                'ribbon_charts', 'pct_cells'),
            filename=paste(organ, 'ribbon', sep='-'),
            save_html=TRUE
        )
    }

    # mean_abs_count
    if ('mean_abs_count' %in% colnames(ribbon_df)) {
        fig <- plot_ribbon_chart(
            df = ribbon_df[(ribbon_df[['organ']]==organ), ],
            x = "group_name", y = "mean_abs_count", subcat = "cell_type",
            title = organ,
            bar_width = 0.7,
            normalize = FALSE,
            show_values = TRUE,
            show_axes = TRUE,
            show_grid = TRUE,
            show_xlabels = TRUE,
            colors = cell_type_to_color,
            group_order = sort(unique(df[['group_name']]))
        )

        if (!troubleshooting) {
            save_fig(
                fig=fig,
                height=400, width=600, scale=4,
                dirpath=file.path(wd, opt[['output-dir']], 'figures',
                    gsub(',', '_', opt[['group-by']]),
                    'ribbon_charts', 'abs_counts'),
                filename=paste(organ, 'ribbon', sep='-'),
                save_html=TRUE
            )
        }
    }

}


# ----------------------------------------------------------------------
# Multiple Comparisons

if (n_distinct(df[['group_name']])==1) {

    log_print(paste(Sys.time(), paste0('Single group found. Skipping multiple comparisons...')))

} else {

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
            gsub(',', '_', opt[['group-by']]), 'statistics')
        if (!dir.exists(dirpath)) { dir.create(dirpath, recursive=TRUE) }
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
                gsub(',', '_', opt[['group-by']]), 'statistics')
            if (!dir.exists(dirpath)) { dir.create(dirpath, recursive=TRUE) }
            filepath = file.path(dirpath, paste0('abs_count-', opt[['stat']], '.csv'))
            write.table(tmp, file = filepath, row.names = FALSE,  sep = ',' )
        }
    }

    # ----------------------------------------------------------------------
    # Plots

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
                ylabel=paste('Percent of', divisor),
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
                gsub(',', '_', opt[['group-by']]),
                'multiple_comparisons', 'pct_cells', organ, opt[['stat']]
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
                        gsub(',', '_', opt[['group-by']]),
                        'multiple_comparisons', 'abs_count', organ, opt[['stat']]
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
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
