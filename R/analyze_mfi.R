## Graphs MFI for each cell population
## There are three different ways of calculating MFI
## 1. Geometric mean: This is the best choice, because the calculation compensates
## for extreme skew, which is standard in flow cytometry data due to the logarithmic
## nature of flow cytometry data.
## 2. Median is a good second choice, however, in situations where the distribution
## is an exponential decay with 0 being the most common value, the median will not
## accurately represent the distribution
## 3. The Arithmetic Mean (ie. the Mean) is a bad choice. It will systemically
## overestimate the representative value due to the logarithmic nature of flow
## cytometry data.
## This script treats all three kinds of MFI equally, however, you can rename it on
## graph using the -y option

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
    'sort_groups_by_metric', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'file_io.R'),
    'import_flowjo_export', 'import_flow_metadata', .character_only=TRUE)

import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'group_by_agg', 'rename_columns',
    .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', 'multiple_replacement', 'move_list_items_to_front',
    .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'math.R'),
    'apply_unpaired_t_test', 'apply_multiple_comparisons',
    'generate_lognormal_data', 'compute_normal_tvd',
    .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'save_fig', 'plot_violin', 'plot_multiple_comparisons',
    'plot_modal_histograms', .character_only=TRUE)
import::from(file.path(wd, 'R', 'config', 'flow.R'),
    'mouse_db_ignore', .character_only=TRUE)
import::from(file.path(wd, 'R', 'config', 'user_input.R'),
    'custom_group_order', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/flow-mfi',
                metavar='data/flow-gmfi', type="character",
                help="input directory, all csv files will be read in"),

    make_option(c("-o", "--output-dir"), default="hfc-output",
                metavar="hfc-output", type="character",
                help="output directory for the data"),

    make_option(c("-d", "--sdev-dir"), default='data/flow-sdev',
                metavar='data/flow-sdev', type="character",
                help="standard deviations"),

    make_option(c("-c", "--counts-dir"), default='data/flow-counts',
                metavar='data/flow-counts', type="character",
                help="counts"),

    make_option(c("-m", "--metadata-dir"), default="data/flow-metadata",
                metavar="data/flow-metadata", type="character",
                help="metadata directory containing csvs mapping fcs to mouse ids"),

    make_option(c("-r", "--ref-dir"), default='data/mice',
                metavar='data/mice', type="character",
                help="directory of files containing all the mouse data"),

    make_option(c("-g", "--group-by"), default='sex,treatment,zygosity',
                metavar='treatment', type="character",
                help="enter a column or comma-separated list of columns, no spaces"),

    make_option(c("-p", "--plot-histograms"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="generate histograms"),

    make_option(c("-y", "--metric"), default="gmfi",
                metavar="gmfi", type="character",
                help="metric name, eg. opt[['metric']] or 'mode'"),

    make_option(c("-s", "--stat"), default='fishers_lsd',
                metavar='fishers_lsd', type="character",
                help="Choose 'fishers_lsd', 't_test', 'tukey', or 'bonferroni'"),

    make_option(c("-n", "--show-numbers"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="in violin with significance, show all numbers instead of stars"),

    make_option(c("-f", "--fluorescence"), default="zygosity",
                metavar="zygosity", type="character",
                help="use this to split out the fluorescence column from the group by"),

    make_option(c("-x", "--xlabel"), default="Comp-GFP-A :: mNeonGreen",
                metavar="Comp-GFP-A :: mNeonGreen", type="character",
                help="x-axis label for histograms"),

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
log <- log_open(paste0("analyze_mfi-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read Data

log_print(paste(Sys.time(), 'Reading data...'))


# import mfi
df <- import_flowjo_export(
    file.path(wd, opt[['input-dir']]),
    metric_name=opt[['metric']],
    include_initial_gates=FALSE
)
if (is.null(df)) {
    stop(paste('No files found in', file.path(wd, opt[['input-dir']])))
}


# left join sdev
sdev_df <- import_flowjo_export(file.path(wd, opt[['sdev-dir']]), metric_name='sdev')
if (!is.null(sdev_df)) {
    log_print(paste(Sys.time(), 'Merging sdev data...'))
    df <- merge(df, sdev_df,
        by=c("fcs_name", "gate", "cell_type"),
        all.x=TRUE, all.y=FALSE, suffixes=c('', '_'))
    df[(is.na(df[['sdev']])), 'sdev'] <- 0  # fillna
}


# left join counts
counts_df <- import_flowjo_export(
    file.path(wd, opt[['counts-dir']]),
    metric_name='num_cells',
    include_initial_gates=TRUE
)
if (!is.null(counts_df)) {
    log_print(paste(Sys.time(), 'Merging counts data...'))

    counts_cols <- c("fcs_name", "gate", "cell_type",  "num_cells")  # don't need pct_cells
    counts_df <- counts_df[, counts_cols]

    ss_df <- import_flowjo_export(
        file.path(wd, opt[['counts-dir']]),
        metric_name='num_cells',
        include_initial_gates=FALSE
    )
    ss_df <- ss_df[
        ((ss_df[['gate']]=='Cells/Single Cells/Single Cells') |
         (ss_df[['gate']]=='Cells/Single Cells/Single Cells/Live Cells')),
        counts_cols
    ]
    counts_df <- rbind(counts_df, ss_df)

    df <- merge(df, counts_df,
        by=c("fcs_name", "gate", "cell_type"),
        all.x=TRUE, all.y=FALSE, suffixes=c('', '_'))
}


# inner join flow metadata
flow_metadata <- import_flow_metadata(file.path(wd, opt[['metadata-dir']]))
df <- merge(df, flow_metadata,
    by="fcs_name", all.x=FALSE, all.y=FALSE, suffixes=c('', '_'))
df <- df[(df[['is_unstained']]==FALSE), ]  # drop unstained cells
if (length(metadata_cols) > 1) {
    df[['group_name']] <- apply( df[ , metadata_cols ] , 1 , paste , collapse = ", " )
    df[['subgroup_name']] <- apply(
        df[ , items_in_a_not_b(metadata_cols, opt[['fluorescence']]) , drop=FALSE],
        1 , paste , collapse = ", " )
} else {
    df[['group_name']] <- df[[metadata_cols]]
}


# left join mouse data
mouse_db <- append_many_csv(
    file.path(wd, opt[['ref-dir']]), recursive=TRUE, include_filepath=FALSE)
if (!is.null(mouse_db)) {
    log_print(paste(Sys.time(), 'Merging mouse data...'))
    df <- merge(df, mouse_db[, items_in_a_not_b(colnames(mouse_db), mouse_db_ignore)],
        by="mouse_id", all.x=TRUE, all.y=FALSE, suffixes=c('', '_'))
    df[['weeks_old']] <- round(df[['age']]/7, 1)
}


# filter and sort
df <- df[order(df[['cell_type']], df[['organ']], df[['group_name']]), ]  # sort rows


# ----------------------------------------------------------------------
# Data Overview

organs <- sort(unique(df[['organ']]))
log_print(paste(Sys.time(), 'Organs included...', paste(organs, collapse = ', ')))

for (organ in sort(organs)) {

    log_print(paste(Sys.time(), 'Plotting', organ, 'overview...'))

    # ------------------------------------------------------------
    # Data

    tmp <- sort_groups_by_metric(
        df[(df[['organ']]==organ),
            unique(intersect(
                c('organ', 'genotype', 'treatment',
                   'group_name', metadata_cols, 'cell_type',
                   'Cells/Single Cells/Single Cells/Live Cells',
                   'num_cells', opt[['metric']], 'viable_cells_conc',
                   'total_vol', 'total_viable_cells', 'abs_count',
                   'mouse_id', 'sex', 'weeks_old', 'fcs_name'),
                colnames(df)
            ))],
        x='cell_type',
        y=opt[['metric']],
        groups=c('group_name')
    )

    if (!troubleshooting) {
        dirpath <- file.path(wd, opt[['output-dir']], 'data',
            gsub(',', '_', opt[['group-by']]), opt[['metric']])
        if (!dir.exists(dirpath)) { dir.create(dirpath, recursive=TRUE) }
        filepath = file.path(dirpath, paste0(opt[['metric']], '-', organ, '.csv'))
        write.table(tmp, file = filepath, row.names = FALSE, sep = ',')
    }

    # ------------------------------------------------------------
    # Plot

    fig <- plot_violin(
        tmp,
        x='cell_type', y=opt[['metric']], group_by='group_name',
        ylabel='MFI', title=organ,
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
                gsub(',', '_', opt[['group-by']]), opt[['metric']], 'interactive'),
            filename=paste(organ, opt[['metric']], sep='-'),
            save_html=TRUE
        )
    }

}


# ----------------------------------------------------------------------
# Total Variation Distance

if (!is.null(sdev_df)) {

    # ------------------------------------------------------------
    # Data

    log_print(paste(Sys.time(), 'Computing total variation distances...'))

    group_cols <- c('organ', 'gate', 'cell_type', metadata_cols)
    val_cols <- c('gmfi', 'sdev')
    tvd_df <- df[
        ((df[opt[['fluorescence']]]!='WT') &
         (df[['sdev']] > 0) & (df[['num_cells']] > 5)),
        c('mouse_id', 'fcs_name', group_cols, 'num_cells', val_cols)]

    # left join WT mean MFI
    wt_df <- group_by_agg(
        df[((df[opt[['fluorescence']]]=='WT') &
            (df[['sdev']] > 0) & (df[['num_cells']] > 5)),
           c(group_cols, val_cols)],
        group_cols, val_cols, mean
    )
    wt_df <- rename_columns(wt_df, c('gmfi'='wt_gmfi', 'sdev'='wt_sdev'))
    
    tvd_df <- merge(
        tvd_df, wt_df,
        by=items_in_a_not_b(group_cols, opt[['fluorescence']]),
        all.x=TRUE, all.y=FALSE, suffixes=c('', '_'))[,
        c('mouse_id', 'fcs_name', group_cols,
          'num_cells', val_cols, 'wt_gmfi', 'wt_sdev')
    ]
    tvd_df[['group_name']] <- apply(
        tvd_df[ , metadata_cols ] , 1 , paste , collapse = ", " )

    tvd_df[['tvd']] <- mapply(
        compute_normal_tvd,
        mean1 = tvd_df[['gmfi']],
        sd1   = tvd_df[['sdev']],
        mean2 = tvd_df[['wt_gmfi']],
        sd2   = tvd_df[['wt_sdev']],
        log_transform=TRUE
    )

    if (!troubleshooting) {
        dirpath <- file.path(wd, opt[['output-dir']], 'data',
            gsub(',', '_', opt[['group-by']]), opt[['metric']])
        # if (!dir.exists(dirpath)) { dir.create(dirpath, recursive=TRUE) }
        write.table(tvd_df, file=file.path(dirpath, paste0('gmfi-tvd', '.csv')),
            row.names=FALSE, sep=',' )
    }

    # ------------------------------------------------------------
    # Plot

    organs <- sort(unique(tvd_df[['organ']]))
    for (organ in sort(organs)) {

        log_print(paste(Sys.time(), 'Plotting', organ, 'tvd...'))

        fig <- plot_violin(
            tvd_df[(tvd_df[['organ']]==organ), ],
            x='cell_type', y='tvd', group_by='group_name',
            ylabel='TVD', title=organ,
            ymin=0, ymax=1,
            hover_data=unique(intersect(
                    c('mouse_id', 'group_name', metadata_cols,
                      'sex', 'treatment', 'weeks_old', 
                      'Cells', 'num_cells', 'fcs_name'),
                    colnames(df)
            ))
        )

        if (!troubleshooting) {
            dirpath <- file.path(wd, opt[['output-dir']], 'figures',
                gsub(',', '_', opt[['group-by']]), opt[['metric']], 'interactive')
            # if (!dir.exists(dirpath)) { dir.create(dirpath, recursive=TRUE) }
            save_fig(
                fig=fig,
                height=opt[['height']], width=opt[['width']],
                dirpath=dirpath,
                filename=paste(organ, 'tvd', sep='-'),
                save_html=TRUE
            )
        }
    }
}


# ----------------------------------------------------------------------
# MFI Multiple Comparisons

log_print(paste(Sys.time(), paste0('Computing ', opt[['stat']], '...')))

if (opt[['stat']]=='t_test') {
    pval_tbl <- apply_unpaired_t_test(
        df,
        index_cols=c('organ', 'cell_type'),
        group_name='group_name',
        metric=opt[['metric']]
    )
} else {
    pval_tbl <- apply_multiple_comparisons(
        df,
        index_cols=c('organ', 'cell_type'),
        group_name='group_name',
        metric=opt[['metric']],
        correction=opt[['stat']]
    )
}

log_print(paste(Sys.time(), paste(
    nrow(pval_tbl), 'populations found across',
    length(organs), 'organs') )
)

if (!troubleshooting) {
    tmp <- pval_tbl
    group_names <- sort(unique( df[['group_name']] ))
    for (col in group_names) {
        tmp[[col]] <- mapply(toJSON, tmp[[col]])
    }

    dirpath <- file.path(wd, opt[['output-dir']], 'data',
        gsub(',', '_', opt[['group-by']]), opt[['metric']])
    # if (!dir.exists(dirpath)) { dir.create(dirpath, recursive=TRUE) }
    filepath = file.path(dirpath, paste0('gmfi-pvals-', opt[['stat']], '.csv'))
    write.table(tmp, file = filepath, row.names = FALSE,  sep = ',' )
}


log_print(paste(Sys.time(), paste('Plotting mulitple comparisons using', opt[['stat']])))
if (opt[['plot-histograms']] && !is.null(sdev_df)) {
    log_print(paste(Sys.time(), 'Plotting fluorescence histograms...'))
}

pbar <- progress_bar$new(total = nrow(pval_tbl))
for (idx in 1:nrow(pval_tbl)) {

    # subset
    organ <- pval_tbl[idx, 'organ'][[1]]
    cell_type <- pval_tbl[idx, 'cell_type'][[1]]
    df_subset <- df[
        (df[['organ']] == organ) & (df[['cell_type']] == cell_type),
        unique(intersect(
                c('mouse_id', 'organ', 'cell_type', 'group_name',
                  'subgroup_name', opt[['fluorescence']],
                   opt[['metric']], 'sdev', 'num_cells', 'abs_count'),
                colnames(df)
        ))
    ]

    # ----------------------------------------------------------------------
    # MFI Violin Plot

    custom_group_order <- move_list_items_to_front(
        unique(df_subset[['group_name']]),
        custom_group_order
    )
    n_groups <- length(unique(custom_group_order))

    withCallingHandlers({
        fig <- plot_multiple_comparisons(
            df_subset[, c("organ", "cell_type", "group_name", opt[['metric']])],
            x='group_name', y=opt[['metric']],
            ymin=min( df_subset[, opt[['metric']]]*1.1, 0 ),
            ylabel=opt[['metric']],
            title=paste(toupper(organ), cell_type),
            test=opt[['stat']],
            show_numbers=opt[['show-numbers']],
            custom_group_order=custom_group_order  # manual input
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
            gsub(',', '_', opt[['group-by']]), opt[['metric']], organ,
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
            if ( any(grepl("containing non-finite values", w),
                     grepl("outside the scale range", w),
                     grepl("fewer than two data points", w),
                     grepl("argument is not numeric or logical: returning NA", w)) ) {
                invokeRestart("muffleWarning")
            }
        })
    }

    # ----------------------------------------------------------------------
    # Histograms

    if (opt[['plot-histograms']] && !is.null(sdev_df)) {

        subgroups <- unique(df_subset[['subgroup_name']])
        for (subgroup in subgroups) {

            df_subsubset <- df_subset[
                (df_subset['subgroup_name'] == subgroup),
                unique(intersect(
                        c('organ', 'mouse_id', 'group_name', 'subgroup_name', 'zygosity',
                          'num_cells', 'gmfi', 'sdev'),
                        colnames(df_subset)
                ))
            ]
            df_subsubset <- df_subsubset [!duplicated(df_subsubset[c('mouse_id', 'group_name')]), ]
            df_subsubset[['color']] <- multiple_replacement(
                df_subsubset[['zygosity']],
                c('homo'='#ACF53B', 'het'='#ACF53B', 'hemi'='#ACF53B', 'WT'='#B0B0B0')
            )

            sim_fluors <- Map(
                generate_lognormal_data,
                n = min(df_subsubset[['num_cells']], 1000, na.rm=TRUE),
                mean = df_subsubset[['gmfi']],
                sd = df_subsubset[['sdev']],
                group_name = apply(
                    df_subsubset[ , c('group_name', 'mouse_id') ],
                    1 , paste , collapse = ", " )
            )
            sim_fluor <- do.call(rbind, sim_fluors)

            fig <- plot_modal_histograms(sim_fluor,
                xlabel = opt[['xlabel']],
                ylabel = "Normalized to Mode",
                title = paste(toupper(organ), cell_type, subgroup),
                spar = 0.4,
                colors = df_subsubset[['color']],  # autogenerate this
                x_ticks = c(0, 250, 500, 1000, 2500, 5000, 10000),
                show_bins=FALSE
            )

            # save
            if (!troubleshooting) {
                dirpath <- file.path(wd, opt[['output-dir']], 'figures',
                    gsub(',', '_', opt[['group-by']]), opt[['metric']],
                    organ, 'histograms'
                )
                if (!dir.exists(dirpath)) { dir.create(dirpath, recursive=TRUE) }
                
                filepath = file.path(dirpath,
                    paste0(gsub(' ', '', tolower(organ)), '-',
                        gsub(' ', '_', tolower(cell_type)), '-',
                        gsub(', ', '_', subgroup), '.png' )
                )
                ggsave(
                    filepath,
                    plot=fig, bg='white',
                    height=1500, width=3000, dpi=200, units='px', scaling=2
                )
            }
        }
    }

    pbar$tick()
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
