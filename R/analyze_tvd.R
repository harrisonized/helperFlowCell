## Graphs TVD for each cell population
## 

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('logr'))
library('optparse')
library('zeallot')  # %<-%
import::from(progress, 'progress_bar')

import::from(file.path(wd, 'R', 'functions', 'file_io.R'),
    'import_flowjo_export_v2', 'import_flow_metadata', .character_only=TRUE)

import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'fillna', 'group_by_agg', 'rename_columns',
    .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'math.R'),
    'compute_normal_tvd', # 'compute_lognormal_tvd',
    .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'save_fig', 'plot_violin', # 'plot_multiple_comparisons',
    'plot_polar_area', .character_only=TRUE)
import::from(file.path(wd, 'R', 'config', 'flow.R'),
    'mouse_db_ignore', .character_only=TRUE)
import::from(file.path(wd, 'R', 'config', 'user_input.R'),
    'composite_gates', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dirs"), default='data/flow-counts,data/flow-mfi,data/flow-sdev',
                metavar='data/flow-counts,data/flow-mfi,data/flow-sdev', type="character",
                help="input directories, enter a comma-separated list in the same order as the example"),

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

    make_option(c("-f", "--fluorescence"), default="zygosity",
                metavar="zygosity", type="character",
                help="use this to split the fluorescence column from the group by"),

    make_option(c("-s", "--stat"), default='fishers_lsd',
                metavar='fishers_lsd', type="character",
                help="Choose 'fishers_lsd', 't_test', 'tukey', or 'bonferroni'"),

    make_option(c("-n", "--show-numbers"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="in violin with significance, show all numbers instead of stars"),

    make_option(c("-p", "--plot-histograms"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="generate histograms"),

    make_option(c("-x", "--xlabel"), default="Comp-GFP-A :: mNeonGreen",
                metavar="Comp-GFP-A :: mNeonGreen", type="character",
                help="x-axis label for histograms"),

    make_option(c("-y", "--metric"), default="gmfi",
                metavar="gmfi", type="character",
                help="metric name, eg. opt[['metric']] or 'mode'"),

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
c(counts_dir, mfi_dir, sdev_dir) %<-% unlist(strsplit(opt[['input-dirs']], ','))
metadata_cols <- unlist(strsplit(opt[['group-by']], ','))

divisor <- opt[['divisor']]


# Start Log
start_time = Sys.time()
log <- log_open(paste0("analyze_mfi-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read Data

log_print(paste(Sys.time(), 'Reading data...'))

# Required files
df <- import_flowjo_export_v2(file.path(wd, mfi_dir), metric_name=opt[['metric']], divisor=NULL)
if (is.null(df)) { stop(paste('No files found in', file.path(wd, dir)))}
flow_metadata <- import_flow_metadata(file.path(wd, opt[['metadata-dir']]))
if (is.null(flow_metadata)) { stop(paste('No files found in', file.path(wd, opt[['metadata-dir']])))}

# Optional files
counts_df <- import_flowjo_export_v2(file.path(wd, counts_dir), metric_name='num_cells', divisor=divisor)
if (is.null(counts_df)) { warning( paste('No files found in', file.path(wd, counts_dir)) )}
sdev_df <- import_flowjo_export_v2(file.path(wd, sdev_dir), metric_name='sdev', divisor=NULL)
if (is.null(sdev_df)) { warning( paste('No files found in', file.path(wd, sdev_dir)) )}
mouse_db <- append_many_csv(file.path(wd, opt[['ref-dir']]), recursive=TRUE, include_filepath=FALSE)
if (is.null(mouse_db)) { warning( paste('No files found in', file.path(wd, opt[['ref-dir']])) )}


# ----------------------------------------------------------------------
# Preprocessing

log_print(paste(Sys.time(), 'Preprocessing...'))

# Left join datasets
if (!is.null(sdev_df)) {
    log_print(paste(Sys.time(), 'Merging sdev data...'))
    df <- merge(df, sdev_df,
        by=c("fcs_name", "gate", "cell_type"),
        all.x=TRUE, all.y=FALSE, suffixes=c('', '_'))
    df[(is.na(df[['sdev']])), 'sdev'] <- 0  # fillna
}
if (!is.null(counts_df)) {
    log_print(paste(Sys.time(), 'Merging counts data...'))

    counts_df[['pct_cells']] <- round(counts_df[['num_cells']] / counts_df[['divisor']] * 100, 4)
    counts_cols <- c("fcs_name", "gate", "cell_type",  "num_cells", "pct_cells")  # remove pregates
    counts_df <- counts_df[, counts_cols]

    df <- merge(df, counts_df,
        by=c("fcs_name", "gate", "cell_type"),
        all.x=TRUE, all.y=FALSE, suffixes=c('', '_'))
}

# inner join flow metadata
df <- merge(df, flow_metadata, by="fcs_name", all.x=FALSE, all.y=FALSE, suffixes=c('', '_'))
if (length(metadata_cols) > 1) {
    df[['group_name']] <- apply( df[ , metadata_cols ] , 1 , paste , collapse = ", " )
    df[['subgroup_name']] <- apply(
        df[ , items_in_a_not_b(metadata_cols, opt[['fluorescence']]) , drop=FALSE],
        1 , paste , collapse = ", " )
} else {
    df[['group_name']] <- df[[metadata_cols]]
}

# left join mouse data
if (!is.null(mouse_db)) {
    log_print(paste(Sys.time(), 'Merging mouse data...'))
    df <- merge(df, mouse_db[, items_in_a_not_b(colnames(mouse_db), mouse_db_ignore)],
        by="mouse_id", all.x=TRUE, all.y=FALSE, suffixes=c('', '_'))
    df[['weeks_old']] <- round(df[['age']]/7, 1)
}

# filter and sort
df <- df[(df[['is_unstained']]==FALSE), ]  # drop unstained cells
if (!is.null(counts_df)) {
    df <- df[(df[['num_cells']] > 50), ]  # filter
}
df <- df[order(df[['cell_type']], df[['organ']], df[['group_name']]), ]  # sort rows


# ----------------------------------------------------------------------
# Compute Total Variation Distance

log_print(paste(Sys.time(), 'Computing TVDs...'))

group_cols <- c('organ', 'gate', 'cell_type', metadata_cols)
val_cols <- c('gmfi', 'sdev')
tvd_df <- df[
    ((df[opt[['fluorescence']]]!='WT') &
     (df[['sdev']] > 0) & (df[['num_cells']] > 5)),
    c('mouse_id', 'fcs_name', group_cols, 'num_cells', 'pct_cells', val_cols)]

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
      'num_cells', 'pct_cells', val_cols, 'wt_gmfi', 'wt_sdev')
]
tvd_df[['group_name']] <- apply(
    tvd_df[ , metadata_cols , drop = FALSE] , 1 , paste , collapse = ", "
)

tvd_df[['tvd']] <- mapply(
    compute_normal_tvd,
    mean1 = tvd_df[['gmfi']],
    sd1   = tvd_df[['sdev']],
    mean2 = tvd_df[['wt_gmfi']],
    sd2   = tvd_df[['wt_sdev']],
    log_transform=TRUE
)


# ------------------------------------------------------------
# Plotly Interactive Overview

organs <- sort(unique(tvd_df[['organ']]))
for (organ in sort(organs)) {

    if (!troubleshooting) {
        dirpath <- file.path(wd, opt[['output-dir']], 'data',
            gsub(',', '_', opt[['group-by']]), opt[['metric']])
        filepath = file.path(dirpath, paste0("tvd-", organ, '.csv'))
        if (!dir.exists(dirpath)) { dir.create(dirpath, recursive=TRUE) }
        write.table(tvd_df[(tvd_df[['organ']]==organ), ], file=filepath, row.names=FALSE, sep=',')
    }

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
            gsub(',', '_', opt[['group-by']]), 'tvd', 'interactive')
        if (!dir.exists(dirpath)) { dir.create(dirpath, recursive=TRUE) }
        save_fig(
            fig=fig,
            height=500, width=2000,
            dirpath=dirpath,
            filename=paste(organ, 'tvd', sep='-'),
            save_html=TRUE
        )
    }
}


# ----------------------------------------------------------------------
# Polar Area Charts

# group by
polar_df <- tvd_df[(tvd_df['zygosity'] == 'homo'), ] %>%
    group_by(organ, cell_type, treatment) %>%
    summarise(
        mean_pct_cells = mean(pct_cells, na.rm = TRUE),
        mean_tvd = mean(tvd, na.rm = TRUE),
        .groups = 'drop'
    )
polar_df <- polar_df[!polar_df[['cell_type']] %in% composite_gates, ]
polar_df <- as.data.frame(polar_df)
polar_df <- fillna(polar_df, 'mean_tvd', 0)

organs <- sort(unique(polar_df[['organ']]))
treatments <- unique(polar_df[['treatment']])
pbar <- progress_bar$new(total = length(organs))
for (organ in sort(organs)) {
    for (treatment in treatments) {
        fig <- plot_polar_area(
            polar_df[(polar_df[['organ']]==organ) & (polar_df[['treatment']]==treatment), ],
            x="mean_pct_cells", y="mean_tvd", group='cell_type',
            min_label_pct=100,
            colors=colors,
            title=organ,
            show_axes=TRUE
        )

        if (!troubleshooting) {
            save_fig(
                fig=fig,
                height=opt[['height']], width=opt[['width']],
                dirpath=file.path(wd, opt[['output-dir']], 'figures', 'polar_area'),
                filename=paste(organ, 'polar_area', treatment, sep='-'),
                save_html=TRUE
            )
        }
    }

    pbar$tick()
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
