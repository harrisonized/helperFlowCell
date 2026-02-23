## Like analyze_mfi, but for data only

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('logr'))
library('optparse')

import::from(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'find_initial_gates', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'file_io.R'),
    'import_flowjo_export', 'import_flow_metadata', .character_only=TRUE)

import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'group_by_agg', 'rename_columns',
    .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'math.R'),
    'compute_normal_tvd', .character_only=TRUE)
import::from(file.path(wd, 'R', 'config', 'flow.R'),
    'mouse_db_ignore', .character_only=TRUE)


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
print(metadata_cols)

last_initial_gate <- 'Live Cells'

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
    last_initial_gate=NA
)
df <- df[, c('fcs_name', 'gate', 'cell_type', 'gmfi')]
if (is.null(df)) {
    stop(paste('No files found in', file.path(wd, opt[['input-dir']])))
}


# left join sdev
sdev_df <- import_flowjo_export(file.path(wd, opt[['sdev-dir']]), metric_name='sdev')
sdev_df <- sdev_df[, c('fcs_name', 'gate', 'cell_type', 'sdev')]
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
    last_initial_gate=last_initial_gate
)
if (!is.null(counts_df)) {
    log_print(paste(Sys.time(), 'Merging counts data...'))

    divisor <- tail(find_initial_gates(colnames(counts_df), last_initial_gate), 1)[[1]]
    counts_df[['pct_cells']] <- round(counts_df[['num_cells']] /
        counts_df[[divisor]] * 100, 4)

    counts_cols <- c("fcs_name", "gate", "cell_type",  "num_cells", "pct_cells")  # don't need pct_cells
    counts_df <- counts_df[, counts_cols]

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
# Total Variation Distance

if (!is.null(sdev_df)) {
    log_print(paste(Sys.time(), 'Computing total variation distances...'))

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

    if (!troubleshooting) {
        dirpath <- file.path(wd, opt[['output-dir']], 'data',
            gsub(',', '_', opt[['group-by']]), opt[['metric']])
        if (!dir.exists(dirpath)) { dir.create(dirpath, recursive=TRUE) }
        write.table(tvd_df, file=file.path(dirpath, paste0('gmfi-tvd', '.csv')),
            row.names=FALSE, sep=',' )
    }
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
