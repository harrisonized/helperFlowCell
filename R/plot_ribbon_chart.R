## Creates a ribbon chart to compare two groups
## Can be used for KO vs. WT comparisons

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('logr'))
library('optparse')
import::from(progress, 'progress_bar')
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'save_fig', 'plot_ribbon_chart', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'group_by_agg', 'fillna', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings


# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='database/ko/steady-state',
                metavar='database/ko/steady-state', type="character",
                help="input directory, all csv files will be read in"),

    make_option(c("-o", "--output-dir"), default="hfc-output",
                metavar="hfc-output", type="character",
                help="output directory"),

    make_option(c("-l", "--height"), default=500,
                metavar="500", type="integer",
                help="height in px"),

    make_option(c("-w", "--width"), default=800,
                metavar="800", type="integer",
                help="width in px"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]
output_dir <- file.path(wd, opt[['output-dir']])
troubleshooting_dir = file.path(output_dir, 'troubleshooting')


# filter
plotly_colors <- c(
    '#1f77b4',  # blue
    '#ff7f0e',  # orange
    '#2ca02c',  # green
    '#d62728',  # red
    '#9467bd',  # purple
    '#8c564b',  # brown
    '#e377c2',  # pink
    '#7f7f7f',  # gray
    '#bcbd22',  # olive
    '#17becf'   # cyan
)
colors <- list(
    'Ly6C-hi Monocytes' = '#2ca02c',  # green
    'Neutrophils' = '#1f77b4',  # blue
    'MoMacs' = '#ff7f0e',  # orange
    'Macrophages' = '#d62728',  # red
    'DCs' = '#17becf',   # cyan
    'NK Cells' = '#bcbd22',  # olive
    'B cells' = '#8c564b',  # brown
    'CD4+ T cells' = '#e377c2',  # pink
    'CD8+ T cells' = '#9467bd',  # purple
    'AM' = '#bcbd22',  # olive
    'IM' = '#7f7f7f'  # gray
)
cell_types <- names(colors)
steady_state <- c("Untreated", "DMSO",  "PBS", "D0")

# Start Log
start_time = Sys.time()
log <- log_open(paste0("analyze_counts-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Preprocess Data

# Read data
raw_counts <- append_many_csv(
    file.path(wd, opt[['input-dir']]), recursive=FALSE, include_filepath=FALSE
)
for (col in c('viable_cells_conc', 'total_vol', 'total_viable_cells', 'abs_count')) {
    if (col %in% colnames(raw_counts)) {
        raw_counts[[col]] <- as.numeric(raw_counts[[col]])
    }
}
raw_counts[(raw_counts[, 'cell_type']=='cDCs'), "cell_type"] <- "DCs"
raw_counts[(raw_counts[, 'cell_type']=='Dendritic Cells'), "cell_type"] <- "DCs"
raw_counts <- raw_counts[raw_counts[['cell_type']] %in% cell_types, ]
raw_counts <- raw_counts[raw_counts[['treatment']] %in% steady_state, ]

# group by
df <- raw_counts %>%
    group_by(organ, cell_type, zygosity) %>%
    summarise(
        mean_pct_cells = mean(pct_cells, na.rm = TRUE),
        mean_abs_count = mean(abs_count, na.rm = TRUE),
        .groups = 'drop'
    )
df <- as.data.frame(df)


# ----------------------------------------------------------------------
# Plot

organs <- unique(df[['organ']])
pbar <- progress_bar$new(total = length(organs))
for (organ in sort(organs)) {

    fig <- plot_ribbon_chart(
        df = df[(df[['organ']]==organ), ],
        x = "zygosity", y = "mean_abs_count", subcat = "cell_type", title="Lung",
        bar_width = 0.5,
        normalize = FALSE,
        show_values=FALSE,
        colors = NULL,
        group_order = c("WT", "KO")
    )

    if (!troubleshooting) {
        save_fig(
            fig=fig,
            height=opt[['height']], width=opt[['width']],
            dirpath=file.path(wd, opt[['output-dir']], 'figures'),
            filename=paste(organ, 'ribbon', sep='-'),
            save_html=TRUE
        )
    }
    pbar$tick()
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()

