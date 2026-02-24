## Creates a polar area plot where angle represents percentage cells and
## radius represents the TVD value
## Run concat_data before running this

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('logr'))
library('optparse')
import::from(progress, 'progress_bar')
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'save_fig', 'plot_polar_area', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'group_by_agg', 'fillna', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings


# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='database/tvd/steady-state',
                metavar='database/tvd/steady-state', type="character",
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
steady_state <- c("untreated", "DMSO",  "PBS")

# Start Log
start_time = Sys.time()
log <- log_open(paste0("analyze_counts-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Preprocess Data

# Read data
raw_tvds <- append_many_csv(
    file.path(wd, opt[['input-dir']]), recursive=FALSE, include_filepath=FALSE
)
raw_tvds[['tvd']] <- as.numeric(raw_tvds[['tvd']])
raw_tvds[(raw_tvds[, 'organ']=='sp'), "organ"] <- "spleen"
raw_tvds[(raw_tvds[, 'organ']=='sp_dc'), "organ"] <- "spleen"
raw_tvds[(raw_tvds[, 'cell_type']=='Dendritic Cells'), "cell_type"] <- "DCs"

raw_tvds <- raw_tvds[raw_tvds[['cell_type']] %in% cell_types, ]
raw_tvds <- raw_tvds[raw_tvds[['treatment']] %in% steady_state, ]

# group by
df <- raw_tvds %>%
    group_by(organ, cell_type) %>%
    summarise(
        mean_pct_cells = mean(pct_cells, na.rm = TRUE),
        mean_tvd = mean(tvd, na.rm = TRUE),
        .groups = 'drop'
    )
df <- as.data.frame(df)
df <- fillna(df, 'mean_tvd', 0)


# ----------------------------------------------------------------------
# Plot


organs <- unique(df[['organ']])
pbar <- progress_bar$new(total = length(organs))
for (organ in sort(organs)) {

    fig <- plot_polar_area(df[(df[['organ']]==organ), ],
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
            dirpath=file.path(wd, opt[['output-dir']], 'figures'),
            filename=paste(organ, 'polar_area', sep='-'),
            save_html=TRUE
        )
    }
    pbar$tick()
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()

