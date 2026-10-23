## Use this script to combine the results of many experiments
## 
## Be sure to  generate the tvds per experiment using
##     Rscript analyze_tvd.R -g treatment,gfp -f gfp
## Use "-g treatment,zygosity" instead to separate male and female
##

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
import::from(file.path(wd, 'R', 'config', 'user_input.R'),
    'excluded_cell_types', 'excluded_zygosities',
    'cell_type_to_color', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/tvd',
                metavar='data/tvd', type="character",
                help="input directory, all csv files will be read in"),

    make_option(c("-o", "--output-dir"), default="hfc-output",
                metavar="hfc-output", type="character",
                help="output directory"),

    make_option(c("-g", "--group-by"), default='treatment,gfp',
                metavar='treatment,gfp', type="character",
                help="enter a column or comma-separated list of columns, no spaces"),
    
    make_option(c("-w", "--make-whole"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="Use this for incomplete panels to add an 'Other' category to the Polar Area plot."),

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
metadata_cols <- unlist(strsplit(opt[['group-by']], ','))


# Start Log
start_time = Sys.time()
log <- log_open(paste0("combine_tvd-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Preprocess Data

# Read data
tvd_df <- append_many_csv(
    file.path(wd, opt[['input-dir']]), recursive=FALSE, include_filepath=FALSE,
    na_strings=c("NA")
)


# ----------------------------------------------------------------------
# Polar Area Charts

log_print(paste0(Sys.time(), ' Plotting polar area charts...'))

tmp_df <- tvd_df[
    (!tvd_df[['cell_type']] %in% excluded_cell_types) &
    (!tvd_df[['zygosity']] %in% excluded_zygosities), ]

polar_df <- tmp_df[(!is.na(tmp_df['cell_type'])), ] %>%
    group_by(across(all_of(c("organ", "cell_type", metadata_cols)))) %>%
    summarise(
        mean_pct_cells = mean(pct_cells, na.rm = TRUE),
        mean_tvd = mean(tvd, na.rm = TRUE),
        .groups = 'drop'
    )
polar_df <- as.data.frame(polar_df)
polar_df <- fillna(polar_df, 'mean_tvd', 0)

organs <- sort(unique(polar_df[['organ']]))
metadata_combinations <- unique(polar_df[, metadata_cols, drop = FALSE])

for (organ in organs) {
    for (i in seq_len(nrow(metadata_combinations))) {

        metadata <- metadata_combinations[i, , drop = FALSE]

        # Filter by organ and all metadata columns
        subset_df <- polar_df[polar_df[['organ']] == organ, ]

        for (col in metadata_cols) {
            subset_df <- subset_df[
                !is.na(subset_df[[col]]) &
                subset_df[[col]] == metadata[[col]][1],
            ]
        }
        
        fig <- plot_polar_area(
            subset_df,
            x="mean_pct_cells", y="mean_tvd", group='cell_type',
            min_label_pct=100,
            colors=cell_type_to_color,
            title=organ,
            show_axes=TRUE,
            fill_remainder=opt[['make-whole']]
        )

        if (!troubleshooting) {
            save_fig(
                fig=fig,
                height=500, width=800,
                dirpath=file.path(path.expand(wd), opt[['output-dir']], 'figures', 'polar_area'),
                filename=paste(organ, 'polar_area', paste(metadata, collapse = "_"), sep='-'),
                save_html=TRUE
            )
        }
    }
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
