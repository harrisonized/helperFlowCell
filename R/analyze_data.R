## Compute cell proportions from each organ using exported flowjo data

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
suppressPackageStartupMessages(library('dplyr'))
library('optparse')
suppressPackageStartupMessages(library('logr'))
import::from(magrittr, '%>%')
import::from(stringr, 'str_detect', 'str_extract')
import::from(tidyr, 'pivot_longer')
import::from(plyr, 'rbind.fill')
import::from(dplyr, 'group_by', 'summarize')
import::from(ggplot2, 'ggsave')
import::from(plotly, 'save_image')
import::from(htmlwidgets, 'saveWidget')  # brew install pandoc

import::from(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'preprocess_flowjo_export', .character_only=TRUE)

import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'rename_columns', 'reset_index', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'dict_zip', 'interleave', 'items_in_a_not_b', 'multiple_replacement',
    .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_dots', 'plot_violin', .character_only=TRUE)
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

    make_option(c("-o", "--output-dir"), default="data/analysis",
                metavar="data/analysis", type="character",
                help="output directory for the data"),

    make_option(c("-f", "--figures-dir"), default="figures/analysis",
                metavar="figures/analysis", type="character",
                help="output directory for the figures"),

    make_option(c("-r", "--ref-dir"), default='ref/mice',
                metavar='ref/mice', type="character",
                help="directory of files containing all the mouse data"),

    make_option(c("-s", "--save-html"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="save html files in addition to PNG"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]
output_dir <- file.path(wd, opt[['output-dir']])
troubleshooting_dir = file.path(output_dir, 'troubleshooting')


# Start Log
start_time = Sys.time()
log <- log_open(paste0("analyze_data-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Raw Data

log_print(paste(Sys.time(), 'Reading data...'))

# reference files
mouse_db <- append_many_csv(file.path(wd, opt[['ref-dir']]), recursive=TRUE, include_filepath=FALSE)
flow_metadata <- append_many_csv(file.path(wd, opt[['metadata-dir']]), recursive=TRUE, include_filepath=FALSE)
# drop duplicates for mouse_db and flow_metadata

# Read counts data exported from flowjo
counts <- append_many_csv(file.path(wd, opt[['input-dir']]), recursive=TRUE)
counts <- preprocess_flowjo_export(counts)
counts <- counts[!str_detect(counts[['fcs_name']], '(U|u)nstained'), ]  # drop unstained cells


# ----------------------------------------------------------------------
# Preprocessing

# Reshape
df <- pivot_longer(counts,
    names_to = "gate", values_to = "num_cells",
    cols=items_in_a_not_b(colnames(counts), c(id_cols, initial_gates)),
    values_drop_na = TRUE
)
df[['cell_type']] <- unlist(lapply(strsplit(df[['gate']], '/'), function(x) x[length(x)]))
df[['cell_type']] <- multiple_replacement(df[['cell_type']], cell_type_spell_check)
df[['pct_cells']] <- round(df[['num_cells']] /
    df[['Cells/Single Cells/Single Cells/Live Cells']] * 100, 2)

# join reference files
df <- merge(df, flow_metadata,
    by="fcs_name", all.x=FALSE, all.y=FALSE)  # left join
df <- merge(df, mouse_db[, items_in_a_not_b(colnames(mouse_db), mouse_db_ignore)],
    by="mouse_id", all.x=TRUE, all.y=FALSE)  # left join
df <- df[(df[['is_unstained']]==FALSE), ]  # drop unstained cells
df[['weeks_old']] <- round(df[['age']]/7, 1)


# unpivot mNeonGreen+ num_cells into its own column
if ( any(str_detect(colnames(counts), 'mNeonGreen')) ) {

    # split out mNeonGreen+ rows
    fp_quant <- df[(df[['cell_type']]=='mNeonGreen+'), c(id_cols, 'gate', 'num_cells')]

    # second-to-last gate
    fp_quant[['gate']] <- unlist(lapply(
        strsplit(fp_quant[['gate']], '/'), function(x) paste0(x[1:length(x)-1], collapse='/')
    ))
    fp_quant[['cell_type']] <- unlist(lapply(
        strsplit(fp_quant[['gate']], '/'), function(x) x[length(x)]
    ))
    fp_quant <- rename_columns(fp_quant, c('num_cells'='num_mneongreen_pos'))
    
    df = merge(
        df[(df[['cell_type']]!='mNeonGreen+'), ],
        fp_quant[, c(id_cols, 'gate', 'num_mneongreen_pos')],
        by=c(id_cols, 'gate'), all.x=TRUE, all.y=FALSE
    )
    df[['pct_mneongreen_pos']] <- round(df[['num_mneongreen_pos']] / df[['num_cells']] * 100, 2)
}


# filters
for (cell_type in cell_type_ignore) {
    df <- df[!str_detect(df[['cell_type']], cell_type), ]
}


# save
if (!troubleshooting) {
    if (!dir.exists(file.path(wd, opt[['output-dir']]))) {
        dir.create(file.path(wd, opt[['output-dir']]), recursive=TRUE)
    }
    filepath = file.path(wd, opt[['output-dir']], 'cell_proportions.csv')  # filename
    write.table(df, file = filepath, row.names = FALSE, sep = ',' )
}


# ----------------------------------------------------------------------
# Plot cell type frequency per organ

log_print(paste(Sys.time(), 'Quantifying cell type frequencies...'))

organs <- sort(unique(df[['organ']]))
log_print(paste(Sys.time(), 'Groups found...', paste(organs, collapse = ', ')))
for (organ in sort(organs)) {

    log_print(paste(Sys.time(), 'Processing...', organ))

    fig <- plot_violin(
        df[(df[['organ']]==organ) & (df[['num_cells']]>10), ],
        x='cell_type', y='pct_cells', group_by='zygosity',
        ylabel='Percent of Live Cells', title=organ,
        ymin=0, ymax=100,
        hover_data=c(
            'mouse_id', 'sex', 'zygosity', 'treatment', 'weeks_old', 
            'Cells', 'num_cells', 'fcs_name'),
        color_discrete_map=c(
            'heterozygous'='#2ca02c',  # green
            'wild type'='#62c1e5'  # blue
        )
    )

    # export
    if (!troubleshooting) {

        if (!dir.exists( file.path(wd, opt[['figures-dir']], 'cell_proportions') )) {
            dir.create( file.path(wd, opt[['figures-dir']], 'cell_proportions'), recursive=TRUE)
        }

        # save PNG
        suppressWarnings(save_image(fig,
            file=file.path(wd, opt[['figures-dir']], 'cell_proportions',
                paste0('violin-pct_cells-', organ, '.png')),  # filename
            height=500, width=800, scale=3
        ))

        # save HTML
        if (opt[['save-html']]) {
            if (!dir.exists( file.path(wd, opt[['figures-dir']], 'cell_proportions', 'html') )) {
                dir.create( file.path(wd, opt[['figures-dir']], 'cell_proportions', 'html'), recursive=TRUE)
            }
            suppressWarnings(saveWidget(
                widget = fig,
                file=file.path(wd, opt[['figures-dir']], 'cell_proportions', 'html',
                    paste0('violin-pct_cells-', organ, '.html')),  # filename
                selfcontained = TRUE
            ))
            unlink(file.path(
                wd, opt[['figures-dir']], 'cell_proportions', 'html',
                paste0('violin-pct_cells-', organ, '_files')
            ), recursive=TRUE)            
        }
    }
}


# ----------------------------------------------------------------------
# Plot mNeonGreen positivity per cell type

if ('pct_mneongreen_pos' %in% colnames(df)) {
    log_print(paste(Sys.time(), 'Quantifying mNeonGreen...'))

    for (organ in sort(organs)) {
        if ( !all(is.na(df[(df[['organ']]==organ), 'pct_mneongreen_pos'])) ) {

            log_print(paste(Sys.time(), 'Processing...', organ))

            fig <- plot_violin(
                df[(df[['organ']]==organ) & (df[['num_cells']]>10), ],
                x='cell_type', y='pct_mneongreen_pos', group_by='zygosity',
                ylabel='Percent mNeonGreen+', title=organ,
                ymin=0, ymax=100,
                hover_data=c('mouse_id', 'sex', 'zygosity', 'treatment', 'weeks_old',
                             'Cells', 'num_cells', 'num_mneongreen_pos', 'fcs_name'),
                color_discrete_map=c(
                    'heterozygous'='#2ca02c',  # green
                    'wild type'='#62c1e5'  # blue
                )
            )

            # export
            if (!troubleshooting) {

                if (!dir.exists( file.path(wd, opt[['figures-dir']], 'mneongreen') )) {
                    dir.create( file.path(wd, opt[['figures-dir']], 'mneongreen'), recursive=TRUE)
                }

                # save PNG
                suppressWarnings(save_image(fig,
                    file=file.path(wd, opt[['figures-dir']], 'mneongreen',
                        paste0('violin-pct_mneongreen_pos-', organ, '.png')),  # filename
                    height=500, width=800, scale=3
                ))

                # save HTML
                if (opt[['save-html']]) {
                    if (!dir.exists( file.path(wd, opt[['figures-dir']], 'mneongreen', 'html') )) {
                        dir.create( file.path(wd, opt[['figures-dir']], 'mneongreen', 'html'), recursive=TRUE)
                    }
                    suppressWarnings(saveWidget(
                        widget = fig,
                        file=file.path(wd, opt[['figures-dir']], 'mneongreen', 'html',
                            paste0('violin-pct_mneongreen_pos-', organ, '.html')),  # filename
                        selfcontained = TRUE
                    ))
                    unlink(file.path(
                        wd, opt[['figures-dir']], 'mneongreen', 'html',
                        paste0('violin-pct_mneongreen_pos-', organ, '_files')
                    ), recursive=TRUE)
                }
            }
        }
    }
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
