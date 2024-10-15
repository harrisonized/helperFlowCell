## Reshape compensation matrix for manual input into FACSDiva

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
library('optparse')
suppressPackageStartupMessages(library('logr'))
import::from(tidyr, 'pivot_longer')
import::from(XML, 'saveXML')

import::from(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'compensation_csv_to_mtx', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'join_many_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'rename_columns', 'reset_index', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', 'move_list_items_to_front', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input"), default='data/compensation/Acquisition-defined.csv',
                metavar='data/compensation/Acquisition-defined.csv', type="character",
                help="specify the directory of json files"),

    make_option(c("-o", "--output-dir"), default="data/compensation",
                metavar="data/compensation", type="character",
                help="set the output directory for the data"),

    make_option(c("-s", "--sort"), default='ref/fluorochrome_order.txt',
                metavar='ref/fluorochrome_order.csv', type="character",
                help="csv file containing a list"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("reshape_compensation-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))

# create output dir
if (!dir.exists(file.path(opt[['output-dir']]))) {
    dir.create(file.path(opt[['output-dir']]), recursive=TRUE)
}


# ----------------------------------------------------------------------
# Main

# read data
compensation_matrix <- read.csv(
    file.path(wd, opt[['input']]),
    row.names=1, header=TRUE, check.names=FALSE
)


# ----------------------------------------------------------------------
# Pivot

# sort
if (file.exists(file.path(wd, opt[['sort']]))) {
    matrix_order <- read.csv(file.path(wd, opt[['sort']]), header=FALSE)[['V1']]
    matrix_order <- move_list_items_to_front(matrix_order, colnames(compensation_matrix))
    compensation_matrix <- compensation_matrix[matrix_order, matrix_order]
}
compensation_matrix <- reset_index(compensation_matrix, index_name="- % Fluorochrome")
channels <- items_in_a_not_b(colnames(compensation_matrix), "- % Fluorochrome")

# pivot
withCallingHandlers({
    cytometer_settings <- as.data.frame(pivot_longer(
        compensation_matrix,
        cols=channels,
        names_to = "Fluorochrome",
        values_to = "Spectral Overlap"
    ))
}, warning = function(w) {
    if ( any(grepl("Using an external vector in selections", w)) ) {
        invokeRestart("muffleWarning")
    }
})


# xml
filename <- basename(tools::file_path_sans_ext(opt[['input']]))
xml <- compensation_csv_to_mtx(cytometer_settings, channels, name=filename)
if (!troubleshooting) {

    invisible(saveXML(xml,
        prefix='<?xml version="1.0" encoding="UTF-8"?>\n',
        file = file.path(wd, opt[['output-dir']], paste0(filename, '.mtx'))
    ))  # cat() to view
}

# save
if (!troubleshooting) {

    # cleanup
    cytometer_settings <- cytometer_settings[
        (cytometer_settings[["Fluorochrome"]] != cytometer_settings[["- % Fluorochrome"]]),
        c("Fluorochrome", "- % Fluorochrome", "Spectral Overlap")
    ]
    cytometer_settings[["Spectral Overlap"]] <- cytometer_settings[["Spectral Overlap"]]*100
    cytometer_settings[["Spectral Overlap"]] <- round(cytometer_settings[["Spectral Overlap"]], 2)

    write.table(
        cytometer_settings,
        file = file.path(wd, opt[['output-dir']], 'cytometer_settings.csv'),
        row.names = FALSE, sep=',',
        na=""
    )
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
