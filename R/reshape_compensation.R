## Reshape compensation matrix for manual input into FACSDiva

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
library('optparse')
suppressPackageStartupMessages(library('logr'))
import::from(tidyr, 'pivot_longer', 'pivot_wider')
import::from(XML, 'saveXML')
import::from(matlib, 'inv')

import::from(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'spillover_to_xml', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'join_many_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'rename_columns', 'reset_index', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', 'list2matrix', 'move_list_items_to_front', .character_only=TRUE)


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
# Preprocessing

spillover_matrix <- read.csv(
    file.path(wd, opt[['input']]),
    row.names=1, header=TRUE, check.names=FALSE
)
channels <- row.names(spillover_matrix)

# sort
if (file.exists(file.path(wd, opt[['sort']]))) {
    matrix_order <- read.csv(file.path(wd, opt[['sort']]), header=FALSE)[['V1']]
    matrix_order <- move_list_items_to_front(matrix_order, colnames(spillover_matrix))
    spillover_matrix <- spillover_matrix[matrix_order, matrix_order]
}

# pivot
withCallingHandlers({
    spillover_table <- as.data.frame(pivot_longer(
        reset_index(spillover_matrix, index_name="- % Fluorochrome"),
        cols=channels,
        names_to = "Fluorochrome",
        values_to = "Spectral Overlap"
    ))
}, warning = function(w) {
    if ( any(grepl("Using an external vector in selections", w)) ) {
        invokeRestart("muffleWarning")
    }
})


# ----------------------------------------------------------------------
# Export Spillover Table

# save
if (!troubleshooting) {

    # cleanup
    diva_spillover <- spillover_table[
        (spillover_table[["Fluorochrome"]] != spillover_table[["- % Fluorochrome"]]),
        c("Fluorochrome", "- % Fluorochrome", "Spectral Overlap")
    ]
    diva_spillover[["Spectral Overlap"]] <- diva_spillover[["Spectral Overlap"]]*100
    diva_spillover[["Spectral Overlap"]] <- round(diva_spillover[["Spectral Overlap"]], 2)

    write.table(
        diva_spillover,
        file = file.path(wd, opt[['output-dir']], 'spillover_table.csv'),
        row.names = FALSE, sep=',',
        na=""
    )
}

# ----------------------------------------------------------------------
# Compensation Matrix
# "The compensation matrix is the inverse of the spillover matrix." - Data File Standard for FCS3.1
# Use this to import compensation matrix directly into FACSDiva
# Experiment > [Import Cytometer Settings] 


compensation_matrix <- inv(as.matrix(spillover_matrix))

# save
if (!troubleshooting) {
    write.table(
        t(as.vector(compensation_matrix)),  # flattened
        file = file.path(wd, opt[['output-dir']], 'compensation_matrix.txt'),
        row.names=FALSE, col.names=FALSE, sep=','
    )
}


# ----------------------------------------------------------------------
# XML file
# Use this to copy directly into Flowjo workspace files

filename <- basename(tools::file_path_sans_ext(opt[['input']]))
spillover_xml <- spillover_to_xml(spillover_table, channels, name=filename)

# save
if (!troubleshooting) {
    invisible(saveXML(spillover_xml,
        prefix='<?xml version="1.0" encoding="UTF-8"?>\n',
        file = file.path(wd, opt[['output-dir']], paste0(filename, '.mtx'))
    ))  # cat() to view
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
