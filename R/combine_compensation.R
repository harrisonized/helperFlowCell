## Combine traditional and autospill
## 1. Filter negative values
## 2. Pick the minimum of the two compensation values
## 3. Output as XML


wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
library('optparse')
suppressPackageStartupMessages(library('logr'))
import::from(tidyr, 'pivot_longer')
import::from(XML, 'saveXML')
import::from(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'spillover_to_xml', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'reset_index', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input1"), default='data/compensation/Traditional.csv',
                metavar='data/compensation/Acquisition-defined.csv', type="character",
                help="specify the directory of json files"),

    make_option(c("-j", "--input2"), default="data/compensation/Autospill.csv",
                metavar="data/compensation", type="character",
                help="set the output directory for the data"),

    make_option(c("-o", "--output-dir"), default="data/compensation",
                metavar="data/compensation", type="character",
                help="set the output directory for the data"),

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

df1 <- read.csv(
    file.path(wd, opt[['input1']]),
    row.names=1, header=TRUE, check.names=FALSE
)
df2 <- read.csv(
    file.path(wd, opt[['input2']]),
    row.names=1, header=TRUE, check.names=FALSE
)

# Combine
df1[(df1 < 0)] <- NA
df2[(df2 < 0)] <- NA
df <- pmin(df1, df2, na.rm=TRUE) 
df[(is.na(df))] <- 0

# save
if (!troubleshooting) {
    write.table(
        df,
        file = file.path(wd, opt[['output-dir']], 'Combined.csv'),
        row.names=TRUE, col.names=NA, sep=','
    )
}


# ----------------------------------------------------------------------
# Export Spillover Table

channels <- row.names(df1)

# pivot
withCallingHandlers({
    spillover_table <- as.data.frame(pivot_longer(
        reset_index(df, index_name="- % Fluorochrome"),
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
# XML file
# Use this to copy directly into Flowjo workspace files

filename <- "combined"
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
