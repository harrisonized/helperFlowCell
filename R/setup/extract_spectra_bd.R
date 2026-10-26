## Extract spectra from svg data downloaded from BD's Spectrum Viewer
## See: https://www.bdbiosciences.com/en-us/resources/bd-spectrum-viewer

wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/helperFlowCell'
library('xml2')
library('optparse')
library('logr')
import::from(stringr, 'str_extract_all', 'str_match_all')


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/raw_spectra',
                metavar='data/raw_spectra', type="character",
                help="specify the directory of svg files"),

    make_option(c("-o", "--output-dir"), default="data/raw_spectra/output",
                metavar="data/raw_spectra/output", type="character",
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
log <- log_open(paste0("extract_spectra-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))

# create output dir
if (!dir.exists(file.path(opt[['output-dir']]))) {
    dir.create(file.path(opt[['output-dir']]), recursive=TRUE)
}


#' Extract 
#' Function to extract x,y coordinates from an SVG path
#' 
extract_xy_coordinates <- function(doc, path_id) {
    # Find path by id
    path <- xml_find_first( doc, sprintf(".//*[@id='%s']", path_id) )
    if (inherits(path, "xml_missing")) {
        stop(paste("Path not found:", path_id))
    }

    # Get the path data
    path_data <- xml_attr(path, "d")
    matches <- str_match_all(path_data,
        "[ML]\\s*(-?[0-9]+\\.?[0-9]*),\\s*(-?[0-9]+\\.?[0-9]*)"
    )[[1]]  # Extract coordinate pairs following M or L commands

    # Create data frame
    df <- data.frame(
        x = as.numeric(matches[, 2]),
        y = as.numeric(matches[, 3])
    )

    return(df)
}


#' Convert the raw data to actual coordinates 
#' 
convert_coordinates <- function(data) {
    data$wavelength_nm <- 300 + (data$x - 60) / 1107 * 600
    data$intensity_percent <- (148 - data$y) / 118
    return(data)
}


# ----------------------------------------------------------------------
# Convert each file

# get files
filepaths <- list.files(file.path(wd, opt[['input-dir']]), full.names=TRUE, recursive=FALSE)
filepaths <- filepaths[!file.info(filepaths)[['isdir']]]
log_print(paste('Number of files found:', length(filepaths)))

# convert
for (filepath in filepaths) {

    log_print( paste('Processing:', basename(filepath)) )
    doc <- read_xml(filepath)

    path_ids <- str_extract_all(as.character(doc), "e[mx]Path_\\d{3}")[[1]]
    for (path_id in path_ids) {
        df <- extract_xy_coordinates(doc, path_id)
        df <- convert_coordinates(df)

        # Save
        if (!troubleshooting) {
            dirpath <- file.path(opt[['output-dir']],
                tools::file_path_sans_ext(basename(filepath)))
            if (!dir.exists(dirpath)) {
                dir.create(dirpath, recursive=TRUE)
            }
            write.table(df,
                file=file.path(dirpath, paste0(path_id, '.csv')),
                row.names=FALSE, sep=','
            )
        }
    }
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
