## Extract spectra from json data downloaded from Biolegend's Spectra Analyzer
## See: https://www.biolegend.com/spectraanalyzer

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
library('rjson')
library('optparse')
library('logr')
import::from(file.path(wd, 'R', 'utils', 'file_io.R'),
    'join_many_csv', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='ref/raw_spectra_json',
                metavar='ref/raw_spectra_json', type="character",
                help="specify the directory of json files"),

    make_option(c("-o", "--output-dir"), default="ref/raw_spectra",
                metavar="ref/raw_spectra", type="character",
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


# ----------------------------------------------------------------------
# Convert each file

# get files
filepaths <- list.files(file.path(wd, opt[['input-dir']]), full.names=TRUE)
log_print(paste('Number of files found:', length(filepaths)))

# convert
for (filepath in filepaths) {

    log_print( paste('Processing:', basename(filepath)) )
    filename_sans_ext <- tools::file_path_sans_ext(basename(filepath))
    
    # convert
    raw_data <- fromJSON(file=filepath)
    df <- as.data.frame(do.call(rbind, raw_data[['data']]))
    colnames(df) <- c('Wavelength', 'Intensity')

    # normalize
    max_intensity <- max(df[2])
    if (max_intensity > 100) {
        df[2] <- df[2] / max_intensity
    } else if (max_intensity > 1) {
        df[2] <- df[2] / 100
    }

    # save
    if (!troubleshooting) {
        write.table(
            df,
            file = file.path(opt[['output-dir']], paste0(filename_sans_ext, '.csv')),
            row.names = FALSE, sep=','
        )
    }
}


# ----------------------------------------------------------------------
# Combine output

log_print(paste('Combining...'))

filepaths <- list.files(file.path(wd, opt[['output-dir']]), full.names=TRUE)
combined <- join_many_csv(
    filepaths,
    index_cols='Wavelength', value_cols='Intensity',
    all.x=TRUE, all.y=TRUE
)

# remove half wavelengths
combined <- combined[(combined[['Wavelength']]%%1 < 0.01), ]

# save
if (!troubleshooting) {
    write.table(
        combined,
        file = file.path(dirname(opt[['output-dir']]), paste0('combined-spectra', '.csv')),
        row.names = FALSE, sep=',',
        na=""
    )
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
