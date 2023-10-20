## Assist in designing flow panels

wd = dirname(this.path::here())  # wd = '~/github/R/harrisonRTools'
library('readxl')
library('tidyr')
library('optparse')
library('logr')
source(file.path(wd, 'R', 'functions', 'text_tools.R'))
source(file.path(wd, 'R', 'functions', 'list_tools.R'))  # find_first_match_index


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-p", "--panel-input"), default='data/flow/panel_input.xlsx',
                metavar='data/flow/panel_input.xlsx', type="character",
                help="specify the antibodies to use in your flow panel"),

    make_option(c("-a", "--antibody-inventory"), default='data/flow/antibody_inventory.xlsx',
                metavar='data/flow/antibody_inventory.xlsx', type="character",
                help="antibody inventory table"),

    make_option(c("-c", "--instrument-config"), default='data/flow/instrument_config.xlsx',
                metavar='data/flow/instrument_config.xlsx', type="character",
                help="instrument configuration file"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("design_flow_panel-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Format reference files

log_print(paste(Sys.time(), 'Reading data...'))


# instrument config
instr_cfg <- read_excel(file.path(wd, opt[['instrument-config']]))
colnames(instr_cfg) <- unlist(lapply(colnames(instr_cfg), title_to_snake_case))  # format columns
instr_cfg <- separate_rows(instr_cfg, 'fluorochrome_detected', sep=', ')  # split comma-separated list


# antibody inventory
ab_inv <- read_excel(file.path(wd, opt[['antibody-inventory']]))
ab_inv <- ab_inv[, 1:(find_first_match_index('\\.{3}\\d{2}', colnames(ab_inv))-1)]  # filter extra columns
colnames(ab_inv) <- unlist(lapply(colnames(ab_inv), title_to_snake_case))  # column names
colnames(ab_inv) <- unlist(lapply(colnames(ab_inv), function(x) gsub('[.]', '', x)))  # column nmaes

# remove 'c2 a0', the "no-break space"
# see: https://stackoverflow.com/questions/68982813/r-string-encoding-best-practice
for (col in colnames(ab_inv)) {
    ab_inv[(!is.na(ab_inv[[col]]) & ab_inv[[col]] == enc2utf8("\u00a0")), col] <- NA
}

# save
if (!troubleshooting) {
    directory = file.path(wd, dirname(opt[['antibody-inventory']]), 'troubleshooting')
    if (!dir.exists(directory)) {
        dir.create(directory, recursive=TRUE)
    }
    filepath = file.path(
        directory,
        paste0('_', tools::file_path_sans_ext(basename(opt[['antibody-inventory']])), '.csv')  # filename
    )
    write.csv(ab_inv, file = filepath, row.names = FALSE)
}




# ----------------------------------------------------------------------
# Do stuff


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
