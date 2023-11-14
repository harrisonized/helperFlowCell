## Standardizes the names in instrument_config and antibody_inventory files

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
library('tidyr')
library('optparse')
library('logr')
source(file.path(wd, 'R', 'preprocessing.R'))
source(file.path(wd, 'R', 'functions', 'file_io.R'))  # read_excel_or_csv


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-c", "--instrument-config"), default='ref/instrument_config.xlsx',
                metavar='ref/instrument_config.xlsx', type="character",
                help="instrument configuration file"),

    make_option(c("-a", "--antibody-inventory"), default='ref/antibody_inventory.xlsx',
                metavar='ref/antibody_inventory.xlsx', type="character",
                help="antibody inventory table"),

    make_option(c("-o", "--output-dir"), default="data/troubleshooting",
                metavar="data/troubleshooting", type="character",
                help="set the output directory for the data")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
output_dir <- file.path(wd, opt[['output-dir']])

if (!dir.exists(file.path(output_dir))) {
    dir.create(file.path(output_dir), recursive=TRUE)
}


# Start Log
start_time = Sys.time()
log <- log_open(paste0("fix_names-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Instrument Config

instr_cfg <- read_excel_or_csv(file.path(wd, opt[['instrument-config']]))
instr_cfg <- preprocess_instrument_config(instr_cfg)
instr_cfg_long <- separate_rows(instr_cfg, 'fluorophore', sep=', ')


# full instrument config
filepath = file.path(output_dir, 
    paste0('_', tools::file_path_sans_ext(basename(opt[['instrument-config']])), '.csv')
)
write.table(instr_cfg, file = filepath, row.names = FALSE, sep=',')


# ----------------------------------------------------------------------
# Antibody inventory

ab_inv <- read_excel_or_csv(file.path(wd, opt[['antibody-inventory']]))
ab_inv <- preprocess_antibody_inventory(ab_inv)


# full antibody inventory
filepath = file.path(output_dir, 
    paste0('_', tools::file_path_sans_ext(basename(opt[['antibody-inventory']])), '.csv')
)
write.table(ab_inv, file = filepath, na = "", row.names = FALSE, sep=',')


# fluorophores in ab_inv not in inst_cfg
unavailable_fluorophores = sort(items_in_a_not_b(
    unique( ab_inv[['fluorophore']] ), unique( instr_cfg_long[['fluorophore']] )
))
filepath = file.path(output_dir, 'unavailable_fluorophores.txt')
write.table(unavailable_fluorophores, filepath,
            row.names = FALSE, col.names = FALSE, quote = FALSE)


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
