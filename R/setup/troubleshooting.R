## troubleshootings files that aid in troubleshooting the data

wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/helperFlowCell'
library('tidyr')
library('optparse')
library('logr')

import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'multiple_replacement', 'items_in_a_not_b', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'read_excel_or_csv', .character_only=TRUE)

import::from(file.path(wd, 'R', 'config', 'replacements.R'),
    'fluorophore_replacements', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'preprocess_instrument_config', 'preprocess_antibody_inventory',
    .character_only=TRUE)


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

    make_option(c("-s", "--spectra-file"), default='ref/spectra.csv',
                metavar='ref/spectra.csv', type="character",
                help="spectra data, main source is FPBase"),

    make_option(c("-o", "--output-dir"), default="ref/troubleshooting",
                metavar="ref/troubleshooting", type="character",
                help="set the output directory for the data")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting_dir <- file.path(wd, opt[['output-dir']])

if (!dir.exists(file.path(troubleshooting_dir))) {
    dir.create(file.path(troubleshooting_dir), recursive=TRUE)
}


# Start Log
start_time = Sys.time()
log <- log_open(paste0("troubleshooting-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Instrument Config

instr_cfg <- read_excel_or_csv(file.path(wd, opt[['instrument-config']]))
instr_cfg <- preprocess_instrument_config(instr_cfg)
instr_cfg_long <- separate_rows(instr_cfg, 'fluorophore', sep=', ')


# full instrument config
filepath = file.path(troubleshooting_dir, 
    paste0('_', tools::file_path_sans_ext(basename(opt[['instrument-config']])), '.csv')
)
write.table(instr_cfg, file = filepath, row.names = FALSE, sep=',')


# ----------------------------------------------------------------------
# Antibody Inventory

ab_inv <- read_excel_or_csv(file.path(wd, opt[['antibody-inventory']]))
ab_inv <- preprocess_antibody_inventory(ab_inv)

# full antibody inventory
filepath = file.path(troubleshooting_dir, 
    paste0('_', tools::file_path_sans_ext(basename(opt[['antibody-inventory']])), '.csv')
)
write.table(ab_inv, file = filepath, na = "", row.names = FALSE, sep=',')


# ----------------------------------------------------------------------
# Spectra file

all_spectra <- read_excel_or_csv(file.path(wd, opt[['spectra-file']]))

# preprocess column names
colnames(all_spectra) <- unname(multiple_replacement(
    colnames(all_spectra), fluorophore_replacements
))
cols <- colnames(all_spectra)
available_fluorophores <- sort(unique(
    gsub( '( EM| EX| AB)', '', cols[2:length(cols)] )
))


# ----------------------------------------------------------------------
# Write files


# all fluorophores
all_fluorophores = sort(unique( ab_inv[['fluorophore']] ))
filepath = file.path(troubleshooting_dir, 'all_fluorophores.txt')
write.table(all_fluorophores, filepath,
            row.names = FALSE, col.names = FALSE, quote = FALSE)


# fluorophores in ab_inv not in inst_cfg
unavailable_fluorophores = sort(items_in_a_not_b(
    unique( ab_inv[['fluorophore']] ), unique( instr_cfg_long[['fluorophore']] )
))
filepath = file.path(troubleshooting_dir, 'unavailable_fluorophores.txt')
write.table(unavailable_fluorophores, filepath,
            row.names = FALSE, col.names = FALSE, quote = FALSE)


# fluorophores in ab_inv not in spectra
unavailable_fluorophores = sort(items_in_a_not_b(
    unique( ab_inv[['fluorophore']] ), available_fluorophores
))
filepath = file.path(troubleshooting_dir, 'unavailable_spectra.txt')
write.table(unavailable_fluorophores, filepath,
            row.names = FALSE, col.names = FALSE, quote = FALSE)


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
