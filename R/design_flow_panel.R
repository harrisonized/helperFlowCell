## Assist in designing flow panels

wd = dirname(this.path::here())  # wd = '~/github/R/harrisonRTools'
library('readxl')
library('tidyr')
library('optparse')
library('logr')
source(file.path(wd, 'R', 'functions', 'text_tools.R'))  # title_to_snake_case
source(file.path(wd, 'R', 'functions', 'list_tools.R'))  # find_first_match_index, multiple_replacement


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-file"), default='data/flow/input_file.xlsx',
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
# Configs

# cleanup antibody names
replacements <- c(

    # lower/capital
    ".*^[Cc][Dd]" = "CD",
    ".*^[Ii][Ll]" = "IL",
    "(IL)([0-9]+)" = "\\1-\\2",
    ".*^Ly-" = "Ly",
    "[Bb][Cc][Ll1]-{0,1}" = "Bcl-",
    ".*^[Oo]nly" = "Only",
    ".*^[Tt][Cc][Rr]" = "TCR",
    ".*^[Tt][Dd][Tt]" = "TdT",

    # animals
    "[Gg][Oo][Aa][Tt]" = "goat",
    "[Mm][Oo][Uu][Ss][Ee]" = "mouse",
    "[Rr][Aa][Tt]" = "rat",
    "[Dd][Oo][Nn][Kk][Ee][Yy]" = "donkey",
    "[Bb][Oo][Vv][Ii][Nn][Ee]" = "bovine",
    "[Ss][Hh][Ee][Ee][Pp]" = "sheep",

    # special characters
    "α" = "a",
    "β" = "b",
    "γ" = "g",
    "™" = "",
    "([A-Za-z0-9],)\\s*([A-Za-z0-9])" = "\\1 \\2",  # comma-separated list

    # specific genes
    " Fixable Viability Kit" = "",
    "FoxP3" = "Foxp3",
    "INFg" = "IFNg",
    "Immunoglobulin " = "Ig",
    "MHC II" = "MHCII",
    "NK cell Pan" = "CD49b",
    ".*^[Nn][Oo]tch" = "Notch",
    "[Rr][Oo][Rr][gγy][yt]" = "RORgt",
    "[Tt][Cc][Rr][Bb-]\\w*" = "TCRb",
    "[Tt][Gg][Ff][Bb-]\\w*" = "TGFb",
    "Vb8.1 Vb8.2" = "Vb8.1, Vb8.2",
    "Vb8.1, 2" = "Vb8.1, Vb8.2"
)


# ----------------------------------------------------------------------
# Format reference files

log_print(paste(Sys.time(), 'Reading data...'))


# instrument config
instr_cfg <- read_excel(file.path(wd, opt[['instrument-config']]))
colnames(instr_cfg) <- unlist(lapply(colnames(instr_cfg), title_to_snake_case))  # format columns
instr_cfg <- separate_rows(instr_cfg, 'fluorochrome_detected', sep=', ')  # split comma-separated list


# antibody inventory
df <- read_excel(file.path(wd, opt[['antibody-inventory']]))
df <- df[, 1:(find_first_match_index('\\.{3}\\d{2}', colnames(df))-1)]  # filter extra columns
colnames(df) <- unlist(lapply(colnames(df), title_to_snake_case))  # column names
colnames(df) <- unlist(lapply(colnames(df), function(x) gsub('[.]', '', x)))  # column nmaes

# remove 'c2 a0', the "no-break space"
# see: https://stackoverflow.com/questions/68982813/r-string-encoding-best-practice
for (col in colnames(df)) {
    df[(!is.na(df[[col]]) & df[[col]] == enc2utf8("\u00a0")), col] <- NA
}

# cleanup antibody names
for (col in c('antibody', 'alternative_name')) {
    df[[col]] = multiple_replacement(df[[col]], replacements, func='gsub')
}

# TODO: cleanup_fluorophore_names

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
    write.table(df, file = filepath, row.names = FALSE, sep=',')
}

# save
if (!troubleshooting) {
    directory = file.path(wd, dirname(opt[['antibody-inventory']]), 'troubleshooting')
    if (!dir.exists(directory)) {
        dir.create(directory, recursive=TRUE)
    }

    filepath = file.path(directory, '_antibodies.txt')
    write.table(sort(unique(df[['antibody']])), filepath,
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    filepath = file.path(directory, '_alternative_names.txt')
    write.table(sort(unique(df[['alternative_name']])), filepath,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
}


# ----------------------------------------------------------------------
# Do stuff


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
