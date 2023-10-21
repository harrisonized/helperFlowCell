## Assist in designing flow panels

wd = dirname(this.path::here())  # wd = '~/github/R/harrisonRTools'
library('readxl')
library('tidyr')
library('optparse')
library('logr')
source(file.path(wd, 'R', 'functions', 'text_tools.R'))  # text_strip
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
# Preprocessing


cleanup_antibody_names <- function(df, columns = c('antibody', 'alternative_name')) {

    for (col in columns) {

        # line endings
        df[[col]] = unlist( lapply(df[[col]], function(x) txt_strip(x, chars='() ') ) )

        # lower/capital
        df[[col]] = unlist( lapply(df[[col]], function(x) sub(".*^[Cc][Dd]", "CD", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) sub(".*^[Ii][Ll]", "IL", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) gsub('(IL)([0-9]+)', '\\1-\\2', x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) sub(".*^Ly-", "Ly", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) gsub("[Bb][Cc][Ll1]-{0,1}", "Bcl-", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) sub(".*^[Oo]nly", "Only", x)) )  # "Only clone name"
        df[[col]] = unlist( lapply(df[[col]], function(x) gsub(".*^[Tt][Cc][Rr]", "TCR", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) gsub(".*^[Tt][Dd][Tt]", "TdT", x)) )

        # animals
        df[[col]] = unlist( lapply(df[[col]], function(x) gsub("[Gg]oat", "goat", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) gsub("[Mm]ouse", "mouse", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) gsub("[Rr]at", "rat", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) gsub("[Dd]onkey", "donkey", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) gsub("[Bb]ovine", "bovine", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) gsub("[Ss]heep", "sheep", x)) )

        # special characters
        df[[col]] = unlist( lapply(df[[col]], function(x) sub("α", "a", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) sub("β", "b", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) sub("γ", "g", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) gsub("™", "", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) gsub(
            "([A-Za-z0-9],)\\s*([A-Za-z0-9])", "\\1 \\2", x)
        ))  # comma-separated list

        # specific genes
        df[[col]] = unlist( lapply(df[[col]], function(x) sub(" Fixable Viability Kit", "", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) sub("FoxP3", "Foxp3", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) sub("INFg", "IFNg", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) gsub("Immunoglobulin ", "Ig", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) sub("MHC II", "MHCII", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) sub("NK cell Pan", "CD49b", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) sub("NK cells", "CD49b", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) sub(".*^[Nn][Oo]tch", "Notch", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) sub("[Rr][Oo][Rr][gγy][yt]", "RORgt", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) sub("[Tt][Cc][Rr][Bb-]\\w*", "TCRb", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) sub("[Tt][Gg][Ff][Bb-]\\w*", "TGFb", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) sub("Vb8.1 Vb8.2", "Vb8.1, Vb8.2", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) sub("Vb8.1, 2", "Vb8.1, Vb8.2", x)) )
        df[[col]] = unlist( lapply(df[[col]], function(x) gsub(
            "\\(Tonegawa nomenclat", "(Tonegawa nomenclat)", x)
        ))  # removed above
    }

    return (df)
}





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

df <- cleanup_antibody_names(df)

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
