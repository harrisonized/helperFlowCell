## Plot the spectra for fluorophores used

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
library('tidyr')
library('ggplot2')
library('optparse')
library('logr')
source(file.path(wd, 'R', 'replacements.R'))
source(file.path(wd, 'R', 'preprocessing.R'))
source(file.path(wd, 'R', 'functions', 'file_io.R'))  # read_excel_or_csv, read_csv_from_text
source(file.path(wd, 'R', 'functions', 'list_tools.R'))  # multiple_replacement


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-file"), default='data/panel.csv',
                metavar='data/panel.csv', type="character",
                help="specify the fluorophores to use in your flow panel"),

    make_option(c("-a", "--antibody-inventory"), default='ref/antibody_inventory.xlsx',
                metavar='ref/antibody_inventory.xlsx', type="character",
                help="antibody inventory table"),

    make_option(c("-s", "--spectra-file"), default='ref/FPBase_Spectra.csv',
                metavar='ref/FPBase_Spectra.csv', type="character",
                help="spectra data downloaded from FPBase"),

    make_option(c("-o", "--output-dir"), default="data",
                metavar="data", type="character",
                help="set the output directory for the data"),

    make_option(c("-f", "--figures-dir"), default="figures",
                metavar="figures", type="character",
                help="set the output directory for the figures"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]
output_dir <- file.path(wd, opt[['output-dir']])
troubleshooting_dir = file.path(output_dir, 'troubleshooting')

# troubleshooting = TRUE

# Start Log
start_time = Sys.time()
log <- log_open(paste0("plot_spectra-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Spectra file

# Read data
all_spectra <- read_excel_or_csv(file.path(wd, opt[['spectra-file']]))

# preprocess column names
cols <- colnames(all_spectra)
colnames(all_spectra) = unname(multiple_replacement(cols, fluorophore_replacements))
available_fluorophores = sort(unique(
    gsub( '( EM| EX| AB)', '', cols[2:length(cols)] )
))


# ----------------------------------------------------------------------
# Read data

panel <- read_excel_or_csv(file.path(wd, opt[['input-file']]))
fluorophores <- panel['fluorophore']

unavailable_fluorophores = sort(items_in_a_not_b(
   unique(panel[['fluorophore']]), available_fluorophores
))

# save
if (!troubleshooting) {
    if (length(unavailable_fluorophores) > 0) {
        if (!dir.exists(file.path(troubleshooting_dir))) {
            dir.create(file.path(troubleshooting_dir), recursive=TRUE)
        }
        filepath = file.path(troubleshooting_dir, 'unavailable_fluorophores.txt')
        write.table(unavailable_fluorophores, filepath,
                    row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
}


# ----------------------------------------------------------------------
# Plot

spectra_subset <- all_spectra[, c('Wavelength', 'Zombie Aqua EX', 'Zombie Aqua EM')]

fig <- ggplot( pivot_longer(spectra_subset, -Wavelength),
               aes(Wavelength, value, fill = name, colour = name) ) + 
    geom_area(alpha = 0.3) + 
    theme_classic()

if (!troubleshooting) {
    ggsave(file.path(wd, opt[['figures-dir']], 'test.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()


# ----------------------------------------------------------------------
# Troubleshooting only

if (FALSE) {

    if (!dir.exists(file.path(troubleshooting_dir))) {
        dir.create(file.path(troubleshooting_dir), recursive=TRUE)
    }

    # find missing fluorophores
    ab_inv <- read_excel_or_csv(file.path(wd, opt[['antibody-inventory']]))
    ab_inv <- preprocess_antibody_inventory(ab_inv)
    all_fluorophores = sort(unique( ab_inv[['fluorophore']] ))
    unavailable_fluorophores = items_in_a_not_b(
       all_fluorophores, available_fluorophores
    )
    filepath = file.path(troubleshooting_dir, 'all_unavailable_fluorophores.txt')
    write.table(unavailable_fluorophores, filepath,
                row.names = FALSE, col.names = FALSE, quote = FALSE)   

}
