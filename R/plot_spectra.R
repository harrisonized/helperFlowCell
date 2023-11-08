## Plot the spectra for fluorophores used
## Based on this tutorial: https://bradyajohnston.github.io/posts/2022-09-03-plotting-fluorescence/

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

    make_option(c("-c", "--instrument-config"), default='ref/instrument_config.xlsx',
                metavar='ref/instrument_config.xlsx', type="character",
                help="instrument configuration file"),

    make_option(c("-s", "--spectra-file"), default='ref/spectra.csv',
                metavar='ref/spectra.csv', type="character",
                help="spectra data, main source is FPBase"),

    make_option(c("-f", "--figures-dir"), default="figures/output",
                metavar="figures/output", type="character",
                help="set the output directory for the figures"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]
figures_dir <- file.path(wd, opt[['figures-dir']])
troubleshooting_dir = file.path(figures_dir, 'troubleshooting')

# Start Log
start_time = Sys.time()
log <- log_open(paste0("plot_spectra-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Configuration files

# Instrument config
instr_cfg <- read_excel_or_csv(file.path(wd, opt[['instrument-config']]))
instr_cfg <- preprocess_instrument_config(instr_cfg)
instr_cfg_long <- separate_rows(instr_cfg, 'fluorophore', sep=', ')


# Spectra file
all_spectra <- read_excel_or_csv(file.path(wd, opt[['spectra-file']]))

# preprocess column names
colnames(all_spectra) = unname(multiple_replacement(
    colnames(all_spectra), fluorophore_replacements
))
cols <- colnames(all_spectra)
available_fluorophores = sort(unique(
    gsub( '( EM| EX| AB)', '', cols[2:length(cols)] )
))


# ----------------------------------------------------------------------
# Get panel data

panel <- read_excel_or_csv(file.path(wd, opt[['input-file']]))
panel <- merge(
    panel,
    instr_cfg_long[, c('fluorophore', 'laser')],
    by='fluorophore', suffixes=c('', '_'),
    all.x=TRUE, all.y=FALSE, na_matches = 'never'
)
fluorophores <- panel[['fluorophore']]

unavailable_fluorophores = sort(items_in_a_not_b(
   unique(fluorophores), available_fluorophores
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

# Subset and reshape spectra
cols <- colnames(all_spectra)
spectra <- all_spectra[ ,
    c('Wavelength', grep(paste0(fluorophores,collapse = " |"), cols, value = TRUE))
]
spectra_long <- pivot_longer(
    spectra, cols=-Wavelength,
    names_to = "trace_name", values_to = "intensity", values_drop_na=TRUE
)
spectra_long[['fluorophore']] <- gsub(' .{2}$', '', spectra_long[['trace_name']])
spectra_long[['spectrum_type']] <- substr_right(spectra_long[['trace_name']], 2)
spectra_long <- merge(
    spectra_long, panel,
    by='fluorophore', suffixes=c('', '_'),
    all.x=TRUE, all.y=FALSE, na_matches = 'never'
)
spectra_long <- spectra_long[with(spectra_long, order(trace_name, Wavelength)), ]
spectra_long <- reset_index(spectra_long, drop=TRUE)

# Plot each laser separately
for ( laser in unique(panel[['laser']]) ) {

    log_print(paste0('Plottting... ', laser, '.png'))

    fig <- spectra_long[(spectra_long['laser']==laser), ] %>%
        ggplot( aes(x = Wavelength, y = intensity, 
                    fill = fluorophore, group = trace_name) ) + 
        geom_area(
            aes(linetype = spectrum_type),
                position = "identity", 
                alpha = 0.3,
                colour = alpha("black", 0.7)
        ) + 
        scale_linetype_manual(
            values = c(
                "AB" = "dotted", 
                "EX" = "dotted", 
                "EM" = "solid"
            )
        ) + 
        ggtitle(paste(laser, 'Laser')) +
        theme_classic()

    # save
    if (!troubleshooting) {
        ggsave(
            file.path(wd, opt[['figures-dir']],
                      paste0(laser, '.png')),  # filename
            plot=fig,
            height=300, width=1000, dpi=200,
            units="px", scaling=0.5
        )
    }
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
