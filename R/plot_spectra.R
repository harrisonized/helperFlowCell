## Plot the spectra for fluorophores used
## Based on this tutorial: https://bradyajohnston.github.io/posts/2022-09-03-plotting-fluorescence/

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
library('optparse')
suppressPackageStartupMessages(library('logr'))
import::from(tidyr, 'separate', 'separate_rows', 'pivot_longer')
import::from(cowplot, 'plot_grid')
import::from(ggplot2, 'ggsave')

import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'multiple_replacement', 'items_in_a_not_b', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'read_excel_or_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'reset_index', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'text_tools.R'),
    'substr_right', .character_only=TRUE)
import::from(file.path(wd, 'R', 'config', 'lasers.R'),
    'lasers', .character_only=TRUE)
import::from(file.path(wd, 'R', 'config', 'replacements.R'),
    'fluorophore_replacements', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'preprocess_instrument_config', 'preprocess_antibody_inventory',
    .character_only=TRUE)

source(file.path(wd, 'R', 'functions', 'plotting.R'))  # plot_spectra_by_each_laser


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-file"), default='data/design/panel.csv',
                metavar='data/design/panel.csv', type="character",
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

    make_option(c("-f", "--figures-dir"), default="figures/spectra",
                metavar="figures/spectra", type="character",
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
# Instrument config

instr_cfg <- read_excel_or_csv(file.path(wd, opt[['instrument-config']]))
instr_cfg <- preprocess_instrument_config(instr_cfg)

instr_cfg_long <- separate_rows(instr_cfg, 'fluorophore', sep=', ')

# detector info
instr_cfg <- separate(
    instr_cfg, bandpass_filter, into = c("emission", "range"), '/', convert=TRUE
)
instr_cfg[['xmin']] <- instr_cfg[['emission']] - instr_cfg[['range']] / 2
instr_cfg[['xmax']] <- instr_cfg[['emission']] + instr_cfg[['range']] / 2


# ----------------------------------------------------------------------
# Panel data

panel <- read_excel_or_csv(file.path(wd, opt[['input-file']]))
panel[['fluorophore']] <- multiple_replacement(panel[['fluorophore']], fluorophore_replacements)
panel <- merge(
    panel[(panel['fluorophore']!= ""), ],
    instr_cfg_long[, c('fluorophore', 'laser')],
    by='fluorophore', suffixes=c('', '_'),
    all.x=TRUE, all.y=FALSE, na_matches = 'never'
)
fluorophores <- panel[['fluorophore']]


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


# Subset and reshape spectra
cols <- colnames(all_spectra)
spectra <- all_spectra[ ,
    c('Wavelength', grep(paste0(fluorophores,collapse = " |"), cols, value = TRUE))
]
spectra_long <- pivot_longer(
    spectra,
    cols= -Wavelength,
    names_to = "trace_name", values_to = "intensity",
    values_drop_na=TRUE
)
spectra_long[['fluorophore']] <- gsub(' .{2}$', '', spectra_long[['trace_name']])
spectra_long[['spectrum_type']] <- substr_right(spectra_long[['trace_name']], 2)
spectra_long <- merge(
    spectra_long, panel,
    by='fluorophore', suffixes=c('', '_'),
    all.x=TRUE, all.y=FALSE, na_matches = 'never'
)
spectra_long <- spectra_long[with(spectra_long, order(trace_name, Wavelength)), ]
spectra_long <- spectra_long[(spectra_long[['fluorophore']] %in% fluorophores), ]  # keep relevant fluorophores only
spectra_long <- reset_index(spectra_long, drop=TRUE)

# fluorophores in panel not found in the spectra file
if (!troubleshooting) {
    unavailable_fluorophores = sort(items_in_a_not_b(
       unique(fluorophores), available_fluorophores
    ))

    if (length(unavailable_fluorophores) > 0) {
        if (!dir.exists(file.path(troubleshooting_dir))) {
            dir.create(file.path(troubleshooting_dir), recursive=TRUE)
        }
        filepath = file.path(troubleshooting_dir, 'unplotted_fluorophores.txt')
        write.table(unavailable_fluorophores, filepath,
                    row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
}


# ----------------------------------------------------------------------
# Plot

# TODO: bugfix for when a fluorophore is available in spectra but not in instr_cfg

# select available lasers
lasers <- Reduce(intersect, list( lasers, unique(panel[['laser']])) ) 
plots <- lapply(lasers, FUN = function(laser)
    plot_spectra_by_each_laser(spectra_long, instr_cfg, laser)
)

# plot all lasers together
# note: title is not added, because it causes the legend to become misaligned
fig <- plot_grid(plotlist = plots,
                 ncol = 1, nrow = length(lasers),
                 align = 'vh', axis = "bt")

# save
if (!troubleshooting) {
    if (!dir.exists(file.path(wd, opt[['figures-dir']]))) {
        dir.create(file.path(wd, opt[['figures-dir']]), recursive=TRUE)
    }

    ggsave(
        file.path(wd, opt[['figures-dir']], 'spectra.png'),  # filename
        plot=fig,
        height=200+200*length(lasers), width=1200, dpi=200,
        units="px", scaling=0.5
    )
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
