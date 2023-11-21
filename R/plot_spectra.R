## Plot the spectra for fluorophores used
## Based on this tutorial: https://bradyajohnston.github.io/posts/2022-09-03-plotting-fluorescence/

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
library('tidyr')
library('ggplot2')
library('cowplot')
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
troubleshooting_dir = file.path(wd, 'figures/troubleshooting')

# Start Log
start_time = Sys.time()
log <- log_open(paste0("plot_spectra-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Script-specific variables and functions

lasers = c('UV', 'Violet', 'Blue', 'Green', 'Red')
color_for_laser <- setNames(
    nm=lasers,  # keys
    object=list('Purple', 'Violet', 'Blue', 'Green', 'Red')  # values
)


#' Msain plotting function
#' df is a dataframe
#' detectors is a dataframe
#' laser is chosen from: c('Red', 'Green', 'Blue', 'Violet', 'UV')
#'
plot_spectra_by_each_laser <- function(df, detectors, laser) {

    excitation <- instr_cfg[(instr_cfg['laser']==laser), ][['excitation']][[1]]
    color <- do.call(switch, c(laser, color_for_laser, "Black"))  # get color, default to "Black"

    fig <- df[(df['laser']==laser), ] %>%
        # base plot
        ggplot( aes(x = .data[['Wavelength']], y = .data[['intensity']], 
                    fill = .data[['fluorophore']], group = .data[['trace_name']]),
                show.legend = FALSE ) +
        # plot lasers
        # note: geom_vline doesn't work well here
        geom_rect(
            aes(xmin=excitation-3, xmax=excitation+3, ymin=0, ymax=1),
            fill=color, alpha=0.8,
            inherit.aes = FALSE
        )  +
        # plot detectors
        geom_rect(
            aes(xmin=xmin, xmax=xmax, ymin=0, ymax=1),
            data = detectors[(detectors['laser']==laser), c('xmin', 'xmax')],
            fill="#A3A3A3", alpha=0.6,  # gray
            inherit.aes = FALSE
        ) +
        # fill curves
        geom_area(
            aes(linetype = .data[['spectrum_type']]),
            position = "identity", 
            alpha = 0.3,
            colour = alpha("black", 0.7)
        ) + 
        facet_wrap(vars(.data[['laser']]), strip.position="right") +
        scale_linetype_manual(
            values = c(
                "AB" = "dotted",
                "EX" = "dotted",
                "EM" = "solid"),
            guide = NULL  # disable group as a legend
        ) +
        guides(fill = guide_legend(order=1)) +  # fluorophore
        theme_classic() +
        theme(legend.box = "horizontal",
              legend.justification = c(0, 1)) +  # align legend
        xlim(min(df[['Wavelength']]), max(df[['Wavelength']])) + 
        ylim(0, 1)

    return(fig)
}


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
        filepath = file.path(troubleshooting_dir, 'unavailable_fluorophores.txt')
        write.table(unavailable_fluorophores, filepath,
                    row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
}


# ----------------------------------------------------------------------
# Plot

# select available lasers
lasers <- Reduce(intersect, list( lasers, unique(panel[['laser']])) ) 

# plot all lasers together
plots <- lapply(
    lasers,
    FUN = function(laser) plot_spectra_by_each_laser(spectra_long, instr_cfg, laser)
)

fig <- plot_grid(
    plotlist = plots,
    ncol = 1, nrow = length(lasers),
    align = 'vh', axis = "bt"
)  # note: title is not added, because it causes the legend to become misaligned

if (!troubleshooting) {
    ggsave(
        file.path(wd, opt[['figures-dir']], 'spectra.png'),  # filename
        plot=fig,
        height=1200, width=1200, dpi=200,
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
