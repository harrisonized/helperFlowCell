import::here(dplyr, 'vars')
import::here(zeallot, '%<-%')
import::here(magrittr, '%>%')
import::here(ggplot2,
    'ggplot', 'aes', 'alpha', 'facet_wrap',
    'geom_rect', 'geom_area', 'geom_segment', 
    'element_blank', 'element_text', 'unit',
    'scale_linetype_manual', 'scale_fill_discrete', 'scale_fill_gradientn',
    'guides', 'guide_legend', 'theme', 'theme_classic',
    'xlim', 'ylim', 'labs')
import::here(scales, 'rescale')
import::here(file.path(wd, 'R', 'config', 'lasers.R'),
    'color_for_laser', .character_only=TRUE)

## Functions
## plot_spectra_by_each_laser


#' Individual Spill Plots
#' 
#' @description Produces a spill plot overlayed with the instrument configuration.
#' These can be stacked using cowplot to produce a spill plot faceted by laser.
#' 
#' @param spectra long-format dataframe of emission/excitation data
#' @param detectors instrument configuration, requiring a column for bandpass in the format of 'center/width', eg. '530/30'
#' @param laser choose from c('Red', 'Green', 'Blue', 'Violet', 'UV')
#' @return Returns a ggplot object
#'
plot_spectra_by_each_laser <- function(spectra, detectors, laser, xmin=300, xmax=900) {

    laser_spectra <- spectra[(spectra[['laser']] == laser) & (!is.na(spectra[['laser']])), ]
    laser_detectors <- detectors[(detectors[['laser']] == laser) & (!is.na(detectors[['laser']])), ]
    excitation <- laser_detectors[['excitation']][[1]]
    laser_color <- do.call(switch, c(laser, color_for_laser, "Black"))  # get color

    fig <- ggplot() +
        # shade detectors
        geom_rect(data=laser_detectors,
            aes(xmin=xmin, xmax=xmax, ymin=0, ymax=1),
            fill="#A3A3A3", alpha=0.6,  # gray
            inherit.aes = FALSE
        ) + { # add detector boundaries
            edges <- list(
                list("xmin", "xmin", 0, 1),  # left
                list("xmin", "xmax", 1, 1),  # top
                list("xmax", "xmax", 0, 1)   # right
            )
            lapply(edges, function(edge) {
                geom_segment(
                    data = laser_detectors,
                    aes(x = .data[[edge[[1]]]], xend = .data[[edge[[2]]]],
                        y = edge[[3]], yend = edge[[4]] ),
                    colour = "#595959",  # dark grey
                    linewidth = 0.5,
                    inherit.aes = FALSE
                )
            })
        } +
        facet_wrap(vars(laser), strip.position = "right") +
        theme_classic() +
        xlim(xmin, xmax) +
        ylim(0, 1) +
        labs(y = "intensity")

    # Add spectra
    if (nrow(laser_spectra) > 0) {

        # Calculate max emission wavelength
        is_em <- grepl("EM$", laser_spectra[["trace_name"]])
        em_indices <- split(
            which(is_em),
            laser_spectra[["fluorophore"]][is_em]
        )
        max_ems <- vapply(em_indices, function(i) {
            laser_spectra[["Wavelength"]][i][
                which.max(laser_spectra[["intensity"]][i])
            ]
        }, numeric(1))
        max_ems <- sort(max_ems)
        laser_spectra[["peak_wavelength"]] <-
            max_ems[laser_spectra[["fluorophore"]]]

        fig <- fig +
            geom_area(
                data = subset(laser_spectra, spectrum_type %in% c("AB", "EX")),
                aes(x = .data[['Wavelength']], y = .data[['intensity']],
                    fill = .data[['peak_wavelength']], group = .data[['trace_name']],
                    linetype = .data[['spectrum_type']]),
                position = "identity",
                alpha = 0.08,
                colour = alpha("black", 0.6)
            ) +
            geom_area(
                data = subset(laser_spectra, spectrum_type == "EM"),
                aes(x = .data[['Wavelength']], y = .data[['intensity']],
                    fill = .data[['peak_wavelength']], group = .data[['trace_name']],
                    linetype = .data[['spectrum_type']]),
                position = "identity",
                alpha = 0.36,
                colour = alpha("black", 0.6)
            ) +
            scale_fill_gradientn(
                limits = c(280, 850),
                breaks = unname(max_ems),
                labels = names(max_ems),
                name = NULL,  # "Fluorophore",
                colours = c(
                    "#4B0082",  # indigo
                    "#8000FF",  # vivid purple
                    "#3b82f6",  # blue 500
                    "#17becf",  # cyan
                    "#2ca02c",  # green
                    "#FFFF00",  # yellow
                    "#ff7f0e",  # orange
                    "#ff7373",  # red
                    "#5C0000"   # deep red
                ),
                values = scales::rescale(c(280, 405, 450, 496, 570, 615, 661, 750, 850)),
                guide = "none") +
            scale_linetype_manual(
                values = c(
                    "AB" = "dotted",
                    "EX" = "dotted",
                    "EM" = "solid"),
                guide = NULL  # disable group as a legend
            ) +
            guides(fill = guide_legend(order = 1)) + 
            theme(
                legend.box = "horizontal",
                legend.justification = c(0, 1),
                legend.text = element_text(size = 8),
                legend.key.size = unit(14, "pt"),
                legend.key.width = unit(14, "pt"),
                legend.key.height = unit(14, "pt"),
            )
    }

    # plot laser excitation
    # note: geom_vline doesn't work well here
    fig <- fig + geom_rect(
            data=data.frame(xmin=excitation-3, xmax=excitation + 3, ymin=0, ymax=1),
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill=laser_color,
            alpha=0.8,
            inherit.aes=FALSE)

    return(fig)
}
