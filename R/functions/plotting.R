library('zeallot')
library('ggplot2')

#' Main plotting function for plot_spectra
#' spectra is a dataframe
#' detectors is a dataframe
#' laser is chosen from: c('Red', 'Green', 'Blue', 'Violet', 'UV')
#'
plot_spectra_by_each_laser <- function(spectra, detectors, laser) {

    excitation <- detectors[(detectors['laser']==laser), ][['excitation']][[1]]
    laser_color <- do.call(switch, c(laser, color_for_laser, "Black"))  # get color
    c(xmin, xmax) %<-% c(min(spectra[['Wavelength']]), max(spectra[['Wavelength']]))

    fig <- spectra[(spectra['laser']==laser), ] %>%
        # base plot
        ggplot( aes(x = .data[['Wavelength']], y = .data[['intensity']], 
                    fill = .data[['fluorophore']], group = .data[['trace_name']]),
                show.legend = FALSE ) +
        # plot lasers
        # note: geom_vline doesn't work well here
        geom_rect(
            aes(xmin=excitation-3, xmax=excitation+3, ymin=0, ymax=1),
            fill=laser_color, alpha=0.8,
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
        xlim(xmin, xmax) + 
        ylim(0, 1)

    return(fig)
}
