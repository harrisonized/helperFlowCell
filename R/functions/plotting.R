import::here(zeallot, '%<-%')
import::here(magrittr, '%>%')
import::here(superb, 'showSignificance')
import::here(ggplot2,
    'ggplot', 'aes', 'aes_string', 'geom_boxplot', 'geom_jitter',
    'stat_summary', 'scale_fill_brewer', 'scale_y_continuous', 'expansion',
    'ggtitle')
import::here(file.path(wd, 'R', 'config', 'lasers.R'),
    'color_for_laser', .character_only=TRUE)

## Functions
## compute_nlevels
## plot_spectra_by_each_laser
## plot_violin_with_significance


#' Calculate the number of non-overlapping bars needed to span every significance level
#' Used for plot_violin_with_significance
#' 
compute_nlevels <- function(n) {

    if ((n %% 2) == 0) {
        m <- n/2
        nlevels <- m*(m+1)-m
    } else {
        m <- floor(n/2)
        nlevels <- m*(m+1)
    }
    return(nlevels)
}


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
plot_spectra_by_each_laser <- function(spectra, detectors, laser) {

    excitation <- detectors[(detectors['laser']==laser), ][['excitation']][[1]]
    laser_color <- do.call(switch, c(laser, color_for_laser, "Black"))  # get color
    # c(xmin, xmax) %<-% c(min(spectra[['Wavelength']]), max(spectra[['Wavelength']]))
    c(xmin, xmax) %<-% c(300, 900)  # standardize this

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


#' Plot violin with significance
#' 
#' @description Produces a violin plot with error bars in between
#' You will need to pre-calculate the pvalues and label them as 1v2, 1v3, etc.
#' Do not include any additional columns in the pval_tbl
#' 
#' @return Returns a ggplot object
#'
plot_violin_with_significance <- function(
    df, pval_tbl,
    x='group_name', y='pct_cells',
    title=NULL
) {

    # compute bar positions
    n_groups <- length(unique(df[[x]]))
    nlevels <- compute_nlevels(n_groups)

    id_combos <- setNames(as.data.frame(t(combn(1:n_groups, 2))), c('a', 'b'))
    id_combos[['diff']] <- id_combos[['b']]-id_combos[['a']]
    id_combos <- id_combos[order(id_combos[['diff']]), ]

    # compute starting height
    h_low <- max(aggregate(df[[y]],
        list(df[[x]]), FUN=function(x) mean(x)+sd(x)
    )[['x']])*1.1
    space <- h_low/10

    fig <- ggplot(df, aes_string(x=x, y=y, fill=x)) + 
        geom_boxplot(alpha=0.7, aes_string(middle=y)) +
        stat_summary(fun=mean, geom="crossbar", width=0.75, linewidth=0.25, linetype = "dashed") +
        geom_jitter() +
        scale_fill_brewer(palette="Dark2") +
        scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
        ggtitle(title) +

        # need a loop for arbitrary numbers
        showSignificance( c(1.1,1.9), h_low, -0.005, round(pval_tbl[['pval_1v2']], 4)) +
        showSignificance( c(2.1,2.9), h_low, -0.005, round(pval_tbl[['pval_2v3']], 4)) +
        showSignificance( c(3.1,3.9), h_low, -0.005, round(pval_tbl[['pval_3v4']], 4)) +

        showSignificance( c(1.1,2.9), h_low + space, -0.005, round(pval_tbl[['pval_1v3']], 4)) +
        showSignificance( c(2.1,3.9), h_low + 2*space, -0.005, round(pval_tbl[['pval_2v4']], 4)) +

        showSignificance( c(1.1,3.9), h_low + 3*space, -0.005, round(pval_tbl[['pval_1v4']], 4))

    return(fig)
}
