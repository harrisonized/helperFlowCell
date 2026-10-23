# import::here(rlang, 'sym')
import::here(magrittr, '%>%')
import::here(dplyr, 'group_by', 'summarize', 'summarise', 'mutate',
    'reframe', 'first', 'n_distinct', 'arrange')
import::here(tidyr, 'pivot_wider', 'unnest')
import::here(tibble, 'tibble')
import::here(ggplot2,
    'ggplot', 'aes', 'labs', 'theme', 'theme_minimal', 'margin',
    'geom_boxplot', 'geom_jitter', 'geom_col', 'geom_line', 'geom_area',
    'element_text', 'stat_summary', 'scale_fill_brewer', 'scale_fill_manual',
    'scale_color_manual', 'scale_x_continuous', 'scale_x_discrete', 'scale_y_continuous',
    'ylim', 'guides', 'guide_legend', 'guide_axis', 'expansion', 'ggtitle')
import::here(cowplot, 'plot_grid', 'get_legend', 'get_plot_component')
import::here(flowCore, 'logicleTransform')
import::here(RColorBrewer, 'brewer.pal')
import::here(ggprism, 'theme_prism')
import::here(superb, 'showSignificance')
import::here(plotly, 'plot_ly', 'add_trace', 'layout', 'save_image')
import::here(htmlwidgets, 'saveWidget')  # brew install pandoc

import::here(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'dropna', 'separate_groups_with_na', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'text_tools.R'),
    'format_number', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'math.R'),
    'unpaired_t_test', 'fishers_lsd', 'tukey_multiple_comparisons',
    'bonferroni_multiple_comparisons', 'generate_lognormal_data',
    .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'math.R'),
    'apply_unpaired_t_test', .character_only=TRUE)
import::here(file.path(wd, 'R', 'config', 'colors.R'),
    'mpl_colors', 'tailwind_colors', .character_only=TRUE)


## Functions
## save_fig
## plot_dots
## plot_scatter
## plot_violin
## compute_nlevels
## generate_base_level
## get_significance_code
## plot_multiple_comparisons
## plot_modal_histograms
## plot_polar_area
## plot_ribbon_chart


#' Save Figure
#'
#' @export
save_fig <- function(fig,
    dirpath='figures',
    filename='fig',
    height=500, width=800, scale=3,
    save_html=FALSE
) {

    if (!dir.exists( file.path(dirpath, 'png') )) {
        dir.create( file.path(dirpath, 'png'), recursive=TRUE)
    }

    # reticulate::py_install(
    #     packages = "kaleido==0.2.1",
    #     envname = "r-reticulate",
    #     method = "virtualenv",
    #     pip = TRUE
    # )
    suppressWarnings(save_image(fig,
        file=file.path(path.expand(dirpath), 'png', paste0(filename, '.png')), 
        height=height, width=width, scale=scale
    ))

    # save HTML
    if (save_html) {
        if (!dir.exists( file.path(dirpath, 'html') )) {
            dir.create( file.path(dirpath, 'html'), recursive=TRUE)
        }
        suppressWarnings(saveWidget(
            widget = fig,
            file=file.path(file.path(dirpath, 'html'), paste0(filename, '.html')),
            selfcontained = TRUE
        ))
        unlink(
            file.path(file.path(dirpath, 'html'), paste0(filename, '_files')),
            recursive=TRUE
        )            
    }
}


#' Plot Dots
#'
#' @description  Thin wrapper around ggplot geom_jitter
#'
plot_dots <- function(df,
    x=NULL, y=NULL, color=NULL,
    xlabel=NULL, ylabel=NULL, title=NULL
) {

    fig <- ggplot(
        df,
        aes(x=reorder(.data[[x]], .data[[y]], decreasing=TRUE),
            y=.data[[y]]), na.rm=TRUE) +
        geom_jitter(aes(colour=.data[[color]])) +
        labs(x=xlabel, y=ylabel, title=title) +
        theme(axis.text.x = element_text(angle = 45, hjust=1))

    return(fig)
}


#' Plot Scatter
#'
#' @description Thin wrapper around Plotly scatter plot
#' 
plot_scatter <- function(
    df, x, y, group_by=NULL, separate_by=NULL, size=NULL,
    xlabel=NULL, ylabel=NULL, title=NULL,
    xmin=NULL, xmax=NULL,
    ymin=NULL, ymax=NULL,
    yaxis_type='linear',
    color='#1f77b4',
    color_discrete_map=NULL,
    hover_data=c(),
    mode='markers'
) {

    if (is.null(ymin)) {
        xmin <- min(df[[x]])
    }
    if (is.null(ymax)) {
        xmax <- max(df[[x]])
    }
    xrange = c(xmin - (xmax-xmin) * 0.1, xmax + (xmax-xmin) * 0.1)  # add padding
    if (is.null(ymin)) {
        ymin <- min(df[[y]])
    }
    if (is.null(ymax)) {
        ymax <- max(df[[y]])
    }
    if (yaxis_type=='linear') {
        yrange = c(ymin - (ymax-ymin) * 0.1, ymax + (ymax-ymin) * 0.1)  # add padding
    }
    if (yaxis_type=='log') {
        yrange = sapply(c(ymin, ymax), function(y) if (y > 0) {log(y, base=10)} else {NULL} )
    }

    if (!is.null(group_by) & !is.null(separate_by)) {
        df <- separate_groups_with_na(df, id=separate_by, order_by=c(separate_by, x), na_cols=c(x, y))
    }

    # color
    if (!is.null(color_discrete_map)) {
        color_discrete_map <- list2env(as.list(color_discrete_map))
    }

    fig <- plot_ly(type='scatter', mode=mode)

    if (is.null(group_by)) {

        hovertext <- ''
        fig <- fig %>%
            add_trace(
                x = df[[x]],
                y = df[[y]],
                color=I(color),
                marker = list(size = 5),
                hovertext = hovertext
            )
    } else {
        for (group in unique(df[[group_by]])) {

            # hoverdata
            hovertext <- ''
            for (field in c(x, y, hover_data)) {
                if (!is.null(field)) {
                    hovertext <- paste0(hovertext, field, "=", df[(df[[group_by]] == group), field], "<br>")
                }
            }
            color <- color_discrete_map[[group]]
            if (!is.null(color)) {
                color <- I(color)
            }

            fig <- fig %>%
                add_trace(
                    x = unlist(df[(df[[group_by]] == group), x]),
                    y = unlist(df[(df[[group_by]] == group), y]),
                    color = color,
                    legendgroup = group,
                    name = group,
                    marker = list(size = 5),
                    hovertext = hovertext
                )
        }
    }

    fig <- fig %>% layout(
        title = list(
            text = title,
            x = 0.5
        ),
        xaxis = list(
            title_text = xlabel,
            range = xrange
        ),
        yaxis = list(
            title_text = ylabel,
            showgrid = TRUE, gridcolor = '#E4EAF2', zeroline = FALSE,
            range = yrange,
            type = yaxis_type
        ),
        plot_bgcolor = 'rgba(0,0,0,0)',
        showlegend = ifelse(is.null(group_by), FALSE, TRUE),
        hovermode = 'closest'

    )

    return(fig)
}


#' Plot Violin
#'
#' @description Thin wrapper arounbd Plotly violin plot
#' 
plot_violin <- function(
    df, x, y, group_by=NULL, size=NULL,
    xlabel=NULL, ylabel=NULL, title=NULL,
    ymin=NULL, ymax=NULL,
    legend_order=NULL,
    xaxis_angle=45,
    yaxis_type='linear',
    color='#1f77b4',
    color_discrete_map=NULL,
    hover_data=c(),
    sort=TRUE,
    descending=TRUE,
    violinmode='group'  # overlay
) {

    if (is.null(ymin)) {
        ymin <- min(df[[y]])
    }
    if (is.null(ymax)) {
        ymax <- max(df[[y]])
    }
    if (yaxis_type=='linear') {
        yrange = c(ymin - (ymax-ymin) * 0.1, ymax + (ymax-ymin) * 0.1)  # add padding
    }
    if (yaxis_type=='log') {
        yrange = sapply(c(ymin, ymax), function(y) if (y > 0) {log(y, base=10)} else {NULL} )
    }

    if (sort) {
        if (descending) {
            x_axis_order <- df %>%
                group_by(.data[[x]]) %>%
                summarize(median = median(.data[[y]], na.rm=TRUE)) %>%
                arrange(-.data[['median']]) %>%
                select(.data[[x]])
        } else {
            x_axis_order <- df %>%
                group_by(.data[[x]]) %>%
                summarize(median = median(.data[[y]], na.rm=TRUE)) %>%
                arrange(.data[['median']]) %>%
                select(.data[[x]])
        }
    }

    # color
    if (!is.null(color_discrete_map)) {
        color_discrete_map <- list2env(as.list(color_discrete_map))
    }

    fig <- plot_ly(type = 'violin')


    # ----------------------------------------------------------------------
    # No groups
    
    if (is.null(group_by)) {
        if (sort) {
            df <- df %>% arrange(match(.data[[x]], unlist(x_axis_order)))
        }

        # hoverdata
        hovertext <- ''
        for (field in c(x, y, hover_data)) {
            if (!is.null(field)) {
                hovertext <- paste0(hovertext, field, "=", df[[field]], "<br>")
            }
        }

        fig <- fig %>%
            add_trace(
                x = df[[x]],
                y = df[[y]],
                color=I(color),
                box = list(visible = TRUE),
                meanline = list(visible = TRUE),
                points = 'all',
                jitter = 0.2,
                pointpos = -1,
                marker = list(size = 5),
                hovertext = hovertext
            )

        return(fig)
    }


    # ----------------------------------------------------------------------
    # Group by

    if (is.null(legend_order)) {
        legend_order <- sort(unique(df[[group_by]]))
    }
    for (group in legend_order) {

        if (sort) {
            df <- df %>% arrange(
                match(.data[[x]], unlist(x_axis_order))
            )
        }

        # hoverdata
        hovertext <- ''
        for (field in c(x, y, hover_data)) {
            if (!is.null(field)) {
                hovertext <- paste0(hovertext, field, "=", df[(df[[group_by]] == group), field], "<br>")
            }
        }
        color <- color_discrete_map[[group]]
        if (!is.null(color)) {
            color <- I(color)
        }

        fig <- fig %>%
            add_trace(
                x = unlist(df[(df[[group_by]] == group), x]),
                y = unlist(df[(df[[group_by]] == group), y]),
                color = color,
                legendgroup = group,
                scalegroup = group,
                name = group,
                box = list(visible = TRUE),
                meanline = list(visible = TRUE),
                points = 'all',
                jitter = 0.2,
                pointpos = -1,
                marker = list(size = 5),
                hovertext = hovertext
            )
    }
    
    fig <- fig %>% layout(
        title = list(
            text = title,
            x = 0.5
        ),
        xaxis = list(title_text = xlabel, tickangle=xaxis_angle),
        yaxis = list(
            title_text = ylabel,
            showgrid = TRUE, gridcolor = '#E4EAF2', zeroline = FALSE,
            range = yrange,
            type = yaxis_type
        ),
        violinmode = violinmode,
        plot_bgcolor = 'rgba(0,0,0,0)',
        showlegend = ifelse(is.null(group_by), FALSE, TRUE),
        hovermode = 'closest'

    )

    return(fig)
}


#' Compute N Levels
#' 
#' @description Calculate the number of non-overlapping bars needed to span
#' every significance level. Didn't end up needing this.
#' Helper function for plot_multiple_comparisons
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


#' Generate Starting Positions
#'
#' @description Generate a sequence used as a floor for each particular level
#' This serves as the correction factor for the bar height of each p_value bar
#' Helper function for plot_multiple_comparisons
#' 
generate_base_level <- function(n_groups) {

    if ( n_groups == 1 ) {
        return(0)
    } else if ( n_groups == 2 ) {
        return(1)
    }

    # compute number of levels within the group
    mid <- floor(n_groups/2)
    if ( (n_groups%%2) == 0 ) {
        nlevels <- c(1:mid, (mid-1):1)  # even
    } else {
        nlevels <- c(1:mid, (mid):1)  # odd
    }

    res <- c()
    value <- 1
    nreps <- n_groups-1
    for (step in nlevels) {
        res <- c(res, rep(value, nreps))
        value <- value + step
        nreps <- nreps - 1
    }
    return(res)
}


#' Get Significance Code
#' 
#' @description Convert p value to significance
#' Helper function for plot_multiple_comparisons
#' Note: Adobe Illustrator may throw the following error:
#' "Certain alternate glyphs are not available in this version of Illustrator.""
#' To circumvent this, first copy the svg into Powerpoint, then from Powerpoint to Adobe Illustrator
#' 
get_significance_code <- function(p, digits=3) {

    if (is.na(p)) {
        return('NA')
    }

    p <- abs(p)
    if (p > 0.20) {
        return('ns')
    } else if (p > 0.05) {
        return( format(round(p, digits), nsmall = digits) )
    } else if (p > 0.01) {
        return('✱')  # U+2731 asterisk
    } else if (p > 0.001) {
        return('✱✱')
    } else if (p > 0.0001 ) {
        return('✱✱✱')
    } else {
        return('✱✱✱✱')
    }
}


#' Plot violin with significance
#' 
#' @description Produces a violin plot with error bars in between
#' You will need to pre-calculate the pvalues and label them as 1v2, 1v3, etc.
#' Do not include any additional columns in the pvals object
#' Will update this function to calculate pvals within the function
#' 
#' @return Returns a ggplot object
#'
plot_multiple_comparisons <- function(
    df,
    x,  # 'group_name'
    y,  #'pct_cells'
    xlabel=NULL, ylabel=NULL, title=NULL,
    ymin=0,
    xaxis_angle=60,
    test='t_test',  # 'fishers_lsd', 't_test', 'tukey', or 'bonferroni'
    show_numbers=FALSE,
    digits=4,
    show_brackets=TRUE,
    custom_group_order=c()
) {

    df <- df[(!is.na(df[[y]])), ]

    if (length(custom_group_order)>=1) {
        group_names <- intersect( custom_group_order, unique(df[[x]]) )
        df <- df[(df[[x]] %in% group_names), ]
        df[[x]] <- factor(df[[x]], levels = group_names)
    } else {
        group_names <- sort(unique( df[[x]] ))
    }

    n_groups <- length(group_names)

    if (show_brackets) {

        # compute pvals
        if (test=='t_test') {
            pvals <- apply_unpaired_t_test(
                df,
                index_cols=items_in_a_not_b(colnames(df), c(x, y)),
                group_name=x,
                metric=y,
                custom_group_order=custom_group_order
            )
        } else if (test=='fishers_lsd') {
            pvals <- fishers_lsd(
                df,
                group=x,
                metric=y,
                custom_group_order=custom_group_order
            )
        } else if (test=='tukey') {
            pvals <- tukey_multiple_comparisons(
                df,
                group=x,
                metric=y,
                custom_group_order=custom_group_order
            )
        } else if (test=='bonferroni') {
            pvals <- bonferroni_multiple_comparisons(
                df,
                group=x,
                metric=y,
                custom_group_order=custom_group_order
            )
        } else {
            stop("Choose from test='t_test', 'fishers_lsd', 'tukey', or 'bonferroni'")
        }

        # compute bar positions
        if (n_groups > 1) {
            bracket_params <- setNames(as.data.frame(t(combn(1:n_groups, 2))), c('left', 'right'))
            bracket_params[['dist']] <- bracket_params[['right']] - bracket_params[['left']]
            bracket_params <- bracket_params[ order(bracket_params[['dist']]), ]
            bracket_params[['base_level']] <- generate_base_level(n_groups)
            bracket_params[['level']] <- bracket_params[['base_level']] +  # base level
                (bracket_params[['left']]-1) %% bracket_params[['dist']]  # alternating correction factor

            # compute starting height
            tryCatch({
                withCallingHandlers({
                    h_low <- max(
                        aggregate(df[[y]], list(df[[x]]),
                        FUN=function(x) mean(x)+sd(x)
                    )[['x']], na.rm=TRUE) * 1.1
                }, warning = function(w) {
                    if ( grepl("no non-missing arguments", w) ) {
                        stop()
                    }
                })
            },
            error = function(condition) {
                h_low <<- max(df[[y]], na.rm=TRUE) * 1.1
            })
            if (h_low < max(df[[y]], na.rm=TRUE)) {
                h_low <- max(df[[y]], na.rm=TRUE) * 1.1
            }
            space <- h_low / 10
        }
    }

    # base plot
    fig <- ggplot(df, aes(x=.data[[x]], y=.data[[y]], fill=.data[[x]])) + 
        geom_boxplot(alpha=0.7, aes(middle=.data[[y]])) +
        stat_summary(fun=mean, geom="crossbar", width=0.75, linewidth=0.25, linetype = "dashed") +
        geom_jitter() +
        (if (n_groups <=8 ) scale_fill_brewer(palette="Dark2")
            else scale_fill_manual(
                values = colorRampPalette(brewer.pal(8, "Dark2"))(n_groups)
            )
        ) +
        scale_x_discrete(guide=guide_axis(angle=xaxis_angle)) +
        scale_y_continuous(
            limits = c(ymin, NA),
            expand = expansion(mult = c(0, 0.1)),
            labels = function(x) format(x, scientific=FALSE)
        ) +
        labs(x=xlabel, y=ylabel, title=title) +
        theme_prism(base_size = round(50/n_groups, 0))

    if (show_brackets) {
        # significance brackets
        if (n_groups > 1) {
            for (row in 1:nrow(bracket_params)) {

                left <- bracket_params[row, "left"]
                right <- bracket_params[row, "right"]
                level <- bracket_params[row, "level"]
                colname <- paste(group_names[[right]], group_names[[left]], sep='-')  # colname

                if (show_numbers) {
                    pval <- toString( format(round(pvals[[colname]], digits), nsmall = digits) )
                } else {
                    pval <- get_significance_code( pvals[[colname]] )
                }

                if (h_low >= 0) {
                    fig <- fig +
                        showSignificance(
                            x=c(left+0.1, right-0.1), y=h_low+(level-1)*space, width=-0.001*h_low,
                            text=pval, textParams=list(size=(if (n_groups <= 6) 3 else 2))
                        )
                }
                if (h_low < 0) {
                    fig <- fig +
                        showSignificance(
                            x=c(left+0.1, right-0.1), y=h_low-(level+1)*space, width=-0.001*h_low,
                            text=pval, textParams=list(size=(if (n_groups <= 6) 3 else 2))
                        )
                }

            }
        }
    }

    return(fig)
}


#' Plot Modal Histograms
#' 
#' Mimics Flowjo's histogram plots
#' Use [generate_gaussian_data()] to generate the input
#' 
plot_modal_histograms <- function(df,
    group='group', value='value',
    xlabel=NULL, ylabel=NULL, title=NULL,
    x_ticks = c(10, 250, 500, 1000, 2000, 4000, 8000),
    colors=c(),  # c("mNeonGreen"="#ACF53B", "gray"="#B0B0B0")
    nbins=100,
    spar=0.33,  # smoothing parameter, 0.33 looks good
    show_bins=FALSE,
    max_scale=8289.7211,  # top of the scale data value
    pos_decades=4.5,  # logicleTransform m param, full width of transformed display
    lin_width=2.25,  # logicleTransform w param, linearization width
    extra_neg=0,  # only works if pos_decades - 2*lin_width > 0
    width_basis=0.25  # determines how much of the left side is shown
) {

    # Clamp lin_width to valid range
    if (lin_width < 0) {
        lin_width <- 0
        warning(sprintf("lin_width (%.2f) is negative; clamping to 0", lin_width))
    }
    if (lin_width == pos_decades/2) {
        lin_width <- pos_decades/2 - 1e-6
    }
    if (lin_width > pos_decades/2) {
        lin_width <- pos_decades/2 - 1e-6
        warning(sprintf("lin_width (%.2f) is too large; clamping to %.2f", lin_width, pos_decades/2))
    }

    # Clamp extra_neg to valid range
    if (extra_neg < 0) {
        extra_neg <- 0
        warning(sprintf("extra_neg (%.2f) is negative; clamping to 0", lin_width))
    }
    max_extra_neg <- pos_decades - 2 * lin_width
    if (extra_neg > max_extra_neg) {
        warning(sprintf("extra_neg (%.2f) is too large; clamping to %.2f", extra_neg, max_extra_neg))
        extra_neg <- max_extra_neg
    }

    # Adjust scale
    raw_min <- max_scale / 10^(pos_decades + extra_neg)  # default positive min
    if (min(x_ticks) < raw_min) {
        raw_min <- min(x_ticks)
    }
    if (max(x_ticks) > max_scale) {
        max_scale <- max(x_ticks)
    } 

    # Logicle transform the data
    logicle <- logicleTransform(w=lin_width, t=max_scale, m=pos_decades, a=extra_neg)  # default options
    df[['logicle_value']] <- logicle(df[[value]])

    # Range calculations
    transformed_min <- logicle(raw_min)
    transformed_max <- logicle(max_scale)
    pad <- width_basis * (transformed_max - transformed_min)  # add padding to min side
    transformed_range <- c(transformed_min - pad, transformed_max)

    # Compute clamped histogram
    breaks <- pretty(transformed_range, n = nbins)
    hist_data <- df %>%
        group_by(group) %>%
        summarise(
            counts = list({
                x_clamped <- pmin(pmax(logicle_value, min(breaks)), max(breaks))
                hist(x_clamped, breaks = breaks, plot = FALSE)$counts
            }),
            mids = list({
                x_clamped <- pmin(pmax(logicle_value, min(breaks)), max(breaks))
                hist(x_clamped, breaks = breaks, plot = FALSE)$mids
            })
        ) %>%
        unnest(cols = c(counts, mids)) %>%
        group_by(group) %>%
        mutate(norm_counts = counts / max(counts))

    # Smoothed and clamped interpolation
    interp_data <- hist_data %>%
        group_by(group) %>%
        reframe({
            fit <- smooth.spline(mids, norm_counts, spar = spar)
            x_vals <- seq(min(mids), max(mids), length.out = 1000)
            y_vals <- predict(fit, x = x_vals)$y
            tibble(
                x = pmin(pmax(x_vals, transformed_range[1]), transformed_range[2]),
                y = pmax(y_vals, 0),
                group = first(group)
            )
        })

    # Generate axis breaks
    if (is.null(x_ticks)) {
        raw_breaks <- 10^seq(floor(log10(abs(raw_min))), ceiling(log10(max_scale)), by = 1)
        raw_breaks <- raw_breaks[(raw_breaks >= raw_min) & (raw_breaks <= max_scale)]
    } else {
        raw_breaks <- x_ticks[(x_ticks >= raw_min) & (x_ticks <= max_scale)]
    }
    transformed_breaks <- logicle(raw_breaks)

    # Expand colors if needed
    n_groups <- n_distinct(df[[group]])
    n_colors <- length(colors)
    if ((n_colors > 0) & (n_colors < n_groups)) {
        colors <- c(colors, rep(colors[length(colors)], n_groups-n_colors))
    }
    if (n_colors > n_groups) {
        colors <- colors[1:n_groups]
    }

    # set group orders
    df[[group]] <- factor(df[[group]], levels = unique(df[[group]]))
    hist_data$group <- factor(hist_data$group, levels = levels(df[[group]]))
    interp_data$group <- factor(interp_data$group, levels = levels(df[[group]]))
    hist_data <- hist_data %>% arrange(group)
    interp_data <- interp_data %>% arrange(group)

    # Plot
    fig <- ggplot()

    groups_ordered <- levels(df[[group]])  # first in front, last in back
    for (g in groups_ordered) {
        hist_sub <- hist_data[(hist_data[[group]]==g), ]
        interp_sub <- interp_data[(interp_data[[group]]==g), ]
        fig <- fig +
            (if (show_bins) {
                # histograms
                geom_col(data = hist_data[(hist_data[[group]]==g), ],
                         aes(x = mids, y = norm_counts, fill = group),
                         position = "identity", alpha = 0.4, width = diff(breaks)[1],
                         inherit.aes=FALSE)
            } else {
                # shaded area
                geom_area(data = interp_data[(interp_data[[group]]==g), ],
                          aes(x = x, y = y, fill = group),
                          alpha = 0.4, position = "identity", inherit.aes=FALSE)
            }) +
            geom_line(data = interp_data[(interp_data[[group]]==g), ],
                      aes(x = x, y = y, color = group),
                      linewidth = 1.2, inherit.aes=FALSE)  # boundary
    }

    fig <- fig +
        (if (length(colors) == n_groups) { scale_fill_manual(values = colors) } else NULL) +  # area
        (if (length(colors) == n_groups) { scale_color_manual(values = colors) } else NULL) +  # line
        scale_x_continuous(limits = transformed_range, breaks = transformed_breaks, labels = raw_breaks) +
        ylim(0, 1) +
        labs(title = title, x = xlabel, y = ylabel) +
        theme_minimal(base_size = 12)

    # split then recombine with define widths
    # legend <- get_plot_component(fig, "guide-box", return_all=TRUE)[[1]]  # silence warning
    legend <- get_legend(fig)  # silence warning
    fig <- fig + theme(legend.position="none")
    fig <- plot_grid(fig, legend, ncol = 2, rel_widths = c(3, 1))

    return(fig)
}


#' Plot Polar Area
#' 
#' Angle represents the number of degrees, radius represents size
#' 
#' Plot Polar Area
#'
#' Angle represents the number of degrees, radius represents size
#'
plot_polar_area <- function(df,
    x = "population_pct", y = "radius_value", group = "category",
    title = NULL,
    show_axes = FALSE,
    colors = NULL,
    normalize = TRUE,
    min_label_pct = 2,  # Minimum percentage to show label
    min_radius = 0.3,   # Minimum radius to show label inside
    min_display_pct = 0 # Minimum percentage for display
) {

    # Default Tailwind palette for categories without an assigned color
    tailwind_colors <- c(
        "#3b82f6",  # blue 500
        "#8b5cf6",  # violet 500
        "#ec4899",  # pink 500
        "#f59e0b",  # amber 500
        "#10b981",  # emerald 500
        "#ef4444",  # red 500
        "#06b6d4",  # cyan 500
        "#f97316",  # orange 500
        "#6366f1",  # indigo 500
        "#14b8a6",  # teal 500
        "#a855f7"   # purple 500
    )

    # Complete the color mapping for categories in the data
    if (is.null(colors)) {
        colors <- list()
    } else {
        colors <- as.list(colors)
    }

    categories <- unique(as.character(df[[group]]))
    categories <- categories[!is.na(categories)]

    # Find categories without explicitly assigned colors
    missing_categories <- sort(setdiff(categories, names(colors)))

    # Assign fallback colors to missing categories
    if (length(missing_categories) > 0) {
        fallback_colors <- rep(
            tailwind_colors,
            length.out = length(missing_categories)
        )

        names(fallback_colors) <- missing_categories

        colors <- c(colors, as.list(fallback_colors))
    }

    # Sort by radius value
    df <- df[order(-df[[y]]), ]

    # Normalize to 100%
    total_pct <- sum(df[[x]])
    if ((normalize && total_pct != 100) || total_pct > 100) {
        df[[x]] <- df[[x]] * 100 / total_pct
    }

    # Adjust small percentages to minimum display threshold
    small_mask <- df[[x]] < min_display_pct & df[[x]] > 0

    if (any(small_mask)) {
        # Total deficit needed to boost small slices
        deficit <- sum(min_display_pct - df[[x]][small_mask])

        # Identify categories that are not being adjusted upward
        large_mask <- !small_mask

        # Total percentage available from the remaining categories
        remaining_total <- sum(df[[x]][large_mask])

        if (remaining_total > 0) {
            # Boost small slices to the minimum
            df[[x]][small_mask] <- min_display_pct

            # Shrink remaining categories proportionally
            df[[x]][large_mask] <-
                df[[x]][large_mask] *
                (remaining_total - deficit) / remaining_total
        }
    }

    df[["angle_degrees"]] <- df[[x]] * (360 / 100)
    df[["theta_start"]] <- c(
        0,
        cumsum(df[["angle_degrees"]])[-nrow(df)]
    )
    df[["theta_center"]] <-
        df[["theta_start"]] + df[["angle_degrees"]] / 2

    # Create labels with truncated values (3 significant figures)
    df[["label"]] <- paste0(
        sprintf("%.3g", df[[x]]), "%\n",
        sprintf("%.3g", df[[y]])
    )

    # Fixed range from -0.1 to 1
    min_radius_range <- -0.1
    max_radius_range <- 1

    # Clamp values to [0.15, 1.0] range
    df[["plot_r"]] <- ifelse(
        df[[y]] <= 0.15,
        0.15,
        pmin(df[[y]], max_radius_range)
    )

    # Determine label positions to avoid overlaps
    df[["label_r"]] <- pmax(0.05, df[["plot_r"]] / 2)
    df[["label_theta"]] <- df[["theta_center"]]
    df[["show_label"]] <- df[[x]] >= min_label_pct

    # Identify slices that need external labels
    small_threshold <- 10

    small_indices <- which(
        df[[x]] < small_threshold |
        df[["plot_r"]] < min_radius
    )

    if (length(small_indices) > 0) {
        # Spread small labels across radial layers
        num_small <- length(small_indices)
        layer_size <- max(1, ceiling(num_small / 8))

        for (i in seq_along(small_indices)) {
            idx <- small_indices[i]
            layer <- ((i - 1) %/% layer_size) + 1

            # Position just outside the chart
            df[idx, "label_r"] <-
                max_radius_range + 0.08 + (layer - 1) * 0.12

            df[idx, "label_theta"] <- df[idx, "theta_center"]
        }
    }

    # Calculate the maximum label position
    max_label_r <- max(df[["label_r"]])

    # Extend axis range when external labels are present
    if (max_label_r > max_radius_range) {
        axis_range <- max_label_r + 0.15
    } else {
        axis_range <- max_radius_range + 0.1
    }

    # Create the polar area diagram
    fig <- plot_ly()

    # Add each sector
    for (i in seq_len(nrow(df))) {

        # Use predefined color or automatically assigned Tailwind color
        sector_color <- colors[[as.character(df[[group]][i])]]

        # Build marker configuration
        marker_config <- list(
            line = list(color = "white", width = 1),
            opacity = 1.0
        )

        if (!is.null(sector_color)) {
            marker_config[["color"]] <- sector_color
        }

        fig <- fig %>%
            add_trace(
                type = "barpolar",
                r = df[i, "plot_r"],
                theta = df[i, "theta_center"],
                width = df[i, "angle_degrees"],
                name = as.character(df[[group]][i]),
                hovertext = paste0(
                    df[[group]][i],
                    "<br>Population: ",
                    sprintf("%.3g", df[[x]][i]),
                    "%<br>Value: ",
                    sprintf("%.3g", df[[y]][i])
                ),
                hoverinfo = "text",
                marker = marker_config
            )
    }

    # Connector lines for small slices with visible external labels
    for (i in seq_len(nrow(df))) {
        if (
            df[i, "show_label"] &&
            df[i, "label_r"] > max_radius_range
        ) {
            fig <- fig %>%
                add_trace(
                    type = "scatterpolar",
                    r = c(df[i, "plot_r"], df[i, "label_r"]),
                    theta = c(
                        df[i, "theta_center"],
                        df[i, "label_theta"]
                    ),
                    mode = "lines",
                    line = list(color = "gray", width = 1),
                    showlegend = FALSE,
                    hoverinfo = "skip"
                )
        }
    }

    # Add text labels for slices above the minimum percentage
    for (i in seq_len(nrow(df))) {
        if (df[i, "show_label"]) {

            label_color <- if (
                df[i, "label_r"] > max_radius_range
            ) {
                "black"
            } else {
                "white"
            }

            # Reduce font size for small slices
            font_size <- if (
                df[[x]][i] < small_threshold ||
                df[i, "plot_r"] < min_radius
            ) {
                8
            } else {
                12
            }

            fig <- fig %>%
                add_trace(
                    type = "scatterpolar",
                    r = df[i, "label_r"],
                    theta = df[i, "label_theta"],
                    text = df[i, "label"],
                    mode = "text",
                    textfont = list(
                        size = font_size,
                        color = label_color
                    ),
                    showlegend = FALSE,
                    hoverinfo = "skip"
                )
        }
    }

    # Radial ticks
    if (show_axes) {
        tick_values <- c(0.5)
    } else {
        tick_values <- NULL
    }

    # Configure layout
    fig <- fig %>%
        layout(
            title = title,
            polar = list(
                radialaxis = list(
                    range = c(0, axis_range),
                    showticklabels = FALSE,
                    ticks = "outside",
                    tickvals = tick_values,
                    ticktext = as.character(tick_values),
                    tickmode = "array",
                    showgrid = show_axes,
                    gridcolor = "lightgray",
                    gridwidth = 1,
                    showline = FALSE
                ),
                angularaxis = list(
                    direction = "clockwise",
                    ticks = "",
                    showticklabels = FALSE,
                    showgrid = FALSE,
                    gridcolor = "lightgray",
                    gridwidth = 1,
                    showline = show_axes
                ),
                bgcolor = "white"
            ),
            hovermode = "closest",
            plot_bgcolor = "white",
            paper_bgcolor = "white"
        )

    return(fig)
}


#' Plot Ribbon chart
#'
#' Function to create connected stacked bar chart for an arbitrary number of groups
#'
plot_ribbon_chart <- function(
    df, x = "group", y = "value",
    subcat = "category",  # ie. what goes in the legend
    xlabel = NULL, ylabel = NULL, title = NULL,
    group_order = NULL,  # Specify order of groups on x-axis
    show_values = FALSE,  # Show/hide numbers inside bars
    show_axes = TRUE,     # Show/hide axes
    show_grid = TRUE,
    normalize = FALSE,
    show_xlabels = TRUE,
    bar_width = 0.4,
    colors = NULL,
    ylim = NULL,  # y-axis range, e.g. c(0, 100). NULL = auto (tozero)
    label_format = "format_number",  # "format_number", "round", or "none"
    num_digits = 3
) {
    df <- dropna(df, y)
    
    # Calculate percentages if needed
    if (normalize) {
        df <- df %>%
            group_by(!!sym(x)) %>%
            mutate(!!sym(y) := !!sym(y) / sum(!!sym(y)) * 100) %>%
            ungroup()
    }
    
    # Get unique groups in order and assign numeric positions
    if (!is.null(group_order)) {
        groups <- group_order
        missing_groups <- setdiff(unique(df[[x]]), groups)
        if (length(missing_groups) > 0) {
            warning("Some groups in data are not in group_order: ", paste(missing_groups, collapse = ", "))
        }
    } else {
        groups <- unique(df[[x]])
    }
    n_groups <- length(groups)
    group_positions <- setNames(0:(n_groups - 1), groups)
    
    # Add numeric positions and formatted labels to df
    df <- df %>%
        mutate(
            x_pos = group_positions[!!sym(x)],
            .label = if (label_format == "format_number") {
                format_number(!!sym(y), sigfigs=num_digits)
            } else if (label_format == "round") {
                round(!!sym(y), num_digits)
            } else {
                !!sym(y)
            }
        )

    # Add stacked bars for each category
    categories <- df %>%
        group_by(!!sym(subcat)) %>%
        summarise(total = sum(!!sym(y), na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(total)) %>%
        pull(!!sym(subcat))

    # Calculate cumulative positions for stacking
    df <- df %>%
        mutate(
            .category_order = match(!!sym(subcat), categories)
        ) %>%
        arrange(!!sym(x), .category_order)
    df <- df %>%
        group_by(!!sym(x)) %>%
        mutate(
            end = cumsum(!!sym(y)),
            start = lag(end, default = 0)
        ) %>%
        ungroup()
    
    # Create base plot
    fig <- plot_ly()
    
    # Handle colors
    if (is.null(colors)) {
        color_map <- setNames( rep(mpl_colors, length.out = length(categories)), categories )
    } else {

        if (is.list(colors)) {
            colors <- unlist(colors)
        }
        if (is.null(names(colors))) {
            color_map <- setNames( rep(colors, length.out = length(categories)), categories )
        } else {

            # fill missing categories with default colors
            color_map <- colors
            missing_categories <- setdiff(categories, names(color_map))
            if (length(missing_categories) > 0) {
                available_colors <- tailwind_colors[ (seq_along(missing_categories) - 1) %% length(tailwind_colors) + 1 ]
                color_map <- c( color_map, setNames(available_colors, missing_categories) )
            }
            color_map <- color_map[categories]
        }
    }

    for (i in seq_along(categories)) {
        cat <- categories[i]
        subset_data <- df[df[[subcat]] == cat, ]
        cat_color <- color_map[[cat]]
        
        fig <- fig %>%
            add_trace(
                data = subset_data,
                x = ~x_pos,
                y = as.formula(paste0("~", y)),
                type = "bar",
                name = cat,
                marker = list(color = cat_color),
                text = if (show_values) ~.label else NULL,
                textposition = if (show_values) "inside" else "none",
                width = bar_width,
                hovertemplate = paste0(cat, "<br>Value: %{y}<extra></extra>")
            )
    }
    
    # Add connecting polygons between consecutive groups
    for (i in seq_along(categories)) {
        cat <- categories[i]
        subset_data <- df[df[[subcat]] == cat, ]
        cat_color <- color_map[[cat]]

        for (j in 1:(n_groups - 1)) {
            left_group <- subset_data[subset_data[[x]] == groups[j], ]
            right_group <- subset_data[subset_data[[x]] == groups[j + 1], ]
            
            if (nrow(left_group) > 0 && nrow(right_group) > 0) {
                left_pos <- left_group[['x_pos']]
                right_pos <- right_group[['x_pos']]
                
                left_edge <- left_pos + bar_width / 2
                right_edge <- right_pos - bar_width / 2
                
                x_coords <- c(left_edge, right_edge, right_edge, left_edge, left_edge)
                y_coords <- c(
                    left_group[['start']], 
                    right_group[['start']], 
                    right_group[['end']], 
                    left_group[['end']],
                    left_group[['start']]
                )
                
                fig <- fig %>%
                    add_trace(
                        x = x_coords,
                        y = y_coords,
                        type = "scatter",
                        mode = "lines",
                        fill = "toself",
                        fillcolor = paste0(cat_color, "40"),
                        line = list(color = cat_color, width = 1),
                        showlegend = FALSE,
                        hoverinfo = "skip"
                    )
            }
        }
    }

    # Build y-axis range config
    yaxis_range <- if (!is.null(ylim)) {
        if (length(ylim) != 2) stop("`ylim` must be a length-2 vector, e.g. c(0, 100)")
        list(range = ylim, rangemode = "normal")
    } else {
        list(rangemode = "tozero")
    }

    fig <- fig %>%
        layout(
            title = title,
            barmode = "stack",
            xaxis = list(
                title = if (!is.null(xlabel)) xlabel else "",
                tickmode = "array",
                tickvals = 0:(n_groups - 1),
                ticktext = groups,
                showticklabels = show_xlabels,
                showline = show_axes
            ),
            yaxis = c(
                list(
                    title = if (!is.null(ylabel)) ylabel else "",
                    showticklabels = show_axes,
                    showline = show_axes,
                    showgrid = show_grid,
                    zeroline = FALSE
                ),
                yaxis_range
            ),
            hovermode = "closest",
            plot_bgcolor = "white",
            paper_bgcolor = "white"
        )
    
    return(fig)
}

