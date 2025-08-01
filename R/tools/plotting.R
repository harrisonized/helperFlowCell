# import::here(rlang, 'sym')
import::here(magrittr, '%>%')
import::here(dplyr, 'group_by', 'summarize')
import::here(tidyr, 'pivot_wider', 'unnest')
import::here(ggplot2,
    'ggplot', 'aes', 'theme', 'labs', 'theme_minimal',
    'geom_boxplot', 'geom_jitter', 'geom_col', 'geom_line', 'element_text',
    'stat_summary', 'scale_fill_brewer', 'scale_fill_manual', 'scale_color_manual',
    'scale_x_continuous', 'scale_x_discrete', 'scale_y_continuous',
    'guide_axis', 'expansion', 'ggtitle')
import::here(flowCore, 'logicleTransform')
import::here(RColorBrewer, 'brewer.pal')
import::here(ggprism, 'theme_prism')
import::here(superb, 'showSignificance')
import::here(plotly, 'plot_ly', 'add_trace', 'layout', 'save_image')
import::here(htmlwidgets, 'saveWidget')  # brew install pandoc

import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'math.R'),
    'unpaired_t_test', 'fishers_lsd', 'tukey_multiple_comparisons',
    'bonferroni_multiple_comparisons', 'generate_gaussian_data',
    .character_only=TRUE)

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

    suppressWarnings(save_image(fig,
        file=file.path(file.path(dirpath, 'png'), paste0(filename, '.png')), 
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
#' @description Thin wrapper arounbd Plotly scatter plot
#' 
plot_scatter <- function(
    df, x, y, group_by=NULL, size=NULL,
    xlabel=NULL, ylabel=NULL, title=NULL,
    xmin=NULL, xmax=NULL,
    ymin=NULL, ymax=NULL,
    yaxis_type='linear',
    color='#1f77b4',
    color_discrete_map=NULL,
    hover_data=c()
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

    # color
    if (!is.null(color_discrete_map)) {
        color_discrete_map <- list2env(as.list(color_discrete_map))
    }

    fig <- plot_ly(type='scatter', mode='markers')

    if (is.null(group_by)) {

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
    descending=TRUE
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
        violinmode = 'group',
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
                        aggregate(df_subset[[y]], list(df_subset[[x]]),
                        FUN=function(x) mean(x)+sd(x)
                    )[['x']], na.rm=TRUE) * 1.1
                }, warning = function(w) {
                    if ( grepl("no non-missing arguments", w) ) {
                        stop()
                    }
                })
            },
            error = function(condition) {
                h_low <<- max(df_subset[[y]], na.rm=TRUE) * 1.1
            })
            if (h_low < max(df[[y]])) {
                h_low <- max(df[[y]]) * 1.1
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


#' Plot Biex Histograms
#' 
#' Mimics Flowjo's histogram function
#' Use [generate_gaussian_data()] to generate the input
#' 
plot_modal_histograms <- function(df,
    group='group', value='value',
    xlabel=NULL, ylabel=NULL, title=NULL,
    colors=c(),  # c("#7F7F7F", "#ACF53B")
    show_bins=FALSE,
    w=2.25, t=10000, m=4.5, a=0  # logicalTransform parameters
) {

    if ( w > m/2 ) {
        w <- m/2
    }
    # Logicle transform
    logicle <- logicleTransform(w=w, t=t, m=m, a=a)  # default options
    df[['logicle_value']] <- logicle(df[[value]])

    # Histogram binning
    bins <- 100
    breaks <- pretty(range(df[['logicle_value']]), n = bins)

    # Compute normalized histogram
    hist_data <- df %>%
        group_by(group) %>%
        summarise(
            counts = list(hist(logicle_value, breaks = breaks, plot = FALSE)$counts),
            mids = list(hist(logicle_value, breaks = breaks, plot = FALSE)$mids)
        ) %>%
        unnest(cols = c(counts, mids)) %>%
        group_by(group) %>%
        mutate(norm_counts = counts / max(counts))
    
    # Smoothed and clamped interpolation
    interp_data <- hist_data %>%
        group_by(group) %>%
        reframe({
            fit <- smooth.spline(mids, norm_counts, spar = 0.7)
            x_vals <- seq(min(mids), max(mids), length.out = 2000)
            y_vals <- predict(fit, x = x_vals)$y
            tibble(x = x_vals, y = pmax(y_vals, 0), group = first(group))
        })

    # Logicle axis ticks
    # This might not be accurate, need to figure this part out
    raw_breaks <- c(0, 10, 50, 100, 200, 500, 1000)
    transformed_breaks <- logicle(raw_breaks)

    # repeat last color value if not enough colors
    n_groups <- n_distinct(df[[group]])
    n_colors <- length(colors)
    if ((n_colors > 0) & (n_colors < n_groups)) {
        colors <- c(colors, rep(colors[length(colors)], n_groups-n_colors))
    }

    # plot
    fig <- ggplot() +
        (if (show_bins) {
            geom_col(data = hist_data, aes(x = mids, y = norm_counts, fill = group),
                     position = "identity", alpha = 0.4, width = diff(breaks)[1])  # area
        } else {
            geom_area(data = interp_data, aes(x = x, y = y, fill = group), alpha = 0.4, position = "identity")
        }) +
        geom_line(data = interp_data, aes(x = x, y = y, color = group), size = 1.2) +  # line
        (if (length(colors) == n_groups) { scale_fill_manual(values = colors) } else NULL) +  # area
        (if (length(colors) == n_groups) { scale_color_manual(values = colors) } else NULL) +  # line
        scale_x_continuous(breaks = transformed_breaks, labels = raw_breaks) +
        labs(title = title, x = xlabel, y = ylabel) +
        theme_minimal(base_size = 14)

    return(fig)
}


# temporary test code
if (FALSE) {
    df1 <- generate_gaussian_data(mean=200, group_name="Mean 200")
    df2 <- generate_gaussian_data(mean=500, group_name="Mean 500")
    df3 <- generate_gaussian_data(mean=700, group_name="Mean 700")
    df <- rbind(df1, df2, df3)

    plot_modal_histograms(df,
        xlabel ="Comp-GFP-A :: mNeonGreen",
        ylabel = "Normalized to Mode",
        title = "Simulated Data",
        colors = c("#7F7F7F", "#ACF53B")
    )
}
