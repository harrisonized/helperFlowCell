# import::here(rlang, 'sym')
import::here(magrittr, '%>%')
import::here(dplyr, 'group_by', 'summarize')
import::here(tidyr, 'pivot_wider')
import::here(ggplot2,
    'ggplot', 'aes', 'theme', 'labs',
    'geom_boxplot', 'geom_jitter', 'element_text',
    'stat_summary', 'scale_fill_brewer', 'scale_x_discrete', 'scale_y_continuous',
    'guide_axis', 'expansion', 'ggtitle')
import::here(superb, 'showSignificance')
import::here(plotly, 'plot_ly', 'add_trace', 'layout', 'save_image')
import::here(htmlwidgets, 'saveWidget')  # brew install pandoc

import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'flatten_matrix', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'math.R'),
    'unpaired_t_test', .character_only=TRUE)

## Functions
## save_fig
## plot_dots
## plot_scatter
## plot_violin
## compute_nlevels
## generate_base_level
## plot_violin_with_significance


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
#' Used for plot_violin_with_significance
#' Calculate the number of non-overlapping bars needed to span every significance level
#' In practice, this is not actually necessary
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
#' Used for plot_violin_with_significance
#' Generate a sequence used as a floor for each particular level
#' This serves as the correction factor for the bar height of each p_value bar
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


#' Plot violin with significance
#' 
#' @description Produces a violin plot with error bars in between
#' You will need to pre-calculate the pvalues and label them as 1v2, 1v3, etc.
#' Do not include any additional columns in the pvals object
#' Will update this function to calculate pvals within the function
#' 
#' @return Returns a ggplot object
#'
plot_violin_with_significance <- function(
    df, pvals,
    x='group_name', y='pct_cells', title=NULL,
    xaxis_angle=60
    # test='t_test'  # unused for now
) {

    # compute bar positions
    # TODO: cache this for quick lookup
    n_groups <- length(unique(df[[x]]))

    if (n_groups > 1) {
        bracket_params <- setNames(as.data.frame(t(combn(1:n_groups, 2))), c('left', 'right'))
        bracket_params[['dist']] <- bracket_params[['right']] - bracket_params[['left']]
        bracket_params <- bracket_params[ order(bracket_params[['dist']]), ]
        bracket_params[['base_level']] <- generate_base_level(n_groups)
        bracket_params[['level']] <- bracket_params[['base_level']] +  # base level
            (bracket_params[['left']]-1) %% bracket_params[['dist']]  # alternating correction factor

        # compute starting height
        h_low <- max(aggregate(df[[y]],
            list(df[[x]]), FUN=function(x) mean(x)+sd(x)
        )[['x']], na.rm=TRUE) * 1.1
        if (h_low < max(df[[y]])) {
            h_low <- max(df[[y]]) * 1.1
        }
        space <- h_low / 10
    }

    # if fewer groups than expected, recompute pvals to get the correct significance labels
    nlevels <- compute_nlevels( length(unique(df[[x]])) )
    if ( (nlevels > 0) & (length(pvals) > nlevels) ) {
        group_names <- sort(unique(df[[x]]))
        n_groups <- length(group_names)
        n_combos <- choose(n_groups, 2)
        id_combos <- flatten_matrix(combn(1:n_groups, 2))  # generate pairs of indexes
        pval_cols <- sapply(id_combos, function(x) paste0('pval_', x[[1]], 'v', x[[2]]))  # colnames for all pairs

        # Collect values into list columns
        pvals <- pivot_wider(
            df[, c(x, y)],
            names_from = c(x),
            values_from = c(y),
            values_fn = list,  # suppress warning
            names_glue = "{.name}"
        )
        pvals <- pvals[, group_names]  # sort cols

        # calculate unpaired t test for all pairs of groups
        for (idx in 1:n_combos) {
            idx1 <- id_combos[[idx]][1]  # 1st col idx
            idx2 <- id_combos[[idx]][2]  # 2nd col idx
            pvals[ pval_cols[idx] ] <- mapply(
                function(x, y) unpaired_t_test(x, y),
                pvals[[ group_names[idx1] ]],  # 1st col
                pvals[[ group_names[idx2] ]]   # 2nd col
            )
        }
    }
    
    fig <- ggplot(df, aes(x=.data[[x]], y=.data[[y]], fill=.data[[x]])) + 
        geom_boxplot(alpha=0.7, aes(middle=.data[[y]])) +
        stat_summary(fun=mean, geom="crossbar", width=0.75, linewidth=0.25, linetype = "dashed") +
        geom_jitter() +
        scale_fill_brewer(palette="Dark2") +
        scale_x_discrete(guide=guide_axis(angle=xaxis_angle)) +
        scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
        ggtitle(title)

    # significance brackets
    if (n_groups > 1) {
        for (row in 1:nrow(bracket_params)) {

            left <- bracket_params[row, "left"]
            right <- bracket_params[row, "right"]
            level <- bracket_params[row, "level"]
            col <- paste0('pval_', left, 'v', right)  # colname

            fig <- fig +
                showSignificance( c(left+0.1, right-0.1), h_low+(level-1)*space, -0.001*h_low,
                    toString(round(pvals[[col]], 4)) )
        }
    }

    return(fig)
}
