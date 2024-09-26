# import::here(rlang, 'sym')
import::here(magrittr, '%>%')
import::here(dplyr, 'group_by', 'summarize')
import::here(ggplot2,
    'ggplot', 'aes', 'theme', 'labs',
    'geom_jitter', 'element_text')
import::here(plotly, 'plot_ly', 'add_trace', 'layout', 'save_image')
import::here(htmlwidgets, 'saveWidget')  # brew install pandoc

## Functions
## save_fig
## plot_dots
## plot_scatter
## plot_violin


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
    if (opt[['save-html']]) {
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
    } else {
        for (group in unique(df[[group_by]])) {

            if (sort) {
                df <- df %>% arrange(
                    .data[[group_by]],
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
    }

    fig <- fig %>% layout(
        title = list(
            text = title,
            x = 0.5
        ),
        xaxis = list(title_text = xlabel),
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
