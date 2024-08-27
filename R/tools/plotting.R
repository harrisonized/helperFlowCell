# import::here(rlang, 'sym')
import::here(magrittr, '%>%')
import::here(dplyr, 'group_by', 'summarize')
import::here(ggplot2,
    'ggplot', 'aes', 'theme', 'labs',
    'geom_jitter', 'element_text')
import::here(plotly, 'plot_ly', 'add_trace', 'layout')

## Functions
## plot_dots
## plot_violin


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


#' Plot Violin
#'
#' @description Thin wrapper arounbd Plotly violin plot
#' 
plot_violin <- function(
    df, x, y, group_by=NULL, size=NULL,
    xlabel=NULL, ylabel=NULL, title=NULL,
    ymin=NULL, ymax=NULL,
    hover_data=c(),
    yaxis_type='linear',
    color_discrete_map=NULL,
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
                summarize(mean = median(.data[[y]])) %>%
                arrange(-.data[['mean']]) %>%
                select(.data[[x]])
        } else {
            x_axis_order <- df %>%
                group_by(.data[[x]]) %>%
                summarize(mean = median(.data[[y]])) %>%
                arrange(.data[['mean']]) %>%
                select(.data[[x]])
        }
    }

    fig <- plot_ly(type = 'violin')

    if (is.null(group_by)) {
        if (sort) {
            df <- df %>% arrange(match(.data[[x]], unlist(x_axis_order)))
        }

        fig <- fig %>%
            add_trace(
                x = df[[x]],
                y = df[[y]],
                box = list(visible = TRUE),
                meanline = list(visible = TRUE),
                points = 'all',
                jitter = 0.2,
                pointpos = -1,
                marker = list(size = 5)
            )
    } else {
        for (group in unique(df[[group_by]])) {

            if (sort) {
                df <- df %>% arrange(
                    .data[[group_by]],
                    match(.data[[x]], unlist(x_axis_order))
                )
            }

            fig <- fig %>%
                add_trace(
                    x = unlist(df[df[[group_by]] == group, x]),
                    y = unlist(df[df[[group_by]] == group, y]),
                    legendgroup = group,
                    scalegroup = group,
                    name = group,
                    box = list(visible = TRUE),
                    meanline = list(visible = TRUE),
                    points = 'all',
                    jitter = 0.2,
                    pointpos = -1,
                    marker = list(size = 5)
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
