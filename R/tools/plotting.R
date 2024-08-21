import::here(ggplot2,
    'ggplot', 'aes', 'theme', 'labs',
    'geom_jitter', 'element_text')
import::here(plotly, 'plot_ly', 'add_trace', 'layout')
import::here(magrittr, '%>%')
# import::here(rlang, 'sym')

## Functions
## plot_dots
## plot_violin


#' Plot Dots
#'
#' @description  Thin wrapper around ggplot geom_jitter
#'
plot_dots <- function(
    df,
    x='cell_type',
    y='pct_cells',
    color='treatment_group',  # must be string
    xlabel='Cell Type',
    ylabel='Percent of Live Cells',
    title=NULL
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
    df, x, y, group_by, size=NULL,
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
            df <- df[order(-df[['pct_cells']]), ]
        } else {
            df <- df[order(df[['pct_cells']]), ]
        }
    }

    fig <- plot_ly(type = 'violin')
    for (group in unique(df[[group_by]])) {
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
        showlegend = TRUE,
        hovermode = 'closest'

    )

    return(fig)
}
