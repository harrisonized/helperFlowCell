import::here(ggplot2,
    'ggplot', 'aes', 'theme', 'labs',
    'geom_jitter', 'element_text')
# import::here(rlang, 'sym')

## Functions
## plot_dots


#' Plot Dots
#'
#' @description  Wrapper around geom_jitter
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
