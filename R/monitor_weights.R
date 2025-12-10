## Monitor survival
## 
## Facilitates survival monitoring by generating the following:
## 1. Kaplan-Meyer Survival curve with confidence intervals
## 2. Percent weight over time and raw weight over time
## 3. Spleen weight comparison if spleen_weight is provided in the weights table

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
library('optparse')
suppressPackageStartupMessages(library('logr'))
suppressPackageStartupMessages(library(dplyr))

import::from(magrittr, '%>%')
import::from(tibble, 'tibble')
import::from(tidyr, 'pivot_longer')
import::from(broom, 'tidy')
import::from(survival, 'survfit', 'Surv')
import::from(ggplot2,
    'ggplot', 'aes', 'geom_step', 'geom_ribbon',
    'labs', 'xlim', 'ylim', 'ggsave')
import::from(ggprism, 'theme_prism')

import::from(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'create_survival_table', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_scatter', 'plot_multiple_comparisons', 'save_fig', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'read_excel_or_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'append_dataframe', 'fillna', 'rename_columns', .character_only=TRUE)

# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-file"), default='data/weights.xlsx',
                metavar='data/weights.xlsx', type="character",
                help="specify input file"),

    make_option(c("-o", "--output-dir"), default="figures/survival",
                metavar="figures/survival", type="character",
                help="set the output directory for the data"),

    make_option(c("-g", "--group-by"), default='treatment,genotype',
                metavar='treatment,genotype', type="character",
                help="enter a column or comma-separated list of columns, no spaces"),

    make_option(c("-w", "--week"), default=FALSE, action="store_true",
                metavar='FALSE', type="logical",
                help="use this to select week for the plots"),

    make_option(c("-p", "--optional-plots"), default=FALSE, action="store_true",
                metavar='FALSE', type="logical",
                help="use this to select optional plots"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]

# args
metadata_cols <- unlist(strsplit(opt[['group-by']], ','))

# Start Log
start_time = Sys.time()
log <- log_open(paste0("monitor_weights-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read file

log_print(paste(Sys.time(), 'Reading data...'))

input_file <- file.path(wd, opt[['input-file']])
weight_tbl <- read_excel_or_csv(input_file)

ext=tools::file_ext(input_file)
if (ext == 'xlsx') {
    date_pattern <- "^\\d+$"
} else if (ext == 'csv') {
    date_pattern <- "^\\d{1,2}/\\d{1,2}/\\d{2,4}$"
}
is_date <- grepl(date_pattern, colnames(weight_tbl))
id_cols <- colnames(weight_tbl)[!is_date]
date_cols <- colnames(weight_tbl)[is_date]
weight_tbl <- weight_tbl[(!is.na(weight_tbl[['mouse_id']])), ]

# Derive group col
metadata_cols <- intersect(metadata_cols, colnames(weight_tbl))
if (length(metadata_cols) > 1) {
    weight_tbl[['group']] <- apply( weight_tbl[ , metadata_cols ] , 1 , paste , collapse = ", " )
} else {
    weight_tbl[['group']] <- weight_tbl[[metadata_cols]]
}

# pivot
df = pivot_longer(
    data=weight_tbl,
    cols=all_of(date_cols),
    names_to='date',
    values_to='weight',
    values_drop_na = FALSE  # important
)

# fix dates
if (ext == 'xlsx') {
    df[['date']] <- as.Date(as.numeric(df[['date']]), origin = "1899-12-30")
}


# ----------------------------------------------------------------------
# Plot survival

log_print(paste(Sys.time(), 'Plotting survival...'))

# extract start and end info
start_info <- df %>%
  group_by(mouse_id) %>%
  slice_min(date, with_ties = FALSE) %>%
  ungroup() %>%
  select(any_of(c('mouse_id', 'spleen_weight', 'group', 'date', 'weight')))

end_info <- df %>%
    filter(!is.na(weight)) %>%
    group_by(mouse_id) %>%
    slice_max(order_by=date, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(all_of(c('mouse_id', 'group', 'date', 'weight')))

survival <- create_survival_table(start_info, end_info, end_date=max(df[['date']]))


# Get confidence intervals
if (opt[['week']]) {
    fit <- survfit(Surv(week, is_dead) ~ group, data = survival)
} else {
    fit <- survfit(Surv(day, is_dead) ~ group, data = survival)
}
surv_data <- tidy(fit)
surv_data[['strata']] <- gsub('group=', '', surv_data[['strata']])


# Plot
fig <- ggplot(surv_data, aes(x = time, y = estimate, color = strata)) +
    geom_step(linewidth = 1) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = strata), alpha = 0.2, color = NA) +
    labs(x = paste("Time", if (opt[['week']]) {'(weeks)'} else {'(days)'}),
         y = "Survival Probability",
         color = "Genotype",
         fill = "Genotype",
         title = "Survival Curve") +
    xlim(0, NA) +
    ylim(0, 1) +
    theme_prism()

# save
if (!troubleshooting) {
    dirpath <- file.path(wd, opt[['output-dir']])
    if (!dir.exists(dirpath)) { dir.create(dirpath, recursive=TRUE) }
    ggsave(
        file.path(dirpath, 'survival-curve.png'),  
        plot=fig,
        height=4000, width=6000, dpi=500, units='px',
        scaling=2
    )
}


# ----------------------------------------------------------------------
# Plot Weights

log_print(paste(Sys.time(), 'Plotting weights...'))


# output_dir for this section
if (!troubleshooting) {
    dirpath <- file.path(wd, opt[['output-dir']], 'mouse-weights')
    if (!dir.exists(dirpath)) { dir.create(dirpath, recursive=TRUE) }
}


# Preprocessing
df <- merge(df,  start_info[, c('mouse_id', 'date', 'weight')],
    by="mouse_id", all.x=TRUE, all.y=TRUE, suffixes=c('', '_start'))

df[['day']] <- as.numeric(difftime(
    as.Date(df[['date']], format = '%m/%d/%y'),
    as.Date(df[['date_start']], format = '%m/%d/%y'),
    units = paste0('day', 's')
))
df[['week']] <- df[['day']]/7
df[['pct_weight']] <- df[['weight']] / df[['weight_start']]
df <- data.frame(df)
df <- append_dataframe(df,
    cbind(weight_tbl[, c('mouse_id', 'group')],
        data.frame(
            day=max(df[['day']]+1),
            week=max(df[['week']]+1/7)
        )
    )
)
df <- df[order(df[['mouse_id']], desc(df[['day']])), ]

# set colors
# TODO: un-hardcode this
color_map <- data.frame(unique(df[c('mouse_id', 'group')]))
color_map[['color']] <- lapply(
    color_map[['group']],
    function(x) if ((x=='WT') | (x=='H')) {
        '#FF7F0E'  # orange
    } else if (x=='M') {
        '#377EB8'  # blue
    } else if (x=='PBS') {
        '#9467bd'  # purple
    } else {
        'rgba(23,190,207,1)'  # sea green
    }
)
color_discrete_map <- setNames(color_map[['color']], color_map[['group']])


# Percent weight
fig <- plot_scatter(
    df,
    x=if (opt[['week']]) {'week'} else {'day'},
    y='pct_weight', group_by='group',
    ymin=0.75,
    xlabel=paste("Time", if (opt[['week']]) {'(weeks)'} else {'(days)'}),
    ylabel='Percent Weight',
    title='Percent Weight Over Time',
    mode='lines+markers',
    color_discrete_map=color_discrete_map,
    hover_data=c('mouse_id', 'sex', 'genotype', 'dob', 'day', 'week')
)

# save
if (!troubleshooting) {
    save_fig(
        fig=fig,
        height=400, width=750,
        dirpath=dirpath,
        filename='pct_weight',
        save_html=TRUE
    )
}

# Raw weight
fig <- plot_scatter(
    df,
    x=if (opt[['week']]) {'week'} else {'day'},
    y='weight', group_by='group',
    xlabel=paste("Time", if (opt[['week']]) {'(weeks)'} else {'(days)'}),
    ylabel='Mouse Weight (g)',
    title='Raw Weight Over Time',
    mode='lines+markers',
    color_discrete_map=color_discrete_map,
    hover_data=c('sex', 'genotype', 'dob', 'day', 'week')
)

# save
if (!troubleshooting) {
    save_fig(
        fig=fig,
        height=400, width=600,
        dirpath=dirpath,
        filename='raw_weight',
        save_html=TRUE
    )
}


# Weight loss at time of death
survival[['weight_ratio']] <- survival[['weight_end']] / survival[['weight_start']] 

fig <- plot_multiple_comparisons(
    survival,
    x='group', y='weight_ratio',
    ymin=0.70,
    ylabel="[Final weight (g)] / [Starting weight (g)]",
    title='Weight Loss at Time of Death',
    show_numbers=TRUE,
    test='fishers_lsd',
    custom_group_order=intersect(c('WT', 'KO'), unique(weight_tbl[['group']]))
)

# save
if (!troubleshooting) {
    withCallingHandlers({
        ggsave(
            file.path(wd, opt[['output-dir']], 'weight-loss.png'),
            plot=fig,
            height=2000, width=2400, dpi=300, units='px',
            scaling=1
        )
    }, warning = function(w) {
        if ( any(grepl("containing non-finite values", w),
                 grepl("outside the scale range", w),
                 grepl("fewer than two data points", w),
                 grepl("argument is not numeric or logical: returning NA", w)) ) {
            invokeRestart("muffleWarning")
        }
    })
}


# ----------------------------------------------------------------------
# Plot spleen weight


if ('spleen_weight' %in% colnames(weight_tbl)) {
    
    # output_dir for this section
    if (!troubleshooting) {
        dirpath <- file.path(wd, opt[['output-dir']], 'spleen-weight')
        if (!dir.exists(dirpath)) { dir.create(dirpath, recursive=TRUE) }
    }

    # ----------------------------------------------------------------------
    # Spleen Weight by Group

    fig <- plot_multiple_comparisons(
        survival,
        x='group', y='spleen_weight',
        ylabel="Spleen Weight (mg)",
        title='Spleen Weight',
        show_numbers=FALSE,
        test='fishers_lsd',
        custom_group_order=intersect(c('WT', 'KO'), unique(weight_tbl[['group']]))
    )

    # save
    if (!troubleshooting) {
        withCallingHandlers({
            ggsave(
                file.path(wd, opt[['output-dir']], 'spleen-weight', 'spleen_weight.svg'),
                plot=fig,
                height=5000, width=5000, dpi=500, units='px',
                scaling=1
            )
        }, warning = function(w) {
            if ( any(grepl("containing non-finite values", w),
                     grepl("outside the scale range", w),
                     grepl("fewer than two data points", w),
                     grepl("argument is not numeric or logical: returning NA", w)) ) {
                invokeRestart("muffleWarning")
            }
        })
    }


    # ----------------------------------------------------------------------
    # Spleen Weight Ratio vs. Treatment Duration

    fig <- plot_scatter(
        survival,
        x=if (opt[['week']]) {'week'} else {'day'},
        y='spleen_weight',
        group_by='group',
        ymin=0,
        xlabel=paste("Time of Death", if (opt[['week']]) {'(weeks)'} else {'(days)'}),
        ylabel='Spleen Weight (mg)',
        title='Spleen Weight vs. Treatment Duration',
        mode='markers'
    )

    # save
    if (!troubleshooting) {
        save_fig(
            fig=fig,
            height=350, width=500,
            dirpath=dirpath,
            filename='spleen_weight_vs_treatment_duration',
            save_html=TRUE
        )
    }


    # ----------------------------------------------------------------------
    # Spleen Weight vs. Mouse Starting Weight

    fig <- plot_scatter(
        survival,
        x='weight_start', y='spleen_weight', group_by='group',
        ymin=0,
        xlabel='Mouse Starting Weight (g)', ylabel='Spleen Weight (mg)',
        title='Spleen Weight vs. Mouse Starting Weight',
        mode='markers',
        # color_discrete_map=color_discrete_map,
        hover_data=c('sex', 'group')
    )

    # save
    if (!troubleshooting) {
        save_fig(
            fig=fig,
            height=350, width=500,
            dirpath=dirpath,
            filename='spleen_vs_mouse_weight',
            save_html=TRUE
        )
    }


    # optional plots
    if (opt[['optional-plots']]) {

        # ----------------------------------------------------------------------
        # Spleen Weight Ratio vs. Treatment Duration

        survival[['spleen_weight_ratio_start']] <- survival[['spleen_weight']] / survival[['weight_start']]

        fig <- plot_scatter(
            survival,
            x=if (opt[['week']]) {'week'} else {'day'},
            y='spleen_weight_ratio_start',
            group_by='group',
            ymin=0,
            xlabel=paste("Time of Death", if (opt[['week']]) {'(weeks)'} else {'(days)'}),
            ylabel='[Spleen Weight / Mouse Starting Weight] (mg/g)',
            title='Spleen Weight Ratio vs. Treatment Duration',
            mode='markers'
        )


        # save
        if (!troubleshooting) {
            save_fig(
                fig=fig,
                height=400, width=600,
                dirpath=dirpath,
                filename='spleen_weight_ratio_start_vs_treatment_duration',
                save_html=TRUE
            )
        }

        # ----------------------------------------------------------------------
        # Spleen Weight Ratio vs. Treatment Duration

        survival[['spleen_weight_ratio_end']] <- survival[['spleen_weight']] / survival[['weight_end']]

        fig <- plot_scatter(
            survival,
            x=if (opt[['week']]) {'week'} else {'day'},
            y='spleen_weight_ratio_end',
            group_by='group',
            ymin=0,
            xlabel=paste("Time of Death", if (opt[['week']]) {'(weeks)'} else {'(days)'}),
            ylabel='[Spleen Weight / Mouse Final Weight] (mg/g)',
            title='Spleen Weight Ratio vs. Treatment Duration',
            mode='markers'
        )

        # save
        if (!troubleshooting) {
            save_fig(
                fig=fig,
                height=400, width=600,
                dirpath=dirpath,
                filename='spleen_weight_ratio_final_vs_treatment_duration',
                save_html=TRUE
            )
        }
    }


}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
