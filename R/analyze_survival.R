## Analyze survival
## Plot Kaplan-Meyer survival curve and weight over time

wd = dirname(this.path::here())  # wd = '~/github/R/helperFlowCell'
library('optparse')
suppressPackageStartupMessages(library('logr'))
suppressPackageStartupMessages(library(dplyr))

import::from(broom, 'tidy')
import::from(magrittr, '%>%')
import::from(tibble, 'tibble')
import::from(tidyr, 'pivot_longer')
import::from(survival, 'survfit', 'Surv')
import::from(ggplot2,
    'ggplot', 'aes', 'geom_step', 'geom_ribbon',
    'labs', 'xlim', 'ylim', 'ggsave'
)
import::from(ggprism, 'theme_prism')

import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_scatter', 'plot_multiple_comparisons', 'save_fig', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'read_excel_or_csv', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-file"), default='data/weights.xlsx',
                metavar='data/weights.xlsx', type="character",
                help="specify input file"),

    make_option(c("-o", "--output-dir"), default="figures/time",
                metavar="figures/time", type="character",
                help="set the output directory for the data"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("analyze_survival-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read file


log_print(paste(Sys.time(), 'Reading data...'))

input_file <- file.path(wd, opt[['input-file']])
weight_tbl <- read_excel_or_csv(input_file)

ext=tools::file_ext(input_file)
if (ext == 'xlsx') {
    pattern <- "^\\d+$"
} else if (ext == 'csv') {
    pattern <- "^\\d{1,2}/\\d{1,2}/\\d{2,4}$"
}
is_date <- grepl(pattern, colnames(weight_tbl))
id_cols <- colnames(weight_tbl)[!is_date]
date_cols <- colnames(weight_tbl)[is_date]
weight_tbl <- weight_tbl[(!is.na(weight_tbl[['mouse_id']])), ]


# ----------------------------------------------------------------------
# Preprocessing

log_print(paste(Sys.time(), 'Reshaping data...'))

df = pivot_longer(
    data=weight_tbl,
    cols=all_of(date_cols),
    names_to='date',
    values_to='weight',
    values_drop_na = FALSE  # important
)

if (ext == 'xlsx') {
    df[['date']] <- as.Date(as.numeric(df[['date']]), origin = "1899-12-30")
}

start_info <- df %>%
  group_by(mouse_id) %>%
  slice_min(date, with_ties = FALSE) %>%
  ungroup() %>%
  select(mouse_id, date, weight)

# start_info <- rename_columns(start_info, c('date'='date_start', 'weight'='weight_start'))

df <- merge(df, start_info,
    by="mouse_id", all.x=TRUE, all.y=TRUE, suffixes=c('', '_start'))

df[['day']] <- as.numeric(difftime(
    as.Date(df[['date']], format = "%m/%d/%y"),
    as.Date(df[['date_start']], format = "%m/%d/%y"),
    units = "days"
))

df[['week']] <- as.numeric(difftime(
    as.Date(df[['date']], format = "%m/%d/%y"),
    as.Date(df[['date_start']], format = "%m/%d/%y"),
    units = "weeks"
))


df[['pct_weight']] <- df[['weight']] / df[['weight_start']]
df[['is_dead']] <- as.numeric(is.na(df[['pct_weight']]))


# label died column
df <- df %>%
    group_by(mouse_id) %>%
    mutate(
        died = if (any(is.na(weight))) {
            as.integer(row_number() == min(row_number()[is.na(weight)]))
        } else { 0 }
    ) %>%
    ungroup()

df[['died']] <- lead(df[['died']])

df <- data.frame(df)


# ----------------------------------------------------------------------
# Plot survival


log_print(paste(Sys.time(), 'Plotting survival...'))


# drop
df_subs <- df %>%
  group_by(mouse_id) %>%
  mutate(first_one_index = match(1, died)) %>%
  filter(row_number() <= first_one_index) %>%
  select(-first_one_index) %>%
  ungroup()


events <- df_subs %>%
  group_by(mouse_id) %>%
  slice_max(day, with_ties = FALSE) %>%
  ungroup() %>%
  select(mouse_id, genotype, week, died)

# Add dummy row with day=0, died=0 for each genotype group
events <- events %>%
  group_by(genotype) %>%
  group_modify(~ bind_rows(tibble(week = 0, died = 0), .x)) %>%
  ungroup()

# Find first death day per genotype group
first_deaths <- events[(events[['died']]==1), ] %>%
  group_by(genotype) %>%
  summarise(first_death = min(week)) %>%
  ungroup()

# Create dummy rows for day = 0 and day just before first death
dummy_rows <- first_deaths %>%
  mutate(
    week_before = first_death - 0.001
  ) %>%
  select(genotype, week_before) %>%
  pivot_longer(cols = c(week_before), names_to = "time_type", values_to = "week") %>%
  # Add day=0 rows for each genotype as well
  bind_rows(
    tibble(
      genotype = unique(first_deaths$genotype),
      day = 0
    )
  ) %>%
  distinct() %>%  # remove duplicates in case day_before = 0
  arrange(genotype, week) %>%
  mutate(
    mouse_id = NA_integer_,
    died = 0
  ) %>%
  select(mouse_id, genotype, week, died)

# Bind dummy rows to events and arrange
events_extended <- events %>%
  bind_rows(dummy_rows) %>%
  arrange(genotype, week)

# Fit survival model
fit <- survfit(Surv(week, died) ~ genotype, data = events_extended)
surv_data <- tidy(fit)
surv_data[['strata']] <- gsub('genotype=', '', surv_data[['strata']])

fig <- ggplot(surv_data, aes(x = time, y = estimate, color = strata)) +
    geom_step(linewidth = 1) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = strata), alpha = 0.2, color = NA) +
    labs(
        x = "Time (weeks)",
        y = "Survival probability",
        color = "Genotype",
        fill = "Genotype",
        title = "Survival Curve"
    ) +
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
        height=2000, width=3000, dpi=250, units='px',
        scaling=2
    )
}



# ----------------------------------------------------------------------
# Plot weights

df <- data.frame(df)
df[['group_name']] <- paste(df[['mouse_id']], df[['genotype']], sep=', ')


# set colors
color_map <- data.frame(unique(df[c('mouse_id', 'genotype', 'group_name')]))
color_map[['color']] <- lapply(
    color_map[['genotype']],
    function(x) if (x=='WT') {
        '#FF7F0E'  # orange
    } else {
        'rgba(23,190,207,1)'  # sea green
    }
)
color_discrete_map <- setNames(color_map[['color']], color_map[['group_name']])

log_print(paste(Sys.time(), 'Plotting weights...'))

fig <- plot_scatter(
    df,
    x='week', y='pct_weight', group_by='group_name',
    ymin=0,
    xlabel='Time (week)', ylabel='Percent Weight', title='Weight Curves',
    mode='lines+markers',
    color_discrete_map=color_discrete_map,
    hover_data=c('sex', 'genotype', 'dob', 'day', 'week')
)

# save
if (!troubleshooting) {
    dirpath <- file.path(wd, opt[['output-dir']], 'weight-curves')
    if (!dir.exists(dirpath)) { dir.create(dirpath, recursive=TRUE) }
    save_fig(
        fig=fig,
        height=350, width=750,
        dirpath=dirpath,
        filename='pct_weight',
        save_html=TRUE
    )
}


fig <- plot_scatter(
    df,
    x='week', y='weight', group_by='group_name',
    ymin=0,
    xlabel='Time (week)', ylabel='Mouse Weight (g)', title='Weight Curves',
    mode='lines+markers',
    color_discrete_map=color_discrete_map,
    hover_data=c('sex', 'genotype', 'dob', 'day', 'week')
)

# save
if (!troubleshooting) {
    save_fig(
        fig=fig,
        height=350, width=750,
        dirpath=dirpath,
        filename='weight',
        save_html=TRUE
    )
}


fig <- plot_multiple_comparisons(
    weight_tbl,
    x='genotype', y='spleen_weight',
    ymin=0,
    ylabel="Spleen Weight (mg)",
    title='Spleen Weight',
    show_numbers=FALSE,
    test='fishers_lsd',
    custom_group_order=c('WT', 'KO')
)


# save
if (!troubleshooting) {
    withCallingHandlers({
        ggsave(
            file.path(wd, opt[['output-dir']], 'spleen_weight.png'),
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



fig <- plot_scatter(
    weight_tbl[(!is.na(weight_tbl[['spleen_weight']])), ],
    x='45884', y='spleen_weight', group_by='genotype',
    ymin=0,
    xlabel='Mouse Starting Weight (g)', ylabel='Spleen Weight (mg)', title='Spleen Weight vs. Mouse Weight',
    mode='markers',
    # color_discrete_map=color_discrete_map,
    hover_data=c('sex', 'genotype')
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


df_subs <- df[(df[['died']]==1) & (!is.na(df[['mouse_id']])), ]

fig <- plot_scatter(
    df_subs,
    x='week', y='spleen_weight', group_by='genotype',
    ymin=0,
    xlabel='Death Date (weeks)', ylabel='Spleen Weight (mg)',
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


df_subs[['spleen_weight_ratio_start']] <- df_subs[['spleen_weight']] / df_subs[['weight_start']]

fig <- plot_scatter(
    df_subs,
    x='week', y='spleen_weight_ratio_start', group_by='genotype',
    ymin=0,
    xlabel='Death Date (weeks)', ylabel='[Spleen Weight / Mouse Starting Weight] (mg/g)',
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


df_subs[['spleen_weight_ratio_last']] <- df_subs[['spleen_weight']] / df_subs[['weight']]

fig <- plot_scatter(
    df_subs,
    x='week', y='spleen_weight_ratio_last', group_by='genotype',
    ymin=0,
    xlabel='Death Date (weeks)', ylabel='[Spleen Weight / Mouse Final Weight] (mg/g)',
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



end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
