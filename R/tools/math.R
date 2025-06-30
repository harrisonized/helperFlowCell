import::here(tidyr, 'pivot_wider')
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'flatten_matrix', .character_only=TRUE)

## Functions
## get_significance_code
## unpaired_t_test
## apply_unpaired_t_test
## tukey_multiple_comparisons
## apply_tukey_multiple_comparisons


#' Get Significance Code
#' 
#' Convert p value to significance
#' 
get_significance_code <- function(p) {

    p <- abs(p)

    if (p > 0.05) {
        return('n.s.')
    } else if (p > 0.01) {
        return('*')
    } else if (p > 0.001) {
        return('**')
    } else if (p > 0.0001 ) {
        return('***')
    } else {
        return('****')
    }
}


#' Unpaired T Test
#' 
unpaired_t_test <- function(x, y) {

    # exit if either side is null or both sides have fewer observations
    if (is.null(x) | is.null(y)) {
        return(NA)
    }
    if ((length(x)==1) & (length(y)==1)) {
        return(NA)
    }

    # two.sided if x and y both have multiple observations
    if ((length(x)>1) & (length(y)>1)) {
        pval <- t.test(x, y, alternative='two.sided', var.equal=FALSE)[['p.value']]
        return(pval)
    }

    # if one side only has a single observation, use one.sided t test 
    if (length(x)==1) {
        ob <- x  # single observation
        obs <- y
    } else if (length(y)==1) {
        obs <- x 
        ob <- y  # single observation
    }
    pval <- t.test(obs, alternative=(if (mean(obs) > ob) {"greater"} else {"less"}), mu=ob)[['p.value']]

    return(pval)
}


#' Apply Unpaired T Test
#' 
#' Apply unpaired t test to dataframe
#' 
apply_unpaired_t_test <- function(
    df,
    index_cols,  # c('organ', 'cell_type')
    group_name,  # 'group_name'
    metric  # 'pct_cells'
) {

    group_names <- sort(unique(df[[group_name]]))
    n_groups <- length(group_names)

    # Collect values into list columns
    res <- pivot_wider(
        df[, c(index_cols, group_name, metric)],
        names_from = group_name,
        values_from = metric,
        values_fn = list,  # suppress warning
        names_glue = "{.name}"
    )
    res[['metric']] <- metric
    res <- res[do.call(order, res[index_cols]), ]  # sort rows
    res <- res[, c(index_cols, 'metric', group_names)]  # sort cols


    # iterate through all pairs of columns
    # collect results directly in pval_col
    id_combos <- flatten_matrix(combn(1:n_groups, 2))  # generate pairs of indexes
    pval_cols <- sapply(id_combos,  
        function(x) paste(group_names[x[[2]]], group_names[x[[1]]], sep='-')
    )  # generate name pairs
    for (idx in 1:choose(n_groups, 2)) {
        idx1 <- id_combos[[idx]][1]  # 1st col
        idx2 <- id_combos[[idx]][2]  # 2nd col

        # t test
        res[ pval_cols[idx] ] <- mapply(
            function(x, y) unpaired_t_test(x, y),
            res[[ group_names[idx1] ]],  # 1st col
            res[[ group_names[idx2] ]]   # 2nd col
        )
    }

    return(res)
}


#' One Way ANOVA with Tukey
#' 
#' Returns a list of p values
#' Uses the Tukey correction (recommended by GraphPad)
#' 
tukey_multiple_comparisons <- function(df, group='group_name', metric='pct_cells') {

    n_groups <- length(unique(df[[group]]))

    # exit if only one group
    if (n_groups==1) {
        return(NaN)
    }

    formula <- as.formula(
        paste(deparse(as.name(metric)), "~", deparse(as.name(group)))  # 'metric ~ group'
    )
    withCallingHandlers({
        fit <- aov(formula, data = df)
        pvals <- TukeyHSD(fit)[[group]][, 'p adj']
    }, warning = function(w) {
        if ( any(grepl("NaNs produced", w)) ) {
            invokeRestart("muffleWarning")
        }
    })

    if (n_groups==2) {
        groups <- unique(df[['group_name']])
        pvals <- setNames(pvals, paste(groups[2],groups[1], sep='-'))
    }

    return(pvals)
}


#' Apply Tukey multiple comparisons
#' 
#' Apply Tukey multiple comparisons to all groups in a dataframe
#' 
apply_tukey_multiple_comparisons <- function(
    df,
    index_cols,  # c('organ', 'cell_type')
    group_name,  # 'group_name'
    metric  # 'pct_cells'
) {

    group_names <- sort(unique(df[[group_name]]))

    # Collect values into list columns
    res <- pivot_wider(
        df[, c(index_cols, group_name, metric)],
        names_from = group_name,
        values_from = metric,
        values_fn = list,  # suppress warning
        names_glue = "{.name}"
    )
    res[['metric']] <- metric
    res <- res[do.call(order, res[rev(index_cols)]), ]  # sort rows
    res <- res[, c(index_cols, 'metric', group_names)]  # sort cols
    

    # split dataframe to fit expected input for tukey_multiple_comparisons
    df_list <- split(
        df[, c(index_cols, group_name, metric)],
        df[, c(index_cols)],
        drop=TRUE
    )
    pvals <- mapply(
        function(x) tukey_multiple_comparisons(x, group='group_name', metric='pct_cells'),
        df_list
    )
    colnames <- unique(unlist(lapply(pvals, names)))
    pvals  <- mapply(function(x) fill_missing_keys(x, colnames), pvals)
    res <- cbind(res, t(pvals))  # res and pvals are sorted for the same order

    res <- res[do.call(order, res[index_cols]), ]    # sort rows in original order
    return(res)
}
