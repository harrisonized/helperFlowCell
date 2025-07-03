import::here(tidyr, 'pivot_wider')
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'collect_matrix_cols', 'matrix2list', .character_only=TRUE)

## Functions
## unpaired_t_test
## apply_unpaired_t_test
## fishers_lsd
## tukey_multiple_comparisons
## bonferroni_multiple_comparisons
## apply_multiple_comparisons


#' Unpaired T Test
#' 
#' @description Use this when you have multiple groups that are not dependent on each other
#' and the variances are clearly unequal between the groups
#' Try this as a second resort if fishers_lsd is giving too many false positives
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
    pval <- t.test(
        obs,
        alternative=(if (mean(obs) > ob) {"greater"} else {"less"}),
        mu=ob,
        var.equal=FALSE
    )[['p.value']]

    return(pval)
}


#' Apply Unpaired T Test
#' 
#' @description Apply unpaired t test to dataframe
#' 
apply_unpaired_t_test <- function(
    df,
    index_cols,  # c('organ', 'cell_type')
    group_name,  # 'group_name'
    metric,  # 'pct_cells'
    custom_group_order=c()
) {

    if (all(is.na( df[[metric]] )) ) {
        return(NA)
    }

    if (length(custom_group_order)>=1) {
        group_names <- custom_group_order
    } else {
        group_names <- sort(unique( df[[group_name]] ))
    }

    n_groups <- length(group_names)

    # exit if only one group
    if (n_groups==1) {
        return(NaN)
    }
    
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
    id_combos <- collect_matrix_cols(combn(1:n_groups, 2))  # generate pairs of indexes
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


#' One Way ANOVA with no correction, ie. Fisher's LSD
#' 
#' @description Use this when you have multiple groups that are not dependent on each other
#' and the variances are roughly equal between the groups.
#' Use this first. If the group with the largest variance is more than 4x the variance of
#' the smallest group, switch to unpaired_t_test
#' 
fishers_lsd <- function(
    df,
    group,  # 'group_name'
    metric,  # 'pct_cells'
    custom_group_order=c()
) {

    if (all(is.na( df[[metric]] )) ) {
        return(NA)
    }

    if (length(custom_group_order)>=1) {
        group_names <- custom_group_order
    } else {
        group_names <- sort(unique( df[[group]] ))
    }

    n_groups <- length(group_names)

    # exit if only one group
    if (n_groups==1) {
        return(NaN)
    }


    # ----------------------------------------------------------------------
    # ANOVA

    # computations
    formula <- as.formula(
        paste(deparse(as.name(metric)), "~", deparse(as.name(group)))  # 'metric ~ group'
    )
    group_means <- aggregate(formula, data = df, FUN = mean)
    group_means <- setNames(group_means[[metric]], group_means[[group]])
    group_sizes <- aggregate(formula, data = df, FUN = length)
    group_sizes <- setNames(group_sizes[[metric]], group_sizes[[group]])

    overall_mean <- mean(df[[metric]])
    ss_total <- sum((df[[metric]] - overall_mean)^2)
    ss_between <- sum(group_sizes * (group_means - overall_mean)^2)
    ss_within <- ss_total - ss_between

    df_between <- length(group_names) - 1  # degrees of freedom between
    df_within <- length(df[[metric]]) - length(group_names)  # degrees of freedom within
    mse <- ss_within / df_within  # mean square error


    # ----------------------------------------------------------------------
    # Fisher's LSD p values

    # iterate through all pairs of columns
    # collect results directly in pval_col
    id_combos <- collect_matrix_cols(combn(1:n_groups, 2))  # generate pairs of indexes
    pval_cols <- sapply(id_combos,  
        function(x) paste(group_names[x[[2]]], group_names[x[[1]]], sep='-')
    )  # generate name pairs
    pvals <- list()

    for (idx in 1:choose(n_groups, 2)) {
        left <- group_names[ id_combos[[idx]][1] ]  # 1st col
        right <- group_names[ id_combos[[idx]][2] ]  # 2st col

        # fishers lsd test
        mean_diff <- group_means[right] - group_means[left]
        se_diff <- sqrt(mse * (1/group_sizes[left] + 1/group_sizes[right]))
        t_stat <- mean_diff / se_diff
        pval <- 2 * pt(-abs(t_stat), df_within)

        pvals[ pval_cols[idx] ] <- pval
    }

    return(unlist(pvals))
}


#' One Way ANOVA with Tukey correction
#' 
#' @description Use this when you have multiple groups that are may be dependent on each other
#' and the variances are roughly equal between the groups. This is supposed to control for false
#' positives due to group effects. In general, this will provide larger p values
#' (ie. less significance) compared with fishers_lsd.
#' 
tukey_multiple_comparisons <- function(
    df,
    group,  # 'group_name',
    metric,  # 'pct_cells'
    custom_group_order=c()
) {

    if (all(is.na( df[[metric]] )) ) {
        return(NA)
    }

    if (length(custom_group_order)>=1) {
        group_names <- custom_group_order
        df[[group]] <- factor(df[[group]], levels = custom_group_order)
    } else {
        group_names <- sort(unique( df[[group]] ))
    }

    n_groups <- length(group_names)

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
        pvals <- setNames(pvals, paste(group_names[2],group_names[1], sep='-'))
    }

    return(pvals)
}


#' One Way ANOVA with Bonferroni correction
#' 
#' @description Use this when you have multiple groups that are heavily dependent on each other
#' and there are a large number of false positives because of group interactions. Since this is
#' much more stringent than tukey_multiple_comparisons, this is not recommended due to the
#' potential to have false negatives, especially for small datasets.
#' 
bonferroni_multiple_comparisons <- function(
    df,
    group,  # 'group_name'
    metric,  # 'pct_cells'
    custom_group_order=c()
) {

    if (all(is.na( df[[metric]] )) ) {
        return(NA)
    }

    if (length(custom_group_order)>=1) {
        group_names <- custom_group_order
        df[[group]] <- factor(df[[group]], levels = custom_group_order)
    } else {
        group_names <- sort(unique( df[[group]] ))
    }

    n_groups <- length(unique( df[[group]] ))

    # exit if only one group
    if (n_groups==1) {
        return(NaN)
    }

    pvals_mat <- pairwise.t.test(
        df[[metric]], df[[group]],
        p.adjust.method = "bonferroni"
    )[['p.value']]
    pvals <- matrix2list(pvals_mat)  # flatten
    pvals <- Filter(Negate(is.na), pvals)  # filter NA

    return(pvals)
}


#' Apply Multiple Comparisons
#' 
#' @description Applies one-way ANOVA across a dataframe
#' Currently, this implementation skips the ANOVA parameters and instead just
#' provides the p_values of the post-hoc analysis for plotting
#' In the future, the ANOVA F-statistic may be provided
#' 
apply_multiple_comparisons <- function(
    df,
    index_cols,  # c('organ', 'cell_type')
    group_name,  # 'group_name'
    metric,  # 'pct_cells' or 'abs_count'
    correction='fishers_lsd',  # 'fishers_lsd', 'tukey', or 'bonferroni'
    custom_group_order=c()
) {
    
    if (length(custom_group_order)>=1) {
        group_names <- custom_group_order
    } else {
        group_names <- sort(unique( df[[group_name]] ))
    }

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
    

    # split dataframe to fit expected input
    df_list <- split(
        df[, c(index_cols, group_name, metric)],
        df[, c(index_cols)],
        drop=TRUE
    )

    if (correction=='fishers_lsd') {
        pvals <- mapply(
            function(x) fishers_lsd(x, group=group_name, metric=metric),
            df_list
        )
    } else if (correction=='tukey') {
        pvals <- mapply(
            function(x) tukey_multiple_comparisons(x, group=group_name, metric=metric),
            df_list
        )
    } else if (correction=='bonferroni') {
        pvals <- mapply(
            function(x) bonferroni_multiple_comparisons(x, group=group_name, metric=metric),
            df_list
        )
    } else {
        stop("Choose correction='fishers_lsd', 'tukey' or 'bonferroni'")
    }

    # convert to matrix
    if (is.null(dim(pvals))) {
        colnames <- unique(unlist(lapply(pvals, names)))
        pvals  <- mapply(function(x) fill_missing_keys(x, colnames), pvals)
    }

    res <- cbind(res, t(pvals))  # res and pvals are sorted for the same order
    res <- res[do.call(order, res[index_cols]), ]    # sort rows in original order

    return(res)
}
