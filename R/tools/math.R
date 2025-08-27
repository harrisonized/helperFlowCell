import::here(tidyr, 'pivot_wider')
import::here(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'pivot_then_collapse',
    .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'collect_matrix_cols', 'matrix2list', 'fill_missing_keys',
    .character_only=TRUE)

## Functions
## unpaired_t_test
## apply_unpaired_t_test
## fishers_lsd
## tukey_multiple_comparisons
## bonferroni_multiple_comparisons
## apply_multiple_comparisons
## generate_normal_data
## generate_lognormal_data
## compute_normal_tvd
## compute_lognormal_tvd


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
        group_names <- intersect( custom_group_order, unique(df[[group_name]]) )
    } else {
        group_names <- sort(unique( df[[group_name]] ))
    }

    n_groups <- length(group_names)

    # exit if only one group
    if (n_groups==1) {
        return(NaN)
    }
    
    # mfi_tbl <- pivot_then_collapse(
    #     df,
    #     index_cols=c('organ', 'cell_type'),
    #     group_name='group_name',
    #     metric=opt[['metric']]
    # )
    res <- pivot_then_collapse(
        df, index_cols, group_name,
        metric, custom_group_order
    )

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
        group_names <- intersect( custom_group_order, unique(df[[group]]) )
    } else {
        group_names <- sort(unique( df[[group]] ))
    }

    n_groups <- length(group_names)

    # exit if only one group
    if (n_groups==1) {
        return(NA)
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

    return(pvals)
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
        group_names <- intersect( custom_group_order, unique(df[[group]]) )
        df[[group]] <- factor(df[[group]], levels = custom_group_order)
    } else {
        group_names <- sort(unique( df[[group]] ))
    }

    n_groups <- length(group_names)

    # exit if only one group
    if (n_groups==1) {
        return(NA)
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
        group_names <- intersect( custom_group_order, unique(df[[group]]) )
        df[[group]] <- factor(df[[group]], levels = custom_group_order)
    } else {
        group_names <- sort(unique( df[[group]] ))
    }

    n_groups <- length(unique( df[[group]] ))

    # exit if only one group
    if (n_groups==1) {
        return(NA)
    }

    pvals_mat <- pairwise.t.test(
        df[[metric]], df[[group]],
        p.adjust.method = "bonferroni"
    )[['p.value']]
    pvals <- matrix2list(pvals_mat)  # flatten
    # pvals <- Filter(Negate(is.na), pvals)  # filter NA

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

    res <- pivot_then_collapse(
        df, index_cols, group_name,
        metric, custom_group_order
    )

    # split dataframe to fit expected input
    df_list <- split(
        df[, c(index_cols, group_name, metric)],
        df[, c(index_cols)],
        drop=TRUE
    )

    if (correction=='fishers_lsd') {
        pval_list <- mapply(
            function(x) fishers_lsd(x, group=group_name, metric=metric),
            df_list, SIMPLIFY=FALSE
        )
    } else if (correction=='tukey') {
        pval_list <- mapply(
            function(x) tukey_multiple_comparisons(x, group=group_name, metric=metric),
            df_list, SIMPLIFY=FALSE
        )
    } else if (correction=='bonferroni') {
        pval_list <- mapply(
            function(x) bonferroni_multiple_comparisons(x, group=group_name, metric=metric),
            df_list, SIMPLIFY=FALSE
        )
    } else {
        stop("Choose correction='fishers_lsd', 'tukey' or 'bonferroni'")
    }

    # convert to matrix
    colnames <- unique(unlist(lapply(pval_list, names)))  # get all column names
    if (length(colnames)==1) {
        pvals <- matrix(pval_list, nrow=1, dimnames=list(colnames, names(pval_list)))
    } else if (is.null(dim(pval_list))) {
        pvals <- mapply(function(x) fill_missing_keys(x, colnames), pval_list)
    }
    pvals <- matrix(
        sapply(pvals, function(x) x[[1]] ),
        nrow = nrow(pvals), ncol = ncol(pvals), dimnames = dimnames(pvals)
    )

    res <- cbind(res, t(pvals))  # res and pvals are sorted for the same order
    res <- res[do.call(order, res[index_cols]), ]  # sort rows in original order

    return(res)
}


#' Generate Normal Data
#' 
#' Generate a dataframe of random values drawn from a normal distribution
#' Use this as an input to plot_modal_histograms
#' 
generate_normal_data <- function(
    n=1000, mean=200, sd=100,
    group_name="group"
) {
    df <- data.frame(group=group_name, value=rnorm(n, mean=mean, sd=sd))
    return(df)
}


#' Generate Lognormal Data
#' 
#' Use this function to simulate data for building fluorescence histograms
#' 
#' Flow intensity data is lognormal
#' See this blog post on how to correctly use the rlnorm function
#' https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r/
#' 
#' Special care had to be used to handle negative MFI values
#' When Flowjo reports a negative MFI, the actual intensity is 1/abs(MFI)
#' To correct for this, I convert the MFI values using the following formula:
#' adj_mean <- 10^(sign(MFI) * log10(1 + abs(MFI)))
#' This formula was chosen so that values in the range -1 < MFI < 1 are smooth and continuous
#' For very large positive values of MFI, the formula approaches adj_mean = MFI
#' For very negative values of MFI, the formula approaches adj_mean <- 1/MFI
#' 
generate_lognormal_data <- function(
  n=1000, mean=200, sd=100, group_name = "group"
) {

    if (mean==0) {
        adj_mean <- 1
    } else {
        # make the mean continuous
        signed_log <- sign(mean) * log10(1 + abs(mean))
        adj_mean <- 10^signed_log
    }
    adj_sd <- abs(sd/mean) * adj_mean

    location <- log(adj_mean^2 / sqrt(adj_sd^2 + adj_mean^2))
    shape <- sqrt(log(1 + (adj_sd / adj_mean)^2))

    # clamp mean to prevent errors
    if (location < -372.5666) {
        location <- -372.5666
    }

    values <- rlnorm(n, meanlog = location, sdlog = shape)
    df <- data.frame(group = group_name, value = values)

    return(df)
}


#' Compute Normal Total Variation Distance
#'
#' Returns a value between 0 and 1 that quantifies 1-overlap. This can be
#' interpreted as dissimilarity between two distributions. To enable this 
#' function handle to negative MFI values, the same log transform used for
#' generate_lognormal_data has been added. This reduces the differences
#' when negative MFI values are compared.
#' 
compute_normal_tvd <- function(mean1, sd1, mean2, sd2, log_transform=FALSE) {

    if (is.vector(mean1) || is.list(mean1)) {
        mean1 <- mean(mean1)
    }
    if (is.vector(sd1) || is.list(sd1)) {
        sd1 <- mean(unlist(sd1))
    }
    if (is.vector(mean2) || is.list(mean2)) {
        mean2 <- mean(unlist(mean2))
    }
    if (is.vector(sd2) || is.list(sd2)) {
        sd2 <- mean(unlist(sd2))
    }

    if (is.na(sd1) & is.na(sd2)) { return(NA) }
    if (is.na(mean2) & is.na(sd2)) { return(NA) }
    if (!is.na(sd1) & is.na(sd2)) { sd2 <- sd1 }
    if (is.na(sd1) & !is.na(sd2)) { sd1 <- sd2 }
    if (sd1==0 | sd2==0) { return(NA) }
    if (is.na(mean2)) { mean2 <- 0 }

    if (log_transform) {
        mean1 <- 10^(sign(mean1)*log10(1+abs(mean1))) /log(10)
        sd1 <- sd1/log(10)
        mean2 <- 10^(sign(mean2)*log10(1+abs(mean2))) /log(10)
        sd2<- sd2/log(10)
    }

    p <- function(x) dnorm(x, mean = mean1, sd = sd1)
    q <- function(x) dnorm(x, mean = mean2, sd = sd2)
    min_pdf <- function(x) pmin(p(x), q(x))  # Select minimum of two curves

    # Get area under the curve of the overlap
    # Set integral limit at 10*sd in each direction
    lower <- min(mean1 - 10 * sd1, mean2 - 10 * sd2)
    upper <- max(mean1 + 10 * sd1, mean2 + 10 * sd2)
    overlap <- integrate(min_pdf, lower, upper)$value

    dissimiliarity <- 1 - overlap
    return(dissimiliarity)
}


#' Compute Lognormal Total Variation Distance
#' 
#' Returns a value between 0 and 1 that quantifies 1-overlap. This can be
#' interpreted as dissimilarity between two distributions. Since this is
#' NOT computed on linear space, this tends to magnify differences when
#' one mean is large and the other is small. Therefore, compute_normal_tvd
#' is a better metric for measuring overlap for fluorescence data.
#' 
compute_lognormal_tvd <- function(mean1, sd1, mean2, sd2) {

    if (is.na(sd1) & is.na(sd2)) { return(NA) }
    if (is.na(mean2) & is.na(sd2)) { return(NA) }
    if (!is.na(sd1) & is.na(sd2)) { sd2 <- sd1 }
    if (is.na(sd1) & !is.na(sd2)) { sd1 <- sd2 }
    if (sd1==0 | sd2==0) { return(NA) }
    if (is.na(mean2)) { mean2 <- 0 }

    means <- c(mean1, mean2)
    sds <- c(sd1, sd2)

    curves <- list()
    for (idx in c(1, 2)) {

        mean <- means[[idx]]
        sd <- sds[[idx]]

        if (mean==0) {
            adj_mean <- 1
        } else {
            # make the mean continuous
            signed_log <- sign(mean) * log10(1 + abs(mean))
            adj_mean <- 10^signed_log
        }
        adj_sd <- abs(sd/mean) * adj_mean

        location <- log(adj_mean^2 / sqrt(adj_sd^2 + adj_mean^2))
        shape <- sqrt(log(1 + (adj_sd / adj_mean)^2))

        curves[[idx]] <- local({
            loc <- location
            shp <- shape
            function(x) dlnorm(x, meanlog = loc, sdlog = shp)
        })
    }

    overlap <- integrate(
        function(x) pmin(curves[[1]](x), curves[[2]](x)),
        lower = 0, upper = max(mean1 + 10 * sd1, mean2 + 10 * sd2),
        rel.tol = 1e-8
    )$value

    dissimiliarity <- 1 - overlap
    return(dissimiliarity)
}
