## Functions
## one_way_anova_with_tukey
## unpaired_t_test


#' One Way ANOVA with Tukey
#' 
#' Returns a list of p values
#' 
one_way_anova_with_tukey <- function(df, group, metric) {

    # exit if only one group
    if (length(unique(df[[group]]))==1) {
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

    return(pvals)
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
