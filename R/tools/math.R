## Functions
## stats_test


#' Switch to select which test to use
#' 
#' Currently only includes t_test
#' Will update in the future to also include ANOVA
#' 
stats_test <- function(x, y, test='t_test') {

    # exit if either side is null or both sides have fewer observations
    if (is.null(x) | is.null(y)) {
        return(NA)
    }
    if ((length(x)==1) & (length(y)==1)) {
        return(NA)
    }

    if (test=='t_test') {

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
    }
}
