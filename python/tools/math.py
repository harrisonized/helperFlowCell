import numpy as np
import pandas as pd
from scipy import stats
from itertools import combinations
from typing import List, Union
from .df_tools import pivot_then_collapse
from .list_tools import collect_matrix_cols, fill_missing_keys


def unpaired_t_test(x, y) -> float:
    """Perform unpaired t-test between two groups."""
    if x is None or y is None:
        return np.nan

    x = np.array(x).flatten() if hasattr(x, '__iter__') else np.array([x])
    y = np.array(y).flatten() if hasattr(y, '__iter__') else np.array([y])

    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]

    if len(x) == 0 or len(y) == 0:
        return np.nan
    if len(x) == 1 and len(y) == 1:
        return np.nan

    if len(x) > 1 and len(y) > 1:
        _, pval = stats.ttest_ind(x, y, equal_var=False)
        return pval

    # One-sided t-test when one group has single observation
    if len(x) == 1:
        ob, obs = x[0], y
    else:
        ob, obs = y[0], x

    alternative = 'greater' if np.mean(obs) > ob else 'less'
    _, pval = stats.ttest_1samp(obs, ob)
    return pval / 2 if alternative == 'less' else pval / 2


def apply_unpaired_t_test(df: pd.DataFrame, index_cols: List[str],
                          group_name: str, metric: str,
                          custom_group_order: List[str] = None) -> pd.DataFrame:
    """Apply unpaired t-test across a dataframe."""
    if df[metric].isna().all():
        return None

    if custom_group_order:
        group_names = [g for g in custom_group_order if g in df[group_name].unique()]
    else:
        group_names = sorted(df[group_name].unique())

    n_groups = len(group_names)
    if n_groups == 1:
        return None

    res = pivot_then_collapse(df, index_cols, group_name, metric, custom_group_order)

    # Generate pairs and compute p-values
    pairs = list(combinations(range(n_groups), 2))
    for idx1, idx2 in pairs:
        col_name = f"{group_names[idx2]}-{group_names[idx1]}"
        pvals = []
        for i in range(len(res)):
            x = res[group_names[idx1]].iloc[i]
            y = res[group_names[idx2]].iloc[i]
            pvals.append(unpaired_t_test(x, y))
        res[col_name] = pvals

    return res


def fishers_lsd(df: pd.DataFrame, group: str, metric: str,
                custom_group_order: List[str] = None) -> dict:
    """One-way ANOVA with Fisher's LSD (no correction)."""
    if df[metric].isna().all():
        return {}

    if custom_group_order:
        group_names = [g for g in custom_group_order if g in df[group].unique()]
    else:
        group_names = sorted(df[group].unique())

    n_groups = len(group_names)
    if n_groups == 1:
        return {}

    # ANOVA computations
    groups_data = [df[df[group] == g][metric].dropna().values for g in group_names]
    group_means = {g: np.mean(groups_data[i]) for i, g in enumerate(group_names)}
    group_sizes = {g: len(groups_data[i]) for i, g in enumerate(group_names)}

    overall_mean = df[metric].mean()
    ss_between = sum(group_sizes[g] * (group_means[g] - overall_mean)**2 for g in group_names)
    ss_total = ((df[metric] - overall_mean)**2).sum()
    ss_within = ss_total - ss_between

    df_within = len(df[metric].dropna()) - n_groups
    mse = ss_within / df_within

    # Fisher's LSD p-values
    pvals = {}
    for idx1, idx2 in combinations(range(n_groups), 2):
        left, right = group_names[idx1], group_names[idx2]
        col_name = f"{right}-{left}"

        mean_diff = group_means[right] - group_means[left]
        se_diff = np.sqrt(mse * (1/group_sizes[left] + 1/group_sizes[right]))
        t_stat = mean_diff / se_diff
        pval = 2 * stats.t.sf(abs(t_stat), df_within)
        pvals[col_name] = pval

    return pvals


def tukey_multiple_comparisons(df: pd.DataFrame, group: str, metric: str,
                               custom_group_order: List[str] = None) -> dict:
    """One-way ANOVA with Tukey HSD correction."""
    from scipy.stats import tukey_hsd

    if df[metric].isna().all():
        return {}

    if custom_group_order:
        group_names = [g for g in custom_group_order if g in df[group].unique()]
        df = df[df[group].isin(group_names)].copy()
    else:
        group_names = sorted(df[group].unique())

    n_groups = len(group_names)
    if n_groups == 1:
        return {}

    groups_data = [df[df[group] == g][metric].dropna().values for g in group_names]

    try:
        result = tukey_hsd(*groups_data)
        pvals = {}

        for i, j in combinations(range(n_groups), 2):
            col_name = f"{group_names[j]}-{group_names[i]}"
            pvals[col_name] = result.pvalue[i, j]

        return pvals
    except Exception:
        return fishers_lsd(df, group, metric, custom_group_order)


def bonferroni_multiple_comparisons(df: pd.DataFrame, group: str, metric: str,
                                    custom_group_order: List[str] = None) -> dict:
    """One-way ANOVA with Bonferroni correction."""
    if df[metric].isna().all():
        return {}

    if custom_group_order:
        group_names = [g for g in custom_group_order if g in df[group].unique()]
        df = df[df[group].isin(group_names)].copy()
    else:
        group_names = sorted(df[group].unique())

    n_groups = len(group_names)
    if n_groups == 1:
        return {}

    # Pairwise t-tests with Bonferroni correction
    from scipy.stats import ttest_ind

    n_comparisons = n_groups * (n_groups - 1) // 2
    pvals = {}

    for i, j in combinations(range(n_groups), 2):
        g1 = df[df[group] == group_names[i]][metric].dropna()
        g2 = df[df[group] == group_names[j]][metric].dropna()

        if len(g1) > 0 and len(g2) > 0:
            _, p = ttest_ind(g1, g2)
            col_name = f"{group_names[j]}-{group_names[i]}"
            pvals[col_name] = min(p * n_comparisons, 1.0)

    return pvals


def apply_multiple_comparisons(df: pd.DataFrame, index_cols: List[str],
                               group_name: str, metric: str,
                               correction: str = 'fishers_lsd',
                               custom_group_order: List[str] = None) -> pd.DataFrame:
    """Apply multiple comparisons across a dataframe."""
    res = pivot_then_collapse(df, index_cols, group_name, metric, custom_group_order)

    # Split dataframe for each unique combination of index columns
    grouped = df.groupby(index_cols)

    if correction == 'fishers_lsd':
        func = fishers_lsd
    elif correction == 'tukey':
        func = tukey_multiple_comparisons
    elif correction == 'bonferroni':
        func = bonferroni_multiple_comparisons
    else:
        raise ValueError("Choose correction='fishers_lsd', 'tukey', or 'bonferroni'")

    pval_list = []
    for name, group_df in grouped:
        pvals = func(group_df, group=group_name, metric=metric,
                     custom_group_order=custom_group_order)
        pval_list.append(pvals)

    # Convert to DataFrame columns
    if pval_list:
        all_cols = set()
        for p in pval_list:
            all_cols.update(p.keys())

        for col in all_cols:
            res[col] = [p.get(col, np.nan) for p in pval_list]

    return res


def generate_normal_data(n: int = 1000, mean: float = 200,
                         sd: float = 100, group_name: str = "group") -> pd.DataFrame:
    """Generate DataFrame with random normal values."""
    return pd.DataFrame({
        'group': [group_name] * n,
        'value': np.random.normal(mean, sd, n)
    })


def generate_lognormal_data(n: int = 1000, mean: float = 200,
                            sd: float = 100, group_name: str = "group") -> pd.DataFrame:
    """Generate DataFrame with lognormal values (for fluorescence simulation)."""
    if mean == 0:
        adj_mean = 1
    else:
        signed_log = np.sign(mean) * np.log10(1 + abs(mean))
        adj_mean = 10 ** signed_log

    adj_sd = abs(sd / mean) * adj_mean if mean != 0 else sd

    location = np.log(adj_mean**2 / np.sqrt(adj_sd**2 + adj_mean**2))
    shape = np.sqrt(np.log(1 + (adj_sd / adj_mean)**2))

    # Clamp location
    if location < -372.5666:
        location = -372.5666

    values = np.random.lognormal(location, shape, n)

    return pd.DataFrame({
        'group': [group_name] * n,
        'value': values
    })


def compute_normal_tvd(mean1, sd1, mean2, sd2, log_transform: bool = False) -> float:
    """Compute Total Variation Distance between two normal distributions."""
    # Handle list inputs
    if hasattr(mean1, '__iter__') and not isinstance(mean1, str):
        mean1 = np.mean(mean1)
    if hasattr(sd1, '__iter__') and not isinstance(sd1, str):
        sd1 = np.mean(sd1)
    if hasattr(mean2, '__iter__') and not isinstance(mean2, str):
        mean2 = np.mean(mean2)
    if hasattr(sd2, '__iter__') and not isinstance(sd2, str):
        sd2 = np.mean(sd2)

    if pd.isna(sd1) and pd.isna(sd2):
        return np.nan
    if pd.isna(mean2) and pd.isna(sd2):
        return np.nan
    if not pd.isna(sd1) and pd.isna(sd2):
        sd2 = sd1
    if pd.isna(sd1) and not pd.isna(sd2):
        sd1 = sd2
    if sd1 == 0 or sd2 == 0:
        return np.nan
    if pd.isna(mean2):
        mean2 = 0

    if log_transform:
        mean1 = 10**(np.sign(mean1) * np.log10(1 + abs(mean1))) / np.log(10)
        sd1 = sd1 / np.log(10)
        mean2 = 10**(np.sign(mean2) * np.log10(1 + abs(mean2))) / np.log(10)
        sd2 = sd2 / np.log(10)

    # Compute overlap using numerical integration
    from scipy.integrate import quad

    def p(x):
        return stats.norm.pdf(x, mean1, sd1)

    def q(x):
        return stats.norm.pdf(x, mean2, sd2)

    def min_pdf(x):
        return min(p(x), q(x))

    lower = min(mean1 - 10*sd1, mean2 - 10*sd2)
    upper = max(mean1 + 10*sd1, mean2 + 10*sd2)

    overlap, _ = quad(min_pdf, lower, upper)
    dissimilarity = np.sign(mean1 - mean2) * (1 - overlap)

    return dissimilarity


def compute_lognormal_tvd(mean1, sd1, mean2, sd2) -> float:
    """Compute TVD between two lognormal distributions."""
    if pd.isna(sd1) and pd.isna(sd2):
        return np.nan
    if pd.isna(mean2) and pd.isna(sd2):
        return np.nan
    if not pd.isna(sd1) and pd.isna(sd2):
        sd2 = sd1
    if pd.isna(sd1) and not pd.isna(sd2):
        sd1 = sd2
    if sd1 == 0 or sd2 == 0:
        return np.nan
    if pd.isna(mean2):
        mean2 = 0

    from scipy.integrate import quad

    params = []
    for mean, sd in [(mean1, sd1), (mean2, sd2)]:
        if mean == 0:
            adj_mean = 1
        else:
            signed_log = np.sign(mean) * np.log10(1 + abs(mean))
            adj_mean = 10 ** signed_log

        adj_sd = abs(sd / mean) * adj_mean if mean != 0 else sd
        location = np.log(adj_mean**2 / np.sqrt(adj_sd**2 + adj_mean**2))
        shape = np.sqrt(np.log(1 + (adj_sd / adj_mean)**2))
        params.append((location, shape))

    def curve1(x):
        return stats.lognorm.pdf(x, params[0][1], scale=np.exp(params[0][0]))

    def curve2(x):
        return stats.lognorm.pdf(x, params[1][1], scale=np.exp(params[1][0]))

    upper = max(mean1 + 10*sd1, mean2 + 10*sd2)
    overlap, _ = quad(lambda x: min(curve1(x), curve2(x)), 0, upper)

    dissimilarity = np.sign(mean1 - mean2) * (1 - overlap)
    return dissimilarity
