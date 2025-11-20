import pandas as pd
import numpy as np
from typing import List, Dict, Any, Union
from .list_tools import items_in_a_not_b, replace_specific_items


def append_dataframe(df1: pd.DataFrame, df2: pd.DataFrame,
                     infront: bool = False, reset_index_flag: bool = True) -> pd.DataFrame:
    """Append df2 to df1, handling missing columns."""
    missing_cols = items_in_a_not_b(list(df1.columns), list(df2.columns))
    for col in missing_cols:
        df2[col] = np.nan

    common_cols = list(set(df1.columns) & set(df2.columns))

    if infront:
        df = pd.concat([df2[common_cols], df1], ignore_index=reset_index_flag)
    else:
        df = pd.concat([df1, df2[common_cols]], ignore_index=reset_index_flag)

    return df


def dataframe_row_from_named_list(items: dict) -> pd.DataFrame:
    """Convert a dictionary into a single-row DataFrame."""
    return pd.DataFrame([items])


def fillna(df: pd.DataFrame, cols: List[str], val: Any = 0) -> pd.DataFrame:
    """Fill NA values in specific columns."""
    df = df.copy()
    for col in cols:
        df[col] = df[col].fillna(val)
    return df


def group_by_agg(df: pd.DataFrame, groups: List[str],
                 values: List[str], agg_func=np.mean) -> pd.DataFrame:
    """Group by and aggregate using base pandas."""
    return df.groupby(groups, as_index=False)[values].agg(agg_func)


def rename_columns(df: pd.DataFrame, columns: dict) -> pd.DataFrame:
    """Rename dataframe columns using a dictionary."""
    df = df.copy()
    df.columns = replace_specific_items(list(df.columns), columns)
    return df


def reset_index(df: pd.DataFrame, index_name: str = 'index', drop: bool = False) -> pd.DataFrame:
    """Reset index, optionally preserving the old index as a column."""
    df = df.copy()
    if drop:
        df = df.reset_index(drop=True)
    else:
        df = df.reset_index()
        if 'index' in df.columns:
            df = df.rename(columns={'index': index_name})
    return df


def stranspose(df: pd.DataFrame, colname: str = None) -> pd.DataFrame:
    """Transpose DataFrame, optionally using a column as new column names."""
    tdf = df.T.copy()
    if colname is not None:
        tdf.columns = tdf.loc[colname]
        tdf = tdf.drop(colname)
    return tdf


def pivot_then_collapse(df: pd.DataFrame, index_cols: List[str],
                        group_name: str, metric: str,
                        custom_group_order: List[str] = None) -> pd.DataFrame:
    """Pivot groups into columns, then collapse metrics into lists."""
    if custom_group_order:
        group_names = [g for g in custom_group_order if g in df[group_name].unique()]
    else:
        group_names = sorted(df[group_name].unique())

    # Pivot and collect values into lists
    subset = df[index_cols + [group_name, metric]].copy()
    pivoted = subset.pivot_table(
        index=index_cols,
        columns=group_name,
        values=metric,
        aggfunc=list
    ).reset_index()

    pivoted['metric'] = metric

    # Reorder columns
    cols = index_cols + ['metric'] + [g for g in group_names if g in pivoted.columns]
    pivoted = pivoted[[c for c in cols if c in pivoted.columns]]

    # Sort rows
    pivoted = pivoted.sort_values(by=index_cols[::-1]).reset_index(drop=True)

    return pivoted
