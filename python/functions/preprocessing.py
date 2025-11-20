import pandas as pd
import numpy as np
import re
from typing import List, Optional
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config.replacements import (fluorophore_replacements, antibody_replacements,
                                  instr_cfg_colreps, ab_inv_colreps)
from config.flow import (id_cols, initial_gates, cell_type_spell_check,
                         cell_type_ignore)
from tools.df_tools import append_dataframe, rename_columns, reset_index
from tools.text_tools import title_to_snake_case
from tools.list_tools import multiple_replacement, find_first_match_index, items_in_a_not_b


def create_survival_table(start_info: pd.DataFrame, end_info: pd.DataFrame,
                          end_date) -> pd.DataFrame:
    """Reshape data for survfit analysis."""
    df = start_info.merge(end_info, on=['mouse_id', 'group'],
                          how='outer', suffixes=('_start', '_end'))

    # Calculate day from date
    df['day'] = (pd.to_datetime(df['date_end']) -
                 pd.to_datetime(df['date_start'])).dt.days
    df['week'] = df['day'] / 7
    df['is_dead'] = (df['date_end'] < end_date).astype(int)

    # Extend line to experiment start for survfit
    dummy_start = pd.DataFrame({
        'group': df['group'].unique(),
        'date_start': df['date_start'].min(),
        'date_end': df['date_start'].min(),
        'day': 0,
        'week': 0,
        'is_dead': 0
    })

    dummy_end = df.groupby('group').apply(
        lambda x: x.loc[x['date_end'].idxmin()]
    ).reset_index(drop=True)[['group', 'date_end', 'day', 'week']].copy()
    dummy_end['date_start'] = df['date_start'].min()
    dummy_end['day'] = dummy_end['day'] - 0.001
    dummy_end['week'] = dummy_end['week'] - 0.001
    dummy_end['is_dead'] = 0

    df = append_dataframe(df, dummy_start)
    df = append_dataframe(df, dummy_end)

    return df


def preprocess_flowjo_export(raw_table: pd.DataFrame, metric_name: str = 'num_cells',
                             include_initial_gates: bool = False) -> pd.DataFrame:
    """Preprocess FlowJo export table."""
    # Rename columns
    col_renames = {
        'X1': 'fcs_name',
        'Count': 'Ungated',
        'Mean (Comp-Alexa Fluor 488-A)': 'Ungated',
        'Mode (Comp-Alexa Fluor 488-A)': 'Ungated',
        'Geometric Mean (Comp-Alexa Fluor 488-A)': 'Ungated'
    }
    raw_table = rename_columns(raw_table, col_renames)

    # Clean column names
    raw_table.columns = [col.split(' | ')[0] for col in raw_table.columns]

    # Drop summary rows
    raw_table = raw_table[~raw_table['fcs_name'].isin(['Mean', 'SD'])]
    raw_table = raw_table[~raw_table['fcs_name'].str.contains('unstained', case=False, na=False)]
    raw_table = raw_table.reset_index(drop=True)

    # Determine columns to pivot
    exclude_cols = id_cols + ['Count']
    if not include_initial_gates:
        exclude_cols.extend(initial_gates)

    value_cols = items_in_a_not_b(list(raw_table.columns), exclude_cols)

    # Reshape to long format
    df = raw_table.melt(
        id_vars=[c for c in ['fcs_name'] if c in raw_table.columns],
        value_vars=value_cols,
        var_name='gate',
        value_name=metric_name
    ).dropna(subset=[metric_name])

    # Extract cell type from gate
    df['cell_type'] = df['gate'].apply(lambda x: x.split('/')[-1])

    # Apply spell check
    df['cell_type'] = df['cell_type'].replace(cell_type_spell_check)

    # Filter ignored cell types
    for cell_type in cell_type_ignore + ['mNeonGreen+']:
        df = df[~df['cell_type'].str.contains(cell_type, regex=False, na=False)]

    return df


def preprocess_antibody_inventory(df: pd.DataFrame) -> pd.DataFrame:
    """Preprocess antibody inventory table."""
    # Find first column with '...' pattern and filter
    pattern_cols = [i for i, col in enumerate(df.columns)
                    if re.search(r'\.{3}\d{2}', str(col))]
    if pattern_cols:
        df = df.iloc[:, :pattern_cols[0]]

    # Clean column names
    df.columns = [title_to_snake_case(col).replace('.', '') for col in df.columns]
    df = rename_columns(df, ab_inv_colreps)

    # Remove non-breaking space
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].replace('\u00a0', '', regex=True)
            df.loc[df[col] == '\u00a0', col] = np.nan

    # Standardize names
    if 'alternative_name' in df.columns:
        df['alternative_name'] = multiple_replacement(
            df['alternative_name'].fillna('').tolist(), antibody_replacements
        )

    df['antibody'] = multiple_replacement(
        df['antibody'].fillna('').tolist(), antibody_replacements
    )
    df['fluorophore'] = multiple_replacement(
        df['fluorophore'].fillna('').tolist(), fluorophore_replacements
    )

    # Fill missing fixable dyes
    fixable_dyes = ['DAPI', 'Zombie UV', 'Zombie Aqua']
    mask = (df['antibody'].isin(fixable_dyes)) & (df['fluorophore'].isna())
    df.loc[mask, 'fluorophore'] = df.loc[mask, 'antibody']

    return df


def preprocess_instrument_config(df: pd.DataFrame) -> pd.DataFrame:
    """Preprocess instrument configuration table."""
    # Clean column names
    df.columns = [title_to_snake_case(col) for col in df.columns]
    df = rename_columns(df, instr_cfg_colreps)

    # Standardize fluorophore names
    df['fluorophore'] = multiple_replacement(
        df['fluorophore'].fillna('').tolist(), fluorophore_replacements
    )

    return df


def sort_groups_by_metric(df: pd.DataFrame, x: str = 'cell_type',
                          y: str = 'pct_cells',
                          groups: List[str] = None) -> pd.DataFrame:
    """Sort rows by groups in decreasing order of the median metric."""
    if groups is None:
        groups = ['group_name']

    # Find median of each group
    x_axis_order = (df.groupby(x)[y]
                    .median()
                    .sort_values(ascending=False)
                    .index.tolist())

    # Sort
    df = df.copy()
    df['_sort_order'] = df[x].apply(lambda val: x_axis_order.index(val)
                                     if val in x_axis_order else len(x_axis_order))
    df = df.sort_values(['_sort_order'] + groups).drop('_sort_order', axis=1)

    return df.reset_index(drop=True)


def spillover_to_xml(spillover_table: pd.DataFrame, channels: list,
                     name: str = 'test') -> str:
    """Convert spillover table to XML format for FlowJo import."""
    from xml.etree.ElementTree import Element, SubElement, tostring
    from xml.dom import minidom

    root = Element('gating:gatingML')

    matrix = SubElement(root, 'transforms:spilloverMatrix')
    matrix.set('spectral', '0')
    matrix.set('weightOptAlgorithmType', 'OLS')
    matrix.set('prefix', 'Comp-')
    matrix.set('name', name)
    matrix.set('editable', '1')
    matrix.set('color', '#00ccff')
    matrix.set('version', 'FlowJo-10.10.0')
    matrix.set('status', 'FINALIZED')
    matrix.set('transforms:id', '')
    matrix.set('suffix', '')

    # Parameters
    params = SubElement(matrix, 'data-type:parameters')
    for channel in channels:
        channel_name = channel.split(' :: ')[0]
        param = SubElement(params, 'data-type:parameter')
        param.set('data-type:name', channel_name)
        param.set('userProvidedCompInfix', f'Comp-{channel_name}')

    # Spillover values
    for channel in channels:
        channel_name = channel.split(' :: ')[0]

        spillover = SubElement(matrix, 'transforms:spillover')
        spillover.set('data-type:parameter', channel_name)
        spillover.set('userProvidedCompInfix', f'Comp-{channel_name}')

        subset = spillover_table[spillover_table['- % Fluorochrome'] == channel]

        for _, row in subset.iterrows():
            coef = SubElement(spillover, 'transforms:coefficient')
            coef.set('data-type:parameter', row['Fluorochrome'].split(' :: ')[0])
            coef.set('transforms:value', str(row['Spectral Overlap']))

    # Pretty print
    xml_str = minidom.parseString(tostring(root)).toprettyxml(indent='  ')
    return xml_str


def parse_flowjo_metadata(df: pd.DataFrame,
                          cols: List[str] = None) -> pd.DataFrame:
    """Parse metadata from FCS filename (deprecated)."""
    if cols is None:
        cols = ['organ', 'mouse_id', 'treatment_group', 'strain']

    df = df.copy()
    df['metadata'] = df['fcs_name'].apply(
        lambda x: '_'.join(x.split('_')[3:-1])
    )

    # Handle sex column
    if 'sex' in df.columns:
        if df['sex'].dtype == bool:
            df['sex'] = df['sex'].map({False: 'F', True: 'T'})

    for i, col in enumerate(cols):
        try:
            df[col] = df['metadata'].apply(lambda x: x.split('-')[i])
        except Exception:
            print(f"Column {col} not in metadata")

    df['strain'] = df['mouse_id'].apply(lambda x: x.split('_')[0])

    return df
