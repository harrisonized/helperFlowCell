import os
import pandas as pd
from typing import List, Union
from .list_tools import filter_list_for_match


def append_many_csv(dir_path: str, sep: str = ',', recursive: bool = False,
                    include_filepath: bool = True, return_list: bool = False,
                    na_strings: List[str] = None) -> Union[pd.DataFrame, dict]:
    """Read all CSV files from a directory and append them into a single DataFrame."""
    if recursive:
        filenames = []
        for root, dirs, files in os.walk(dir_path):
            for f in files:
                if f.endswith('.csv'):
                    filenames.append(os.path.join(root, f))
    else:
        filenames = [os.path.join(dir_path, f) for f in os.listdir(dir_path)
                     if f.endswith('.csv')]

    if not filenames:
        return {} if return_list else None

    dfs = {}
    for filepath in filenames:
        try:
            df = pd.read_csv(filepath, sep=sep, na_values=na_strings)

            # Remove all-NA columns and rows
            df = df.dropna(axis=1, how='all')
            df = df.dropna(axis=0, how='all')

            # Repair empty column names
            empty_cols = [i for i, col in enumerate(df.columns) if not col or col.startswith('Unnamed')]
            for i, idx in enumerate(empty_cols):
                cols = list(df.columns)
                cols[idx] = f'X{i+1}'
                df.columns = cols

            # Add filename info
            if include_filepath:
                rel_path = os.path.relpath(filepath, dir_path)
                df.insert(0, 'filepath', rel_path)
                df.insert(1, 'filename', os.path.basename(filepath))

            dfs[filepath] = df
        except Exception as e:
            print(f"Error reading {filepath}: {e}")

    if return_list:
        return dfs
    else:
        if dfs:
            combined = pd.concat(dfs.values(), ignore_index=True)
            return combined
        return None


def list_files(dir_path: str, ext: str = None, recursive: bool = True) -> List[str]:
    """List all files with a specific extension."""
    all_files = []

    if recursive:
        for root, dirs, files in os.walk(dir_path):
            for f in files:
                all_files.append(os.path.join(root, f))
    else:
        all_files = [os.path.join(dir_path, f) for f in os.listdir(dir_path)
                     if os.path.isfile(os.path.join(dir_path, f))]

    if ext:
        all_files = [f for f in all_files if f.endswith(f'.{ext}')]

    return all_files


def join_many_csv(paths: Union[str, List[str]], index_cols: List[str],
                  value_cols: List[str] = None, all_x: bool = False,
                  all_y: bool = False, recursive: bool = True) -> pd.DataFrame:
    """Read CSV/TSV files and left join them column-wise."""
    # Get file list
    if isinstance(paths, str):
        if os.path.isdir(paths):
            paths = list_files(paths, recursive=recursive)

    if not paths or not any(os.path.exists(p) for p in paths):
        raise ValueError("No files found!")

    # Split by file type
    csv_paths = filter_list_for_match(paths, 'csv')
    tsv_paths = filter_list_for_match(paths, 'tsv')

    df_list = []
    filenames = []

    for fp in csv_paths:
        df_list.append(pd.read_csv(fp, sep=','))
        filenames.append(os.path.splitext(os.path.basename(fp))[0])

    for fp in tsv_paths:
        df_list.append(pd.read_csv(fp, sep='\t'))
        filenames.append(os.path.splitext(os.path.basename(fp))[0])

    if not df_list:
        raise ValueError("No valid CSV/TSV files found!")

    # Merge all dataframes
    cols_to_use = index_cols + (value_cols or [])
    result = df_list[0][cols_to_use].copy()

    for i, df in enumerate(df_list[1:], 1):
        how = 'outer' if (all_x and all_y) else ('left' if all_x else ('right' if all_y else 'inner'))
        result = result.merge(df[cols_to_use], on=index_cols, how=how,
                              suffixes=('', f'_{filenames[i]}'))

    # Rename columns with filename suffixes
    if value_cols:
        new_cols = index_cols.copy()
        for val_col in value_cols:
            for fname in filenames:
                new_cols.append(f"{val_col}-{fname}")
        if len(result.columns) == len(new_cols):
            result.columns = new_cols

    return result


def read_excel_or_csv(filepath: str) -> pd.DataFrame:
    """Read Excel or CSV file based on extension."""
    ext = os.path.splitext(filepath)[1].lower()

    if ext == '.xlsx' or ext == '.xls':
        return pd.read_excel(filepath)
    elif ext == '.csv':
        return pd.read_csv(filepath)
    else:
        raise ValueError('Please enter a xlsx or csv file.')
