import pandas as pd
from typing import Optional
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tools.file_io import append_many_csv
from tools.df_tools import reset_index
from functions.preprocessing import preprocess_flowjo_export


def import_flowjo_export(dirpath: str, metric_name: str = 'num_cells',
                         include_initial_gates: bool = False) -> Optional[pd.DataFrame]:
    """Import FlowJo export data (counts, gmfi, or sdev)."""
    raw_tables = append_many_csv(dirpath, recursive=True,
                                  na_strings=['n/a'], include_filepath=False,
                                  return_list=True)
    if not raw_tables:
        return None

    dfs = []
    for key, table in raw_tables.items():
        processed = preprocess_flowjo_export(
            table, metric_name=metric_name,
            include_initial_gates=include_initial_gates
        )
        dfs.append(processed)

    if not dfs:
        return None

    df = pd.concat(dfs, ignore_index=True)
    return df


def import_flow_metadata(dirpath: str) -> pd.DataFrame:
    """Import flow metadata from directory."""
    flow_metadata = append_many_csv(dirpath, recursive=True, include_filepath=False)

    if flow_metadata is None or len(flow_metadata) == 0:
        raise ValueError(f"No metadata found. Please check {dirpath}...")

    # Handle 'sex' column interpreted as boolean
    if 'sex' in flow_metadata.columns:
        if flow_metadata['sex'].dtype == bool:
            flow_metadata['sex'] = flow_metadata['sex'].map({False: 'F', True: 'T'})

    return flow_metadata
