#!/usr/bin/env python3
"""
Graphs percentage of live cells for each cell population.
Expects the first four gates to be Cells/Single Cells/Single Cells/Live Cells.
If abs_count is included as a column in your metadata file, this script
will also calculate absolute counts.
"""

import os
import sys
import argparse
import json
import logging
from datetime import datetime
from tqdm import tqdm

import pandas as pd
import numpy as np

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from functions.preprocessing import sort_groups_by_metric
from functions.file_io import import_flowjo_export, import_flow_metadata
from tools.file_io import append_many_csv
from tools.list_tools import items_in_a_not_b, move_list_items_to_front
from tools.math import apply_unpaired_t_test, apply_multiple_comparisons
from tools.plotting import save_fig, plot_violin, plot_multiple_comparisons
from config.flow import mouse_db_ignore
from config.user_input import custom_group_order


def main():
    parser = argparse.ArgumentParser(description='Analyze flow cytometry counts')
    parser.add_argument('-i', '--input-dir', default='data/flow-counts',
                        help='Input directory with CSV files')
    parser.add_argument('-o', '--output-dir', default='hfc-output',
                        help='Output directory')
    parser.add_argument('-m', '--metadata-dir', default='data/flow-metadata',
                        help='Metadata directory')
    parser.add_argument('-r', '--ref-dir', default='data/mice',
                        help='Mouse reference data directory')
    parser.add_argument('-g', '--group-by', default='sex,treatment,zygosity',
                        help='Grouping columns (comma-separated)')
    parser.add_argument('-s', '--stat', default='fishers_lsd',
                        choices=['fishers_lsd', 't_test', 'tukey', 'bonferroni'],
                        help='Statistical test')
    parser.add_argument('-n', '--show-numbers', action='store_true',
                        help='Show p-values instead of stars')
    parser.add_argument('-l', '--height', type=int, default=500)
    parser.add_argument('-w', '--width', type=int, default=2000)
    parser.add_argument('-t', '--troubleshooting', action='store_true')

    args = parser.parse_args()

    # Setup
    wd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(wd, args.output_dir)
    metadata_cols = args.group_by.split(',')

    # Logging
    start_time = datetime.now()
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(message)s')
    logging.info(f'Script started at: {start_time}')

    # Read data
    logging.info('Reading data...')
    df = import_flowjo_export(
        os.path.join(wd, args.input_dir),
        metric_name='num_cells', include_initial_gates=True
    )

    if df is None or len(df) == 0:
        raise ValueError(f"No data found in {os.path.join(wd, args.input_dir)}")

    # Calculate percent cells
    live_col = 'Cells/Single Cells/Single Cells/Live Cells'
    if live_col in df.columns:
        df['pct_cells'] = (df['num_cells'] / df[live_col] * 100).round(4)

    # Load metadata
    flow_metadata = import_flow_metadata(os.path.join(wd, args.metadata_dir))
    mouse_db = append_many_csv(os.path.join(wd, args.ref_dir),
                               recursive=True, include_filepath=False)

    # Preprocessing
    logging.info('Preprocessing...')
    df = df.merge(flow_metadata, on='fcs_name', how='inner', suffixes=('', '_'))

    if 'is_unstained' in df.columns:
        df = df[df['is_unstained'] == False]

    # Create group name
    if len(metadata_cols) > 1:
        df['group_name'] = df[metadata_cols].apply(lambda x: ', '.join(x.astype(str)), axis=1)
    else:
        df['group_name'] = df[metadata_cols[0]]

    # Calculate absolute counts
    if 'total_viable_cells' in df.columns:
        df['abs_count'] = (df['total_viable_cells'] * df['pct_cells'] / 100).round(0)

    # Merge mouse data
    if mouse_db is not None and len(mouse_db) > 0:
        logging.info('Merging mouse data...')
        mouse_cols = items_in_a_not_b(list(mouse_db.columns), mouse_db_ignore)
        df = df.merge(mouse_db[mouse_cols], on='mouse_id', how='left', suffixes=('', '_'))
        if 'age' in df.columns:
            df['weeks_old'] = (df['age'] / 7).round(1)

    # Sort and filter
    df = df.sort_values(['cell_type', 'organ', 'group_name'])
    group_names = sorted(df['group_name'].unique())
    df = df[df['num_cells'] > 10]

    # Compute statistics
    logging.info(f'Computing {args.stat}...')
    if args.stat == 't_test':
        pval_tbl = apply_unpaired_t_test(
            df, index_cols=['organ', 'cell_type'],
            group_name='group_name', metric='pct_cells'
        )
    else:
        pval_tbl = apply_multiple_comparisons(
            df, index_cols=['organ', 'cell_type'],
            group_name='group_name', metric='pct_cells',
            correction=args.stat
        )

    # Save p-values
    if not args.troubleshooting and pval_tbl is not None:
        dirpath = os.path.join(output_dir, 'data',
                               args.group_by.replace(',', '_'), 'pct_cells')
        os.makedirs(dirpath, exist_ok=True)

        pval_tbl.to_csv(os.path.join(dirpath, f'pct_cells-{args.stat}.csv'), index=False)

    # Process each organ
    organs = sorted(df['organ'].unique())
    logging.info(f'Organs found: {", ".join(organs)}')

    for organ in organs:
        logging.info(f'Processing {organ}...')

        # Get data for this organ
        organ_df = df[df['organ'] == organ].copy()

        # Sort by metric
        organ_df = sort_groups_by_metric(
            organ_df, x='cell_type', y='pct_cells', groups=['group_name']
        )

        if not args.troubleshooting:
            # Save data
            dirpath = os.path.join(output_dir, 'data',
                                   args.group_by.replace(',', '_'), 'pct_cells')
            organ_df.to_csv(os.path.join(dirpath, f'{organ}.csv'), index=False)

            # Plot percent cells
            fig = plot_violin(
                organ_df, x='cell_type', y='pct_cells', group_by='group_name',
                ylabel='Percent of Live Cells', title=organ,
                ymin=0, ymax=100
            )
            save_fig(
                fig, height=args.height, width=args.width,
                dirpath=os.path.join(output_dir, 'figures',
                                     args.group_by.replace(',', '_'),
                                     'pct_cells', 'interactive'),
                filename=f'{organ}-pct_cells',
                save_html=True
            )

    # Multiple comparisons plots
    if pval_tbl is not None:
        logging.info(f'Plotting multiple comparisons using {args.stat}')

        for idx in tqdm(range(len(pval_tbl)), desc='Multiple comparisons'):
            organ = pval_tbl.iloc[idx]['organ']
            cell_type = pval_tbl.iloc[idx]['cell_type']

            df_subset = df[(df['organ'] == organ) & (df['cell_type'] == cell_type)]

            if len(df_subset) == 0:
                continue

            # Reorder groups
            local_group_order = move_list_items_to_front(
                df_subset['group_name'].unique().tolist(),
                custom_group_order
            )

            try:
                fig = plot_multiple_comparisons(
                    df_subset[['organ', 'cell_type', 'group_name', 'pct_cells']],
                    x='group_name', y='pct_cells',
                    ylabel='Percent of\nLive Cells',
                    title=f'{organ.upper()} {cell_type}',
                    test=args.stat,
                    show_numbers=args.show_numbers,
                    custom_group_order=local_group_order
                )

                if not args.troubleshooting:
                    dirpath = os.path.join(output_dir, 'figures',
                                           args.group_by.replace(',', '_'),
                                           'pct_cells', organ, args.stat)
                    os.makedirs(dirpath, exist_ok=True)

                    filename = f"{organ.lower().replace(' ', '')}-{cell_type.lower().replace(' ', '_')}.svg"
                    fig.savefig(os.path.join(dirpath, filename),
                                dpi=500, bbox_inches='tight')
                    plt.close(fig)
            except Exception as e:
                logging.warning(f"Error plotting {organ} {cell_type}: {e}")

    end_time = datetime.now()
    logging.info(f'Script ended at: {end_time}')
    logging.info(f'Script completed in: {end_time - start_time}')


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    main()
