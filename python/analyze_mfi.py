#!/usr/bin/env python3
"""
Graphs MFI (Mean Fluorescence Intensity) for each cell population.
Supports geometric mean (recommended), median, or arithmetic mean.
"""

import os
import sys
import argparse
import logging
from datetime import datetime
from tqdm import tqdm

import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from functions.preprocessing import sort_groups_by_metric
from functions.file_io import import_flowjo_export, import_flow_metadata
from tools.file_io import append_many_csv
from tools.df_tools import group_by_agg, rename_columns
from tools.list_tools import items_in_a_not_b, move_list_items_to_front
from tools.math import (apply_unpaired_t_test, apply_multiple_comparisons,
                        generate_lognormal_data, compute_normal_tvd)
from tools.plotting import (save_fig, plot_violin, plot_multiple_comparisons,
                            plot_modal_histograms)
from config.flow import mouse_db_ignore
from config.user_input import custom_group_order


def main():
    parser = argparse.ArgumentParser(description='Analyze flow cytometry MFI')
    parser.add_argument('-i', '--input-dir', default='data/flow-mfi',
                        help='Input MFI directory')
    parser.add_argument('-o', '--output-dir', default='hfc-output',
                        help='Output directory')
    parser.add_argument('-d', '--sdev-dir', default='data/flow-sdev',
                        help='Standard deviations directory')
    parser.add_argument('-c', '--counts-dir', default='data/flow-counts',
                        help='Counts directory')
    parser.add_argument('-m', '--metadata-dir', default='data/flow-metadata',
                        help='Metadata directory')
    parser.add_argument('-r', '--ref-dir', default='data/mice',
                        help='Mouse reference directory')
    parser.add_argument('-g', '--group-by', default='sex,treatment,zygosity',
                        help='Grouping columns')
    parser.add_argument('-p', '--plot-histograms', action='store_true',
                        help='Generate histograms')
    parser.add_argument('-y', '--metric', default='gmfi',
                        help='Metric name (gmfi, mode)')
    parser.add_argument('-s', '--stat', default='fishers_lsd',
                        choices=['fishers_lsd', 't_test', 'tukey', 'bonferroni'])
    parser.add_argument('-n', '--show-numbers', action='store_true')
    parser.add_argument('-f', '--fluorescence', default='zygosity',
                        help='Fluorescence grouping column')
    parser.add_argument('-x', '--xlabel', default='Comp-GFP-A :: mNeonGreen',
                        help='X-axis label for histograms')
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
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
    logging.info(f'Script started at: {start_time}')

    # Read MFI data
    logging.info('Reading data...')
    df = import_flowjo_export(
        os.path.join(wd, args.input_dir),
        metric_name=args.metric, include_initial_gates=False
    )

    if df is None:
        raise ValueError(f'No files found in {os.path.join(wd, args.input_dir)}')

    # Join sdev
    sdev_df = import_flowjo_export(
        os.path.join(wd, args.sdev_dir), metric_name='sdev'
    )
    if sdev_df is not None:
        logging.info('Merging sdev data...')
        df = df.merge(sdev_df, on=['fcs_name', 'gate', 'cell_type'],
                      how='left', suffixes=('', '_'))
        df['sdev'] = df['sdev'].fillna(0)

    # Join counts
    counts_df = import_flowjo_export(
        os.path.join(wd, args.counts_dir),
        metric_name='num_cells', include_initial_gates=True
    )
    if counts_df is not None:
        logging.info('Merging counts data...')
        df = df.merge(counts_df[['fcs_name', 'gate', 'cell_type', 'num_cells']],
                      on=['fcs_name', 'gate', 'cell_type'],
                      how='left', suffixes=('', '_'))

    # Join metadata
    flow_metadata = import_flow_metadata(os.path.join(wd, args.metadata_dir))
    df = df.merge(flow_metadata, on='fcs_name', how='inner', suffixes=('', '_'))

    if 'is_unstained' in df.columns:
        df = df[df['is_unstained'] == False]

    # Create group names
    if len(metadata_cols) > 1:
        df['group_name'] = df[metadata_cols].apply(
            lambda x: ', '.join(x.astype(str)), axis=1
        )
        subgroup_cols = [c for c in metadata_cols if c != args.fluorescence]
        df['subgroup_name'] = df[subgroup_cols].apply(
            lambda x: ', '.join(x.astype(str)), axis=1
        )
    else:
        df['group_name'] = df[metadata_cols[0]]

    # Join mouse data
    mouse_db = append_many_csv(os.path.join(wd, args.ref_dir),
                               recursive=True, include_filepath=False)
    if mouse_db is not None:
        logging.info('Merging mouse data...')
        mouse_cols = items_in_a_not_b(list(mouse_db.columns), mouse_db_ignore)
        df = df.merge(mouse_db[mouse_cols], on='mouse_id',
                      how='left', suffixes=('', '_'))
        if 'age' in df.columns:
            df['weeks_old'] = (df['age'] / 7).round(1)

    # Sort
    df = df.sort_values(['cell_type', 'organ', 'group_name'])

    # Process each organ
    organs = sorted(df['organ'].unique())
    logging.info(f'Organs included: {", ".join(organs)}')

    for organ in organs:
        logging.info(f'Plotting {organ} overview...')

        organ_df = df[df['organ'] == organ].copy()
        organ_df = sort_groups_by_metric(
            organ_df, x='cell_type', y=args.metric, groups=['group_name']
        )

        if not args.troubleshooting:
            dirpath = os.path.join(output_dir, 'data',
                                   args.group_by.replace(',', '_'), args.metric)
            os.makedirs(dirpath, exist_ok=True)
            organ_df.to_csv(os.path.join(dirpath, f'{args.metric}-{organ}.csv'),
                            index=False)

            fig = plot_violin(
                organ_df, x='cell_type', y=args.metric, group_by='group_name',
                ylabel='MFI', title=organ, ymin=0
            )
            save_fig(
                fig, height=args.height, width=args.width,
                dirpath=os.path.join(output_dir, 'figures',
                                     args.group_by.replace(',', '_'),
                                     args.metric, 'interactive'),
                filename=f'{organ}-{args.metric}',
                save_html=True
            )

    # Compute statistics
    logging.info(f'Computing {args.stat}...')

    if args.stat == 't_test':
        pval_tbl = apply_unpaired_t_test(
            df, index_cols=['organ', 'cell_type'],
            group_name='group_name', metric=args.metric
        )
    else:
        pval_tbl = apply_multiple_comparisons(
            df, index_cols=['organ', 'cell_type'],
            group_name='group_name', metric=args.metric,
            correction=args.stat
        )

    if pval_tbl is not None:
        logging.info(f'{len(pval_tbl)} populations found across {len(organs)} organs')

        if not args.troubleshooting:
            dirpath = os.path.join(output_dir, 'data',
                                   args.group_by.replace(',', '_'), args.metric)
            pval_tbl.to_csv(os.path.join(dirpath, f'gmfi-pvals-{args.stat}.csv'),
                            index=False)

        # Plot multiple comparisons
        logging.info(f'Plotting multiple comparisons using {args.stat}')

        for idx in tqdm(range(len(pval_tbl)), desc='MFI comparisons'):
            organ = pval_tbl.iloc[idx]['organ']
            cell_type = pval_tbl.iloc[idx]['cell_type']

            df_subset = df[(df['organ'] == organ) & (df['cell_type'] == cell_type)]

            if len(df_subset) == 0:
                continue

            local_group_order = move_list_items_to_front(
                df_subset['group_name'].unique().tolist(),
                custom_group_order
            )

            try:
                fig = plot_multiple_comparisons(
                    df_subset[['organ', 'cell_type', 'group_name', args.metric]],
                    x='group_name', y=args.metric,
                    ylabel=args.metric,
                    title=f'{organ.upper()} {cell_type}',
                    test=args.stat,
                    show_numbers=args.show_numbers,
                    custom_group_order=local_group_order
                )

                if not args.troubleshooting:
                    dirpath = os.path.join(output_dir, 'figures',
                                           args.group_by.replace(',', '_'),
                                           args.metric, organ, args.stat)
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
