#!/usr/bin/env python3
"""
Monitor survival and weight tracking.

Generates:
1. Kaplan-Meier survival curve with confidence intervals
2. Percent weight over time and raw weight over time
3. Spleen weight comparison if provided
"""

import os
import sys
import argparse
import logging
from datetime import datetime

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from functions.preprocessing import create_survival_table
from tools.plotting import plot_scatter, plot_multiple_comparisons, save_fig
from tools.file_io import read_excel_or_csv
from tools.df_tools import append_dataframe, fillna, rename_columns


def main():
    parser = argparse.ArgumentParser(description='Monitor weights and survival')
    parser.add_argument('-i', '--input-file', default='data/weights.xlsx',
                        help='Input weights file')
    parser.add_argument('-o', '--output-dir', default='figures/survival',
                        help='Output directory')
    parser.add_argument('-g', '--group-by', default='treatment,genotype',
                        help='Grouping columns')
    parser.add_argument('-w', '--week', action='store_true',
                        help='Use weeks instead of days')
    parser.add_argument('-p', '--optional-plots', action='store_true',
                        help='Generate optional plots')
    parser.add_argument('-t', '--troubleshooting', action='store_true')

    args = parser.parse_args()

    # Setup
    wd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    metadata_cols = args.group_by.split(',')

    # Logging
    start_time = datetime.now()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
    logging.info(f'Script started at: {start_time}')

    # Read file
    logging.info('Reading data...')
    input_file = os.path.join(wd, args.input_file)
    weight_tbl = read_excel_or_csv(input_file)

    ext = os.path.splitext(input_file)[1].lower()
    if ext == '.xlsx':
        date_pattern = r'^\d+$'
    else:
        date_pattern = r'^\d{1,2}/\d{1,2}/\d{2,4}$'

    # Identify date columns
    is_date = weight_tbl.columns.str.match(date_pattern)
    id_cols = weight_tbl.columns[~is_date].tolist()
    date_cols = weight_tbl.columns[is_date].tolist()

    weight_tbl = weight_tbl[~weight_tbl['mouse_id'].isna()]

    # Derive group column
    metadata_cols = [c for c in metadata_cols if c in weight_tbl.columns]
    if len(metadata_cols) > 1:
        weight_tbl['group'] = weight_tbl[metadata_cols].apply(
            lambda x: ', '.join(x.astype(str)), axis=1
        )
    else:
        weight_tbl['group'] = weight_tbl[metadata_cols[0]]

    # Pivot to long format
    df = weight_tbl.melt(
        id_vars=[c for c in weight_tbl.columns if c not in date_cols],
        value_vars=date_cols,
        var_name='date',
        value_name='weight'
    )

    # Fix dates for Excel
    if ext == '.xlsx':
        df['date'] = pd.to_datetime(df['date'].astype(float), unit='D',
                                     origin='1899-12-30')
    else:
        df['date'] = pd.to_datetime(df['date'])

    # Plot survival
    logging.info('Plotting survival...')

    # Extract start and end info
    start_info = df.groupby('mouse_id').apply(
        lambda x: x.loc[x['date'].idxmin()]
    ).reset_index(drop=True)
    start_info = start_info[['mouse_id', 'group', 'date', 'weight']].copy()
    if 'spleen_weight' in df.columns:
        start_info = start_info.merge(
            weight_tbl[['mouse_id', 'spleen_weight']], on='mouse_id', how='left'
        )

    end_info = df[~df['weight'].isna()].groupby('mouse_id').apply(
        lambda x: x.loc[x['date'].idxmax()]
    ).reset_index(drop=True)
    end_info = end_info[['mouse_id', 'group', 'date', 'weight']].copy()

    survival = create_survival_table(start_info, end_info, end_date=df['date'].max())

    # Kaplan-Meier plot
    fig, ax = plt.subplots(figsize=(10, 6))

    time_col = 'week' if args.week else 'day'
    for group in survival['group'].unique():
        group_data = survival[survival['group'] == group]

        kmf = KaplanMeierFitter()
        kmf.fit(group_data[time_col], group_data['is_dead'], label=group)
        kmf.plot_survival_function(ax=ax, ci_show=True)

    ax.set_xlabel(f'Time ({"weeks" if args.week else "days"})')
    ax.set_ylabel('Survival Probability')
    ax.set_title('Survival Curve')
    ax.set_xlim(0, None)
    ax.set_ylim(0, 1)
    ax.legend(title='Group')

    if not args.troubleshooting:
        dirpath = os.path.join(wd, args.output_dir)
        os.makedirs(dirpath, exist_ok=True)
        fig.savefig(os.path.join(dirpath, 'survival-curve.png'),
                    dpi=300, bbox_inches='tight')
        plt.close(fig)

    # Plot weights
    logging.info('Plotting weights...')

    if not args.troubleshooting:
        dirpath = os.path.join(wd, args.output_dir, 'mouse-weights')
        os.makedirs(dirpath, exist_ok=True)

    # Merge start info for calculations
    df = df.merge(
        start_info[['mouse_id', 'date', 'weight']].rename(
            columns={'date': 'date_start', 'weight': 'weight_start'}
        ),
        on='mouse_id', how='left'
    )

    df['day'] = (df['date'] - df['date_start']).dt.days
    df['week'] = df['day'] / 7
    df['pct_weight'] = df['weight'] / df['weight_start']

    # Percent weight plot
    color_map = {}
    for group in df['group'].unique():
        if group in ['WT', 'H']:
            color_map[group] = '#FF7F0E'
        else:
            color_map[group] = 'rgba(23,190,207,1)'

    fig = plot_scatter(
        df,
        x='week' if args.week else 'day',
        y='pct_weight', group_by='group',
        ymin=0.75,
        xlabel=f'Time ({"weeks" if args.week else "days"})',
        ylabel='Percent Weight',
        title='Percent Weight Over Time',
        mode='lines+markers',
        color_discrete_map=color_map
    )

    if not args.troubleshooting:
        save_fig(fig, height=400, width=750, dirpath=dirpath,
                 filename='pct_weight', save_html=True)

    # Raw weight plot
    fig = plot_scatter(
        df,
        x='week' if args.week else 'day',
        y='weight', group_by='group',
        xlabel=f'Time ({"weeks" if args.week else "days"})',
        ylabel='Mouse Weight (g)',
        title='Raw Weight Over Time',
        mode='lines+markers',
        color_discrete_map=color_map
    )

    if not args.troubleshooting:
        save_fig(fig, height=400, width=600, dirpath=dirpath,
                 filename='raw_weight', save_html=True)

    # Weight loss at death
    survival['weight_ratio'] = survival['weight_end'] / survival['weight_start']

    try:
        fig = plot_multiple_comparisons(
            survival, x='group', y='weight_ratio',
            ymin=0.70,
            ylabel='[Final weight] / [Starting weight]',
            title='Weight Loss at Time of Death',
            show_numbers=True,
            test='fishers_lsd',
            custom_group_order=['WT', 'KO']
        )

        if not args.troubleshooting:
            fig.savefig(os.path.join(wd, args.output_dir, 'weight-loss.png'),
                        dpi=300, bbox_inches='tight')
            plt.close(fig)
    except Exception as e:
        logging.warning(f"Could not plot weight loss: {e}")

    # Spleen weight
    if 'spleen_weight' in weight_tbl.columns:
        logging.info('Plotting spleen weight...')

        if not args.troubleshooting:
            dirpath = os.path.join(wd, args.output_dir, 'spleen-weight')
            os.makedirs(dirpath, exist_ok=True)

        try:
            fig = plot_multiple_comparisons(
                survival, x='group', y='spleen_weight',
                ylabel='Spleen Weight (mg)',
                title='Spleen Weight',
                show_numbers=False,
                test='fishers_lsd',
                custom_group_order=['WT', 'KO']
            )

            if not args.troubleshooting:
                fig.savefig(os.path.join(dirpath, 'spleen_weight.svg'),
                            bbox_inches='tight')
                plt.close(fig)
        except Exception as e:
            logging.warning(f"Could not plot spleen weight: {e}")

    end_time = datetime.now()
    logging.info(f'Script ended at: {end_time}')
    logging.info(f'Script completed in: {end_time - start_time}')


if __name__ == '__main__':
    main()
