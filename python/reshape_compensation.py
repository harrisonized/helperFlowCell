#!/usr/bin/env python3
"""Reshape compensation matrix for FACSDiva import."""

import os
import sys
import argparse
import logging
from datetime import datetime

import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from functions.preprocessing import spillover_to_xml
from tools.df_tools import rename_columns, reset_index
from tools.list_tools import move_list_items_to_front


def main():
    parser = argparse.ArgumentParser(description='Reshape compensation matrix')
    parser.add_argument('-i', '--input', default='data/compensation/Acquisition-defined.csv',
                        help='Input spillover matrix CSV')
    parser.add_argument('-o', '--output-dir', default='data/compensation',
                        help='Output directory')
    parser.add_argument('-s', '--sort', default='ref/fluorochrome_order.txt',
                        help='Fluorochrome order file')
    parser.add_argument('-t', '--troubleshooting', action='store_true')

    args = parser.parse_args()

    # Setup
    wd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # Logging
    start_time = datetime.now()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
    logging.info(f'Script started at: {start_time}')

    # Create output directory
    os.makedirs(os.path.join(wd, args.output_dir), exist_ok=True)

    # Read spillover matrix
    spillover_matrix = pd.read_csv(
        os.path.join(wd, args.input),
        index_col=0
    )
    channels = spillover_matrix.index.tolist()

    # Sort if order file exists
    order_file = os.path.join(wd, args.sort)
    if os.path.exists(order_file):
        logging.info('Fluorochrome order file found')
        with open(order_file) as f:
            fluorochrome_order = [line.strip() for line in f]

        fluorochrome_order = [f for f in fluorochrome_order
                              if f in spillover_matrix.columns]
        fluorochrome_order = move_list_items_to_front(
            spillover_matrix.columns.tolist(), fluorochrome_order
        )
        spillover_matrix = spillover_matrix.loc[fluorochrome_order, fluorochrome_order]
        channels = spillover_matrix.index.tolist()

    # Pivot to long format
    spillover_matrix_reset = spillover_matrix.reset_index()
    spillover_matrix_reset.columns = ['- % Fluorochrome'] + channels

    spillover_table = spillover_matrix_reset.melt(
        id_vars='- % Fluorochrome',
        value_vars=channels,
        var_name='Fluorochrome',
        value_name='Spectral Overlap'
    )

    # Export spillover table
    if not args.troubleshooting:
        diva_spillover = spillover_table[
            spillover_table['Fluorochrome'] != spillover_table['- % Fluorochrome']
        ].copy()
        diva_spillover['Spectral Overlap'] = (
            diva_spillover['Spectral Overlap'] * 100
        ).round(2)

        filepath = os.path.join(wd, args.output_dir, 'spillover_table.csv')
        diva_spillover.to_csv(filepath, index=False)

    # Compute compensation matrix (inverse of spillover)
    compensation_matrix = np.linalg.inv(spillover_matrix.values)

    if not args.troubleshooting:
        # Flatten and save
        filepath = os.path.join(wd, args.output_dir, 'compensation_matrix.txt')
        with open(filepath, 'w') as f:
            f.write(','.join(str(round(v, 8)) for v in compensation_matrix.flatten()))

    # Export XML
    filename = os.path.splitext(os.path.basename(args.input))[0]
    xml_content = spillover_to_xml(spillover_table, channels, name=filename)

    if not args.troubleshooting:
        filepath = os.path.join(wd, args.output_dir, f'{filename}.mtx')
        with open(filepath, 'w') as f:
            f.write(xml_content)

    end_time = datetime.now()
    logging.info(f'Script ended at: {end_time}')
    logging.info(f'Script completed in: {end_time - start_time}')


if __name__ == '__main__':
    main()
