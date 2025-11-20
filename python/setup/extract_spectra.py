#!/usr/bin/env python3
"""Extract spectra from JSON data downloaded from Biolegend's Spectra Analyzer."""

import os
import sys
import json
import argparse
import logging
from datetime import datetime

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tools.file_io import join_many_csv


def main():
    parser = argparse.ArgumentParser(description='Extract spectra from JSON')
    parser.add_argument('-i', '--input-dir', default='data/raw_spectra',
                        help='Input directory with JSON files')
    parser.add_argument('-o', '--output-dir', default='data/raw_spectra/output',
                        help='Output directory')
    parser.add_argument('-t', '--troubleshooting', action='store_true')

    args = parser.parse_args()

    # Setup
    wd = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    # Logging
    start_time = datetime.now()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
    logging.info(f'Script started at: {start_time}')

    # Create output directory
    os.makedirs(os.path.join(wd, args.output_dir), exist_ok=True)

    # Get JSON files
    input_dir = os.path.join(wd, args.input_dir)
    filepaths = [os.path.join(input_dir, f) for f in os.listdir(input_dir)
                 if f.endswith('.json')]
    logging.info(f'Number of files found: {len(filepaths)}')

    # Convert each file
    for filepath in filepaths:
        logging.info(f'Processing: {os.path.basename(filepath)}')
        filename = os.path.splitext(os.path.basename(filepath))[0]

        with open(filepath) as f:
            raw_data = json.load(f)

        df = pd.DataFrame(raw_data['data'], columns=['Wavelength', 'Intensity'])

        # Normalize
        max_intensity = df['Intensity'].max()
        if max_intensity > 100:
            df['Intensity'] = df['Intensity'] / max_intensity
        elif max_intensity > 1:
            df['Intensity'] = df['Intensity'] / 100

        if not args.troubleshooting:
            df.to_csv(os.path.join(wd, args.output_dir, f'{filename}.csv'),
                      index=False)

    # Combine output
    logging.info('Combining...')
    output_files = [os.path.join(wd, args.output_dir, f)
                    for f in os.listdir(os.path.join(wd, args.output_dir))
                    if f.endswith('.csv')]

    if output_files:
        combined = join_many_csv(
            output_files, index_cols=['Wavelength'], value_cols=['Intensity'],
            all_x=True, all_y=True
        )

        # Remove half wavelengths
        combined = combined[combined['Wavelength'] % 1 < 0.01]

        if not args.troubleshooting:
            parent_dir = os.path.dirname(os.path.join(wd, args.output_dir))
            combined.to_csv(os.path.join(parent_dir, 'combined-spectra.csv'),
                            index=False)

    end_time = datetime.now()
    logging.info(f'Script ended at: {end_time}')
    logging.info(f'Script completed in: {end_time - start_time}')


if __name__ == '__main__':
    main()
