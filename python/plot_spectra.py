#!/usr/bin/env python3
"""Plot spectra for fluorophores used in flow panel."""

import os
import sys
import argparse
import logging
from datetime import datetime

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from tools.list_tools import multiple_replacement, items_in_a_not_b
from tools.file_io import read_excel_or_csv
from tools.df_tools import reset_index
from tools.text_tools import substr_right
from config.lasers import lasers
from config.replacements import fluorophore_replacements
from functions.preprocessing import preprocess_instrument_config
from functions.plotting import plot_spectra_by_each_laser


def main():
    parser = argparse.ArgumentParser(description='Plot fluorophore spectra')
    parser.add_argument('-i', '--input-file', default='data/design/panel.csv',
                        help='Panel input file')
    parser.add_argument('-a', '--antibody-inventory', default='ref/antibody_inventory.xlsx',
                        help='Antibody inventory')
    parser.add_argument('-c', '--instrument-config', default='ref/instrument_config.xlsx',
                        help='Instrument configuration')
    parser.add_argument('-s', '--spectra-file', default='ref/spectra.csv',
                        help='Spectra data file')
    parser.add_argument('-f', '--figures-dir', default='figures/spectra',
                        help='Output directory for figures')
    parser.add_argument('-t', '--troubleshooting', action='store_true')

    args = parser.parse_args()

    # Setup
    wd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    figures_dir = os.path.join(wd, args.figures_dir)

    # Logging
    start_time = datetime.now()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
    logging.info(f'Script started at: {start_time}')

    # Load instrument config
    instr_cfg = read_excel_or_csv(os.path.join(wd, args.instrument_config))
    instr_cfg = preprocess_instrument_config(instr_cfg)

    # Expand fluorophores
    instr_cfg_long = instr_cfg.copy()
    if 'fluorophore' in instr_cfg_long.columns:
        instr_cfg_long = instr_cfg_long.assign(
            fluorophore=instr_cfg_long['fluorophore'].str.split(', ')
        ).explode('fluorophore')

    # Parse bandpass filter
    if 'bandpass_filter' in instr_cfg.columns:
        bp_split = instr_cfg['bandpass_filter'].str.split('/', expand=True)
        instr_cfg['emission'] = bp_split[0].astype(float)
        instr_cfg['range'] = bp_split[1].astype(float)
        instr_cfg['xmin'] = instr_cfg['emission'] - instr_cfg['range'] / 2
        instr_cfg['xmax'] = instr_cfg['emission'] + instr_cfg['range'] / 2

    # Load panel
    panel = read_excel_or_csv(os.path.join(wd, args.input_file))
    panel['fluorophore'] = multiple_replacement(
        panel['fluorophore'].fillna('').tolist(), fluorophore_replacements
    )
    panel = panel[panel['fluorophore'] != ''].merge(
        instr_cfg_long[['fluorophore', 'laser']],
        on='fluorophore', how='left'
    )
    fluorophores = panel['fluorophore'].tolist()

    # Load spectra
    all_spectra = read_excel_or_csv(os.path.join(wd, args.spectra_file))

    # Preprocess column names
    all_spectra.columns = multiple_replacement(
        list(all_spectra.columns), fluorophore_replacements
    )
    cols = all_spectra.columns.tolist()
    available_fluorophores = sorted(set(
        col.replace(' EM', '').replace(' EX', '').replace(' AB', '')
        for col in cols[1:]
    ))

    # Subset spectra
    pattern_cols = [col for col in cols
                    if any(f in col for f in fluorophores)]
    spectra = all_spectra[['Wavelength'] + pattern_cols]

    # Reshape to long format
    spectra_long = spectra.melt(
        id_vars='Wavelength',
        var_name='trace_name',
        value_name='intensity'
    ).dropna()

    spectra_long['fluorophore'] = spectra_long['trace_name'].str.replace(r' ..$', '', regex=True)
    spectra_long['spectrum_type'] = spectra_long['trace_name'].apply(lambda x: substr_right(x, 2))

    spectra_long = spectra_long.merge(panel, on='fluorophore', how='left')
    spectra_long = spectra_long.sort_values(['trace_name', 'Wavelength'])
    spectra_long = spectra_long[spectra_long['fluorophore'].isin(fluorophores)]

    # Find unavailable fluorophores
    unavailable = sorted(set(fluorophores) - set(available_fluorophores))
    if unavailable and not args.troubleshooting:
        troubleshooting_dir = os.path.join(figures_dir, 'troubleshooting')
        os.makedirs(troubleshooting_dir, exist_ok=True)
        with open(os.path.join(troubleshooting_dir, 'unplotted_fluorophores.txt'), 'w') as f:
            f.write('\n'.join(unavailable))

    # Plot
    available_lasers = [l for l in lasers if l in panel['laser'].unique()]

    fig = plt.figure(figsize=(12, 3 * len(available_lasers)))
    gs = GridSpec(len(available_lasers), 1, figure=fig)

    for i, laser in enumerate(available_lasers):
        ax = fig.add_subplot(gs[i])
        laser_spectra = spectra_long[spectra_long['laser'] == laser]

        if 'excitation' in instr_cfg.columns:
            laser_detectors = instr_cfg[instr_cfg['laser'] == laser]
            if len(laser_detectors) > 0:
                excitation = laser_detectors['excitation'].iloc[0]
                ax.axvspan(excitation - 3, excitation + 3, alpha=0.8, color='gray')

                for _, row in laser_detectors.iterrows():
                    if 'xmin' in row and 'xmax' in row:
                        ax.axvspan(row['xmin'], row['xmax'], alpha=0.3, color='#A3A3A3')

        for fluor in laser_spectra['fluorophore'].unique():
            fluor_data = laser_spectra[laser_spectra['fluorophore'] == fluor]
            for st in fluor_data['spectrum_type'].unique():
                st_data = fluor_data[fluor_data['spectrum_type'] == st]
                ls = ':' if st in ['AB', 'EX'] else '-'
                ax.fill_between(st_data['Wavelength'], st_data['intensity'],
                                alpha=0.3, label=f'{fluor} {st}')
                ax.plot(st_data['Wavelength'], st_data['intensity'],
                        color='black', alpha=0.7, linestyle=ls, linewidth=0.8)

        ax.set_xlim(300, 900)
        ax.set_ylim(0, 1)
        ax.set_ylabel(laser)
        ax.legend(loc='upper right', fontsize=6)

    plt.tight_layout()

    if not args.troubleshooting:
        os.makedirs(figures_dir, exist_ok=True)
        fig.savefig(os.path.join(figures_dir, 'spectra.png'),
                    dpi=200, bbox_inches='tight')
        plt.close(fig)

    end_time = datetime.now()
    logging.info(f'Script ended at: {end_time}')
    logging.info(f'Script completed in: {end_time - start_time}')


if __name__ == '__main__':
    main()
