import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config.lasers import color_for_laser


def plot_spectra_by_each_laser(spectra: pd.DataFrame, detectors: pd.DataFrame,
                               laser: str):
    """Produce a spill plot overlayed with instrument configuration for one laser."""
    excitation = detectors[detectors['laser'] == laser]['excitation'].iloc[0]
    laser_color = color_for_laser.get(laser, 'Black')

    xmin, xmax = 300, 900

    fig, ax = plt.subplots(figsize=(12, 3))

    # Plot laser excitation line
    ax.axvspan(excitation - 3, excitation + 3, alpha=0.8, color=laser_color)

    # Plot detector ranges
    detector_subset = detectors[detectors['laser'] == laser]
    for _, row in detector_subset.iterrows():
        ax.axvspan(row['xmin'], row['xmax'], alpha=0.3, color='#A3A3A3')

    # Plot spectra
    spectra_subset = spectra[spectra['laser'] == laser]
    fluorophores = spectra_subset['fluorophore'].unique()

    for fluor in fluorophores:
        fluor_data = spectra_subset[spectra_subset['fluorophore'] == fluor]

        for spectrum_type in fluor_data['spectrum_type'].unique():
            type_data = fluor_data[fluor_data['spectrum_type'] == spectrum_type]

            linestyle = ':' if spectrum_type in ['AB', 'EX'] else '-'
            alpha = 0.3

            ax.fill_between(type_data['Wavelength'], type_data['intensity'],
                            alpha=alpha, label=f"{fluor} {spectrum_type}")
            ax.plot(type_data['Wavelength'], type_data['intensity'],
                    color='black', alpha=0.7, linestyle=linestyle, linewidth=0.8)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(0, 1)
    ax.set_ylabel(laser)
    ax.legend(loc='upper right', fontsize=8)

    return fig
