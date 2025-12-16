import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
from astropy.coordinates import SkyCoord

os.environ['PYSYN_CDBS'] =  '/Users/rstreet/cdbs'

from astropy import units as u
from astropy.coordinates import SkyCoord
from Spyctres import Spyctres


def run(args):

    # Load spectrum data
    raw_data = np.loadtxt(args.data_file, skiprows=3)
    spectrum = np.c_[raw_data[:,0], raw_data[:,1], raw_data[:,1]*0.01]

    # Load model data if any
    if 'none' not in str(args.model_file).lower():
        model_data = np.loadtxt(args.model_file, skiprows=1)
    else:
        model_data = np.array([])

    # Load telluric lines
    telluric_lines, telluric_mask = Spyctres.load_telluric_lines(0.90)

    momo = Spyctres.star_spectrum_new(5000, 0, 4.5, catalog='k93models')
    wave = momo._model.points[0]
    mask = (wave < 10000) & (wave > 4000)
    wave_ref = wave[mask]

    # Plot
    fig, ax = plt.subplots(1,1, figsize=(10,8))
    #ax.errorbar(spectrum[:,0], spectrum[:,1], spectrum[:,2], fmt='.', ls='-', label='Data')
    ax.plot(spectrum[:,0], spectrum[:,1], ls='-', label='Data')
    if len(model_data) > 0:
        ax.plot(model_data[:,0], model_data[:,1], ls='-.', c='k', label='Model')
    ax.fill_between(wave_ref,0,10000,where=telluric_mask(wave_ref),color='grey',alpha=0.25, label='Tellurics')
    xmin = spectrum[:,0].min() * 0.99
    xmax = spectrum[:,0].max() * 1.01
    ymin = 100
    ymax = 1000
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.set_yscale('log')
    ax.set_xlabel(r'$\lambda [\AA]$', fontsize=18)
    ax.set_ylabel(r'$F_\lambda [erg/s/cm^2/\AA]$', fontsize=18)

    xticks = ax.get_xticks()
    xticklabels = ax.get_xticklabels()
    ax.set_xticks(xticks, labels=xticklabels, fontsize=16)
    yticks = ax.get_yticks()
    yticklabels = ax.get_yticklabels()
    ax.set_yticks(yticks, labels=yticklabels, fontsize=16)
    ax.legend(loc='lower right')
    ax.grid()
    plt.tight_layout()
    plt.savefig(args.output_file)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('data_file', help='Path to input spectrum file')
    parser.add_argument('model_file', help='Path to input spectrum file or None')
    parser.add_argument('output_file', help='Path to output spectrum plot file')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    run(args)