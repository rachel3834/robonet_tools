from os import path, remove
from astropy.io import fits
import argparse
import numpy as np

def convert_spectrum(args):

    spectrum, dateobs = load_nres_fits_spectrum(args)

def load_nres_fits_spectrum(args):
    if not path.isfile(args.fits_file):
        raise IOError('Cannot find input FITS file')

    with fits.open(args.fits_file) as hdul:
        header = hdul[0].header
        data = hdul[0].data

    bunit = 1e-20
    dateobs = header['DATE-OBS']
    dwave = (header['XMAX'] - header['XMIN']) / header['NAXIS1']

    spectrum = np.zeros((header['NAXIS1'],2))
    spectrum[:,0] = np.arange(0,header['NAXIS1'],1) * dwave + header['XMIN']
    spectrum[:,1] = data[0,0,:] * bunit

    return spectrum, dateobs

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('fits_file', help='Path to input FITS extracted spectrum file')
    parser.add_argument('csv_file', help='Path to output CSV-format spectrum file')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    convert_spectrum(args)