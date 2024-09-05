from os import path, remove
from astropy.io import fits
import argparse
import numpy as np

def convert_spectrum(args):

    spectrum, dateobs = load_floyds_fits_spectrum(args)

    output_csv_spectrum(args, spectrum, dateobs)

def output_csv_spectrum(args, spectrum, dateobs):

    if path.isfile(args.csv_file):
        remove(args.csv_file)

    with open(args.csv_file, 'w') as f:
        f.write('# FACILITY: LCO\n')
        f.write('# DATE-OBS: ' + dateobs.split('T')[0] + '\n')
        f.write('wavelength flux\n')

        for i in range(0,spectrum.shape[0],1):
            f.write(str(spectrum[i,0]) + ' ' + str(spectrum[i,1]) + '\n')

        f.close()


def load_floyds_fits_spectrum(args):
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