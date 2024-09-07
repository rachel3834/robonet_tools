from os import path, remove
from astropy.io import fits
from astropy.modeling.polynomial import Chebyshev1D, Polynomial1D
import argparse
import numpy as np

def convert_spectrum(args):

    spectrum, dateobs = load_goodman_fits_spectrum(args)

    output_csv_spectrum(args, spectrum, dateobs)

def output_csv_spectrum(args, spectrum, dateobs):

    if path.isfile(args.csv_file):
        remove(args.csv_file)

    with open(args.csv_file, 'w') as f:
        f.write('# FACILITY: SOAR\n')
        f.write('# DATE-OBS: ' + dateobs.split('T')[0] + '\n')
        f.write('wavelength flux\n')

        for i in range(0,spectrum.shape[0],1):
            f.write(str(spectrum[i,0]) + ' ' + str(spectrum[i,1]) + '\n')

        f.close()


def load_goodman_fits_spectrum(args):
    if not path.isfile(args.fits_file):
        raise IOError('Cannot find input FITS file')

    with fits.open(args.fits_file) as hdul:
        header = hdul[0].header
        data = hdul[0].data

    dateobs = header['DATE-OBS']

    # Coefficients of the model used for the wavelength solution
    model = header['GSP_FUNC']
    order = header['GSP_ORDR']
    degree = 4
    C = [header['GSP_C00'+str(i)] for i in range(0, order+1, 1)]
    T = [header['GSP_TC0'+str(i)] for i in range(0, 3, 1)]

    spectrum = np.zeros((header['NAXIS1'], 2))
    spectrum[:,0] = np.arange(1, header['NAXIS1']+1, 1)
    wmodel = Polynomial1D(header['GSP_TORD'])
    fmodel = Chebyshev1D(degree)
    norm = fmodel.evaluate(data, *np.array(T))
    spectrum[:,0] = wmodel.evaluate(spectrum[:,0], *np.array(C))
    spectrum[:,1] = data

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