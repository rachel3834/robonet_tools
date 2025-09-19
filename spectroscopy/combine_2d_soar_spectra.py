from os import path
import argparse
from astropy.io import fits
import copy
import numpy as np
import datetime

def coadd_spectra(args):

    # Sanity check that these are SOAR spectra, and that they have had
    # the instrument signatures removed etc.
    sanity_checks(args)

    # Build the data lists
    image_paths = [args.spectrum1, args.spectrum2]
    if str(args.spectrum3).lower() != 'none':
        image_paths.append(args.spectrum3)
    image_files = [path.basename(f) for f in image_paths]

    # Load the 2D spectra images
    images = load_2d_spectra(image_paths)

    # Build data array and coadd the data
    image_data = [images[f]['data'] for f in image_files]

    coadd_data = np.array(image_data).sum(axis=0)

    # Combine header keyword data
    coadd_hdr = combine_header_data(args, image_files, images)

    # Output combined spectrum
    output_coadded_spectrum(args, coadd_hdr, coadd_data, image_files, images)

def output_coadded_spectrum(args, coadd_hdr, coadd_data, image_files, images):

    hdu0 = fits.PrimaryHDU(data=coadd_data, header=coadd_hdr)
    hdulist = [hdu0]

    if 'mask_header' in images[image_files[0]].keys():
        mask_hdr = images[image_files[0]]['mask_header']
        mask_data = images[image_files[0]]['mask_data']
        hdu1 = fits.ImageHDU(data=mask_data, header=mask_hdr)
        hdulist.append(hdu1)

    if 'uncert_header' in images[image_files[0]].keys():
        uncert_hdr = images[image_files[0]]['uncert_header']
        uncert_data = np.sqrt(images[image_files[0]]['uncert_data'] ** 2 + images[image_files[1]]['uncert_data'] ** 2)
        hdu2 = fits.ImageHDU(data=uncert_data, header=uncert_hdr)
        hdulist.append(hdu2)

    new_hdul = fits.HDUList(hdulist)
    new_hdul.writeto(args.output_file, overwrite=True)

def combine_header_data(args, image_files, images):

    # Most of the header is a copy of the header of the first spectrum.
    # Specific keywords will be updated below.
    hdr1 = images[image_files[0]]['header']
    hdr2 = images[image_files[1]]['header']
    if len(images) == 3:
        hdr3 = images[image_files[1]]['header']

    hdr = copy.deepcopy(hdr1)

    # Update the timestamps.
    # Due to the formatting of the SOAR headers, the DATE-OBS stamp represents the start of the first exposure
    t1 = datetime.datetime.strptime(hdr1['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
    t2 = datetime.datetime.strptime(hdr2['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
    if len(images) == 3:
        t3 = datetime.datetime.strptime(hdr3['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')

    start = min(t1, t2)
    if len(images) == 3:
        start = min(start, t3)
    hdr['DATE-OBS'] = start.strftime('%Y-%m-%dT%H:%M:%S.%f')

    # The TIME keyword is a string in the format: "%H:%M:%S.%f to %H:%M:%S.%f" representing the start and end time of the exposures.
    # This is approximate.  Here we take the start from the first exposure and the end from the second, assuming they were back-to-back
    t1 = datetime.datetime.strptime(str(hdr1['TIME']).split()[0], '%H:%M:%S.%f')
    t2 = datetime.datetime.strptime(str(hdr2['TIME']).split()[0], '%H:%M:%S.%f')
    if len(images) == 3:
        t3 = datetime.datetime.strptime(str(hdr3['TIME']).split()[0], '%H:%M:%S.%f')
    start = min(t1, t2)
    if len(images) == 3:
        start = min(start, t3)

    t1 = datetime.datetime.strptime(str(hdr1['TIME']).split()[-1], '%H:%M:%S.%f')
    t2 = datetime.datetime.strptime(str(hdr2['TIME']).split()[-1], '%H:%M:%S.%f')
    if len(images) == 3:
        t3 = datetime.datetime.strptime(str(hdr3['TIME']).split()[-1], '%H:%M:%S.%f')
    end = max(t1, t2)
    if len(images) == 3:
        end = max(start, t3)

    hdr['TIME'] = start.strftime('%H:%M:%S.%f') + ' to ' + end.strftime('%H:%M:%S.%f')

    # Similarly, update the UT timestamp
    t1 = datetime.datetime.strptime(str(hdr1['UT']).split()[0], '%H:%M:%S.%f')
    t2 = datetime.datetime.strptime(str(hdr2['UT']).split()[0], '%H:%M:%S.%f')
    if len(images) == 3:
        t3 = datetime.datetime.strptime(str(hdr2['UT']).split()[0], '%H:%M:%S.%f')
    ut = min(t1, t2)
    if len(images) == 3:
        ut = min(start, t3)
    hdr['UT'] = ut.strftime('%H:%M:%S.%f')

    # Update filename stored in header
    hdr['GSP_FNAM'] = path.basename(args.output_file)

    # Update the exposure time
    hdr['EXPTIME'] = hdr1['EXPTIME'] + hdr2['EXPTIME']
    if len(images) == 3:
        hdr['EXPTIME'] += hdr3['EXPTIME']

    # Correct the readnoise for multiple exposures
    hdr['RDNOISE'] = hdr1['RDNOISE'] + hdr2['RDNOISE']
    if len(images) == 3:
        hdr['RDNOISE'] += hdr3['RDNOISE']

    return hdr

def load_2d_spectra(image_paths, debug=False):
    images = {}
    for fpath in image_paths:
        with fits.open(fpath) as hdul:
            if debug: print(hdul)

            entry = {
                'header': hdul[0].header,
                'data': hdul[0].data
            }
            if len(hdul) > 1:
                entry['mask_header'] = hdul[1].header
                entry['mask_data'] = hdul[1].data
                entry['uncert_header'] = hdul[2].header
                entry['uncert_data'] = hdul[2].data

            images[path.basename(fpath)] = entry

        hdul.close()

    return images

def sanity_checks(args):

    # Check that the data are available
    if not path.isfile(args.spectrum1) or not path.isfile(args.spectrum2):
        raise IOError('Cannot find one or more of the first two spectrum files')

    if str(args.spectrum3).lower() != 'none' and not path.isfile(args.spectrum3):
        raise IOError('Cannot find the third spectrum file')

    # Check that the spectra provided are SOAR spectra and have the instrument signature
    # removed
    file_list = [args.spectrum1, args.spectrum2]
    if str(args.spectrum3).lower() != 'none':
        file_list.append(args.spectrum3)

    for f in file_list:
        hdr = fits.getheader(f)
        if hdr['TELESCOP'] != 'SOAR 4.1m':
            raise IOError('Input spectrum ' + path.basename(f) + ' is not a SOAR spectrum')

        if 'cfzst_' not in path.basename(f) and 'cfsto_' not in path.basename(f):
            raise IOError('Input spectrum ' + path.basename(f)
                          + ' is not at the expected stage of reduction')
def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('spectrum1', help='Path to first 2D spectrum file')
    parser.add_argument('spectrum2', help='Path to second 2D spectrum file')
    parser.add_argument('spectrum3', help='Path to third 2D spectrum file or None')
    parser.add_argument('output_file', help='Path to output combined spectrum')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    coadd_spectra(args)