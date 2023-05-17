from os import path, makedirs, remove
import numpy as np
from pyDANDIA import hd5_utils
from pyDANDIA import crossmatch
from pyDANDIA import field_photometry
import crossmatch_to_ipactable
import aws_utils
import argparse
from datetime import datetime
from astropy.table import Table, Column
from astropy.io import fits

VERSION = 0.1
FILTER_LIST = ['gp', 'rp', 'ip']

def convert_to_ipac_lightcurves(args):

    # Load the AWS configuration
    aws_config = aws_utils.AWSConfig(args.aws_config_file)

    # Load the field data products for this quadrant
    xmatch = crossmatch.CrossMatchTable()
    xmatch.load(args.crossmatch_file, log=None)
    quad_phot = hd5_utils.read_phot_from_hd5_file(args.phot_hdf_file,
    										  return_type='array')

    # Select stars in the quadrant
    select = np.where(xmatch.field_index['quadrant'] == args.qid)[0]
    for j,field_id in enumerate(select):
        field_idx = field_id - 1
        quad_id = xmatch.field_index['quadrant_id'][field_idx]
        quad_idx = quad_id - 1

        # Extract the photometry for this star, in a dictionary of arrays
        # of the lightcurve data for each filter
        lc = fetch_photometry_for_star(field_idx, quad_idx, xmatch, quad_phot)

        # Create the IPAC-format multi-extention FITS lightcurve table
        lc_file_path = output_to_ipac_lightcurve(args, field_id, quad_id, xmatch, lc)

        # Upload the lightcurve to the AWS Bucket
        aws_utils.upload_file(aws_config, lc_file_path, args.output_dir,
                              args.aws_root)

        # Remove the temporary local copy of the lightcurve to save space:
        remove(lc_file_path)

        if j%1000 == 0.0:
            print('-> Completed output of lightcurve '+str(j))

def output_to_ipac_lightcurve(args, field_id, quad_id, xmatch, lc):
    """Function to create a lightcurve file for one star in multi-extension
    FITS format, with the lightcurves of the star in different filters contained
    as separate binary tables."""

    quad_idx = quad_id - 1
    field_idx = field_id - 1

    # Primary file header will contain basic identification information
    # for the star, plus information on the data available from the
    # lightcurve tables
    hdr = fits.Header()
    hdr['NAME'] = args.field_name+'_'+str(field_id)
    hdr['FIELD'] = args.field_name
    hdr['FIELD_ID'] = field_id
    hdr['QUADRANT'] = args.qid
    hdr['QUAD_ID'] = quad_id
    hdr['RA'] = xmatch.stars['ra'][field_idx]
    hdr['DEC'] = xmatch.stars['dec'][field_idx]
    hdr['GAIA_ID'] = xmatch.stars['gaia_source_id'][field_idx]
    hdr['GAIACAT'] = 'Gaia_EDR3'
    for f in FILTER_LIST:
        hdr['NDATA_'+str(f).replace('p','').upper()] = len(lc[f])

    # Add the lightcurves in each filter as a binary table extention.
    # This will create zero-length table if no data is available for a given filter.
    hdu_list = [fits.PrimaryHDU(header=hdr)]
    for f in FILTER_LIST:
        fname = str(f).replace('p','').upper()
        hdu_list.append(fits.BinTableHDU(lc[f], name='LIGHTCURVE_'+fname))

    hdu_list = fits.HDUList(hdu_list)

    # Output to disk:
    file_path = crossmatch_to_ipactable.get_lc_file_path(args, field_id)
    file_path = file_path.replace('./','')
    file_path = path.join(args.output_dir,file_path)
    if not path.isdir(path.dirname(file_path)):
        makedirs(path.dirname(file_path))
    hdu_list.writeto(file_path, overwrite=True)

    return file_path

def fetch_photometry_for_star(field_idx, quad_idx, xmatch, quad_phot):

    lc = {}
    (mag_col1, merr_col1) = field_photometry.get_field_photometry_columns('instrumental')
    (mag_col2, merr_col2) = field_photometry.get_field_photometry_columns('calibrated')
    (mag_col3, merr_col3) = field_photometry.get_field_photometry_columns('normalized')
    qc_col = 16

    lc = {}
    for f in FILTER_LIST:
        lc[f] = {'photometry': [], 'dataids': []}

    for f in FILTER_LIST:
        datalist = np.where(xmatch.datasets['dataset_filter'] == f)[0]

        for ddx in datalist:
            dataset = xmatch.datasets[ddx]

            # Extract the photometry of this object for the images from this dataset,
            # if the field index indicates that the object was measured in this dataset
            if xmatch.field_index[dataset['dataset_code']+'_index'][field_idx] != 0:
                shortcode = xmatch.get_dataset_shortcode(dataset['dataset_code'])
                # Select those images from the HDF5 pertaining to this dataset,
                # then select valid measurements for this star
                idx1 = np.where(xmatch.images['dataset_code'] == dataset['dataset_code'])[0]
                idx2 = np.where(quad_phot[quad_idx,:,0] > 0.0)[0]
                idx3 = np.where(quad_phot[quad_idx,:,mag_col1] > 0.0)[0]
                idx = set(idx1).intersection(set(idx2))
                idx = list(idx.intersection(set(idx3)))

                # Store the photometry
                if len(idx) > 0:
                    data = lc[f]

                    # Extract the photometry into holding arrays, handling
                    # string-based data separately from numerical
                    lc_data = np.zeros((len(idx),8))
                    lc_data[:,0] = quad_phot[quad_idx,idx,0]
                    lc_data[:,1] = quad_phot[quad_idx,idx,mag_col1]
                    lc_data[:,2] = quad_phot[quad_idx,idx,merr_col1]
                    lc_data[:,3] = quad_phot[quad_idx,idx,mag_col2]
                    lc_data[:,4] = quad_phot[quad_idx,idx,merr_col2]
                    lc_data[:,5] = quad_phot[quad_idx,idx,mag_col3]
                    lc_data[:,6] = quad_phot[quad_idx,idx,merr_col3]
                    lc_data[:,7] = quad_phot[quad_idx,idx,qc_col]

                    lc_dataids = np.array([dataset['dataset_code']]*len(idx))

                    # Concatenate this dataset's photometry with others using
                    # the same filter:
                    if len(data['photometry']) == 0:
                        data['photometry'] = lc_data
                        data['dataids'] = lc_dataids
                    else:
                        data['photometry'] = np.concatenate((data['photometry'],lc_data))
                        data['dataids'] = np.concatenate((data['dataids'], lc_dataids))

                    lc[f] = data

        # Having accumulated all of the data from different datasets, now
        # convert the arrays into a single astropy Table for easier handling
        data = lc[f]
        if len(data['photometry']) > 0:
            lc[f] = Table([
            Column(name='HJD', data=data['photometry'][:,0]),
            Column(name='inst_mag', data=data['photometry'][:,1]),
            Column(name='inst_mag_error', data=data['photometry'][:,2]),
            Column(name='calib_mag', data=data['photometry'][:,3]),
            Column(name='calib_mag_error', data=data['photometry'][:,4]),
            Column(name='norm_mag', data=data['photometry'][:,5]),
            Column(name='norm_mag_error', data=data['photometry'][:,6]),
            Column(name='qc_flag', data=data['photometry'][:,7]),
            Column(name='dataset', data=data['dataids']),
            ])
        else:
            lc[f] = Table([
            Column(name='HJD', data=np.array([])),
            Column(name='inst_mag', data=np.array([])),
            Column(name='inst_mag_error', data=np.array([])),
            Column(name='calib_mag', data=np.array([])),
            Column(name='calib_mag_error', data=np.array([])),
            Column(name='norm_mag', data=np.array([])),
            Column(name='norm_mag_error', data=np.array([])),
            Column(name='qc_flag', data=np.array([])),
            Column(name='dataset', data=np.array([])),
            ])

    return lc


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('crossmatch_file', type=str,
                    help='Path to crossmatch file')
    parser.add_argument('phot_hdf_file', type=str,
                    help='Path to quadrant HDF file')
    parser.add_argument('qid', type=int,
                    help='Quadrant number to process')
    parser.add_argument('output_dir', type=str,
                    help='Path to output directory')
    parser.add_argument('field_name', type=str,
                    help='Name of the field')
    parser.add_argument('aws_config_file', type=str,
                    help='Path to AWS configuration file')
    parser.add_argument('aws_root', type=str,
                    help='Path to top-level AWS output directory')

    args = parser.parse_args()

    if 'quad'+str(args.qid) not in str(args.phot_hdf_file):
        raise IOError('Requested quadrant ID ('+str(args.qid)
                +') does not match the HDF file given ('+args.phot_hdf_file+')')

    return args


if __name__ == '__main__':
    args = get_args()
    convert_to_ipac_lightcurves(args)
