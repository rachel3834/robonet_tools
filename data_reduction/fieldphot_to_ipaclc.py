from os import path
import numpy as np
from pyDANDIA import hd5_utils
from pyDANDIA import crossmatch
from pyDANDIA import field_photometry
from pyDANDIA import crossmatch_to_ipactable
import argparse
from datetime import datetime
from astropy.table import Table, Column

VERSION = 0.1

def convert_to_ipac_lightcurves(args):

    # Load the field data products for this quadrant
    xmatch = crossmatch.CrossMatchTable()
	xmatch.load(args.crossmatch_file, log=log)
    quad_phot = hd5_utils.read_phot_from_hd5_file(args.phot_hdf_file,
												  return_type='array')

    # Select stars in the quadrant
    select = np.where(xmatch.field_index['quadrant'] == args.qid)[0]
    for field_id in select:
        field_idx = field_id - 1
        quad_id = xmatch.field_index['quadrant_id'][field_idx]

        # Extract the photometry for this star, in a dictionary of arrays
        # of the lightcurve data for each filter
        lc = fetch_photometry_for_star(quad_idx, xmatch, quad_phot)

        # Create the IPAC-format multi-extention FITS lightcurve table

def output_to_ipac_lightcurve(args, field_id, quad_id, lc, filter_list):
    """Function to create a lightcurve file for one star in multi-extension
    FITS format, with the lightcurves of the star in different filters contained
    as separate binary tables."""

    # Primary file header will contain basic identification information
    # for the star, plus information on the data available from the
    # lightcurve tables
    hdr = fits.Header()
    hdr['NAME'] = args.field_name+'_'+str(field_id)
    hdr['FIELD_ID'] = args.field_name
    hdr['QUAD_ID'] = quad_id
    hdr['RA'] = ra
    hdr['DEC'] = dec
    hdr['GAIA_ID'] = gaia_id
    hdr['GAIACAT'] = 'Gaia_EDR3'
    for f in filter_list:
        hdr['NDATA_'+str(f).replace('p','').upper()] = len(lc[f])

    # Add the lightcurves in each filter as a binary table extention.
    # This will create zero-length table if no data is available for a given filter.
    hdu_list = [fits.PrimaryHDU(header=hdr)]
    for f in filter_list:
        fname = str(f).replace('p','').upper()
        hdu_list.append(fits.BinTableHDU(lc[f], name='LIGHTCURVE_'+fname))

    hdu_list = fits.HDUList(hdu_list)

    # Output to disk:
    file_path = path.join(args.output_dir, '')
    hdu_list.writeto(file_path, overwrite=True)

def fetch_photometry_for_star(quad_idx, xmatch, quad_phot):

	lc = {}
	(mag_col1, merr_col1) = field_photometry.get_field_photometry_columns('instrumental')
	(mag_col2, merr_col2) = field_photometry.get_field_photometry_columns('calibrated')
	(mag_col3, merr_col3) = field_photometry.get_field_photometry_columns('normalized')
	qc_col = 16

    # Initialize holding lists for lightcurve data
    filter_list = ['gp', 'rp', 'ip']
    lc = {}
    for f in filter_list:
        lc[f] = []

    for f in filter_list:
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
    			idx3 = np.where(quad_phot[quad_idx,:,mag_col] > 0.0)[0]
    			idx = set(idx1).intersection(set(idx2))
    			idx = list(idx.intersection(set(idx3)))

    			# Store the photometry
    			if len(idx) > 0:
                    data = lc[f]

                    lc_data = quad_phto

                    data.append([
                            quad_phot[quad_idx,idx,0],          # HJD
                            quad_phot[quad_idx,idx,mag_col1],   # Instrumental mag
                            quad_phot[quad_idx,idx,merr_col1],
                            quad_phot[quad_idx,idx,mag_col2],   # Calibrated mag
                            quad_phot[quad_idx,idx,merr_col2],
                            quad_phot[quad_idx,idx,mag_col3],   # Normalized mag
                            quad_phot[quad_idx,idx,merr_col3],
                            quad_phot[quad_idx,idx,qc_col],     # QC
                            [dataset['dataset_code']]*len(idx)             # Dataset ID
                            ])

                    lc[f] = data
    	lc[f] = np.array(lc[f])

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

    args = parser.parse_args()

    if 'quad_'+str(qid) not in args.phot_hdf_file:
        raise IOError('Requested quadrant ID ('+str(qid)
                +') does not match the HDF file given ('+args.phot_hdf_file+')')

    return args


if __name__ == '__main__':
    args = get_args()
    convert_to_ipac_lightcurves(args)
