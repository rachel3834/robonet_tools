from os import path, makedirs, remove
import numpy as np
from pyDANDIA import hd5_utils
from pyDANDIA import crossmatch
from pyDANDIA import field_photometry
from pyDANDIA import logs
import crossmatch_to_ipactable
import aws_utils
import argparse
import json
from astropy.table import Table, Column
from astropy.io import fits

VERSION = 0.1
FILTER_LIST = ['gp', 'rp', 'ip']
upload_aws = False

def convert_to_ipac_lightcurves(args):

    log = logs.start_stage_log(args.output_dir, 'fieldphot_to_ipaclc_Q'+str(args.qid))

    # Load the AWS configuration
    if upload_aws:
        aws_config = aws_utils.AWSConfig(args.aws_config_file)

    # Load the field data products for this quadrant
    xmatch = crossmatch.CrossMatchTable()
    xmatch.load(args.crossmatch_file, log=None)
    quad_phot = hd5_utils.read_phot_from_hd5_file(args.phot_hdf_file,
    										  return_type='array')
    log.info('Loaded field photometry data products for quadrant ' + str(args.qid))

    # Holding dictionary for the count of datapoints per star
    datacounts = {}

    # Select stars in the quadrant
    log.info('Converting lightcurve to IPAC format: ')
    MAXSTAR = 100
    quadrant_rows = np.where(xmatch.field_index['quadrant'] == args.qid)[0]
    for j,row in enumerate(quadrant_rows):
        field_id = xmatch.field_index['field_id'][row]
        field_idx = field_id - 1
        quad_id = xmatch.field_index['quadrant_id'][row]
        quad_idx = quad_id - 1

        # Extract the photometry for this star, in a dictionary of arrays
        # of the lightcurve data for each filter
        lc, star_datacount = fetch_photometry_for_star(quad_idx, xmatch, quad_phot)
        datacounts[int(field_id)] = star_datacount

        # Create the IPAC-format multi-extention FITS lightcurve table
        if star_datacount['gp'] > 0 or star_datacount['rp'] > 0 or star_datacount['ip'] > 0:
            lc_file_path = output_to_ipac_lightcurve(args, field_id, quad_id, xmatch, lc)

        # Upload the lightcurve to the AWS Bucket
        if upload_aws:
            aws_utils.upload_file(aws_config, lc_file_path, args.output_dir,
                              args.aws_root)

        # Remove the temporary local copy of the lightcurve to save space:
        #remove(lc_file_path)

        if j%1000 == 0.0:
            log.info('-> Completed output of lightcurve ' + str(j) + ', row='
                     + str(row) + ', field_id=' + str(field_id) + ', quad_id=' + str(quad_id))

        # TEMPORARY CAP FOR TESTING
        #if j%10 == 0.0:
        #    log.info('-> Completed output of lightcurve '+str(j))
        #if j==MAXSTAR:
        #    break

    log.info('Finished converting lightcurves')

    # Output summary of the number of datapoints per star:
    json_data = json.dumps(datacounts, indent=4)
    file_path = path.join(args.output_dir, args.field_name + '_starcounts_Q'+str(args.qid)+'.json')
    with open(file_path, 'w') as write_file:
        write_file.write(json_data)
        write_file.close()
    log.info('Output star datapoint counts')

    log.info('End of processing')
    logs.close_log(log)

def output_to_ipac_lightcurve(args, field_id, quad_id, xmatch, lc):
    """Function to create a lightcurve file for one star in multi-extension
    FITS format, with the lightcurves of the star in different filters contained
    as separate binary tables."""

    field_idx = field_id - 1

    # Primary file header will contain basic identification information
    # for the star, plus information on the data available from the
    # lightcurve tables
    hdr = fits.Header()
    hdr['NAME'] = args.field_name+'_'+str(crossmatch_to_ipactable.bin_num_string(field_id))
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
        hdu_list.append(fits.BinTableHDU(lc[f], name='LIGHTCURVE_SDSS_'+fname))

    hdu_list = fits.HDUList(hdu_list)

    # Output to disk:
    file_path = crossmatch_to_ipactable.get_lc_file_path(args, field_id)
    file_path = file_path.replace('./','')
    file_path = path.join(args.output_dir,file_path)
    if not path.isdir(path.dirname(file_path)):
        makedirs(path.dirname(file_path))
    hdu_list.writeto(file_path, overwrite=True)

    return file_path

def fetch_photometry_for_star(quad_idx, xmatch, quad_phot):

    (mag_col1, merr_col1) = field_photometry.get_field_photometry_columns('instrumental')
    (mag_col2, merr_col2) = field_photometry.get_field_photometry_columns('calibrated')
    (mag_col3, merr_col3) = field_photometry.get_field_photometry_columns('normalized')
    qc_col = 16

    lc = {}
    datacount = {}
    for f in FILTER_LIST:
        # Identify valid photometry array entries for all images of the current star in the current filter
        # The second index is designed to eliminate datapoints with HJD=0.0, which is a flag that the
        # DATE-OBS keyword from the original image was mis-formed due to a temporary issue at the telescope,
        # so no accurate HJD can be calculated.
        idx1 = np.where(xmatch.images['filter'] == f)[0]
        idx2 = np.where(quad_phot[quad_idx, :, 0] > 0.0)[0]
        idx3 = np.where(quad_phot[quad_idx, :, mag_col1] > 0.0)[0]
        idx = set(idx1).intersection(set(idx2))
        idx = list(idx.intersection(set(idx3)))

        lc_data = []
        if len(idx) > 0:
            # Convert the dataset codes for these images from long to short form
            dataset_codes = [xmatch.get_dataset_shortcode(x) for x in xmatch.images['dataset_code'][idx]]

            # Compile the lightcurve table data for this filter
            lc_data.append(Column(name='HJD', data=quad_phot[quad_idx, idx, 0]))
            lc_data.append(Column(name='inst_mag', data=quad_phot[quad_idx, idx, mag_col1]))
            lc_data.append(Column(name='inst_mag_error', data=quad_phot[quad_idx, idx, merr_col1]))
            lc_data.append(Column(name='calib_mag', data=quad_phot[quad_idx, idx, mag_col2]))
            lc_data.append(Column(name='calib_mag_error', data=quad_phot[quad_idx, idx, merr_col2]))
            lc_data.append(Column(name='norm_mag', data=quad_phot[quad_idx, idx, mag_col3]))
            lc_data.append(Column(name='norm_mag_error', data=quad_phot[quad_idx, idx, merr_col3]))
            lc_data.append(Column(name='qc_flag', data=quad_phot[quad_idx, idx, qc_col]))
            lc_data.append(Column(name='dataset', data=dataset_codes))
            lc_data.append(Column(name='airmass', data=xmatch.images['airmass'][idx]))
            lc_data.append(Column(name='moon_frac', data=xmatch.images['moon_fraction'][idx]))
            lc_data.append(Column(name='moon_sep', data=xmatch.images['moon_ang_separation'][idx]))
            lc_data.append(Column(name='sky_bkgd', data=xmatch.images['sky'][idx]))
            lc_data.append(Column(name='fwhm', data=xmatch.images['fwhm'][idx]))

        # If there are no valid photometry measurements in this filter for this star,
        # return zero-length table in the same format
        else:
            lc_data.append(Column(name='HJD', data=np.array([])))
            lc_data.append(Column(name='inst_mag', data=np.array([])))
            lc_data.append(Column(name='inst_mag_error', data=np.array([])))
            lc_data.append(Column(name='calib_mag', data=np.array([])))
            lc_data.append(Column(name='calib_mag_error', data=np.array([])))
            lc_data.append(Column(name='norm_mag', data=np.array([])))
            lc_data.append(Column(name='norm_mag_error', data=np.array([])))
            lc_data.append(Column(name='qc_flag', data=np.array([])))
            lc_data.append(Column(name='dataset', data=np.array([])))
            lc_data.append(Column(name='airmass', data=np.array([])))
            lc_data.append(Column(name='moon_frac', data=np.array([])))
            lc_data.append(Column(name='moon_sep', data=np.array([])))
            lc_data.append(Column(name='sky_bkgd', data=np.array([])))
            lc_data.append(Column(name='fwhm', data=np.array([])))

        lc[f] = Table(lc_data)
        lc[f].sort(['HJD'])
        datacount[f] = len(lc[f])

    return lc, datacount


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
